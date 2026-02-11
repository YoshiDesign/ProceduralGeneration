package duals2

import (
	"image/color"
	"math"
	"procedural_generation/terrain_generation/duals2/core"
	"procedural_generation/terrain_generation/duals2/hydro"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

// DualsDemo is a visualization module for the Delaunay terrain system.
type DualsDemo struct {
	Manager       *ChunkManager   // Chunk manager with caching
	Chunks        []*core.TerrainChunk // Active chunks for rendering
	ScreenW       int
	ScreenH       int
	Scale         float64 // pixels per world unit
	OffsetX       float64
	OffsetZ       float64
	ShowVoronoi   bool
	ShowNormals   bool
	ShowHydrology bool // Toggle hydrology visualization

	// Debug UI for parameter tuning
	DebugUI   *DebugUI
	chunkCfg  core.ChunkConfig // Store for regeneration
}

func exampleHeightFunc(x, z float64) float64 {
	// Multi-octave sine waves for terrain-like appearance
	h := 0.0
	h += 15.0 * math.Sin(x*0.02) * math.Cos(z*0.02)
	h += 8.0 * math.Sin(x*0.05+1.0) * math.Cos(z*0.04+0.5)
	h += 4.0 * math.Sin(x*0.1+2.0) * math.Sin(z*0.08+1.0)
	return h
}

func FractalNoise(x, z float64, octaves int, persistence, lacunarity float64) float64 {
    total := 0.0
    amplitude := 1.0
    frequency := 1.0
    maxValue := 0.0
    
    for i := 0; i < octaves; i++ {
        total += Simplex2D(x*frequency, z*frequency) * amplitude
        maxValue += amplitude
        amplitude *= persistence   // Each octave is quieter
        frequency *= lacunarity    // Each octave is higher frequency
    }
    
    return total / maxValue  // Normalize to [-1, 1]
}

// FractalNoiseV2 computes fractal noise using simplex noise.
// Good defaults: amplitude 1.0, persistence 0.5, lacunarity 2.0
func FractalNoiseV2(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64 {
	sum := 0.0
	for i := 0; i < octaves; i++ {
		sum += Simplex2D(x*frequency, z*frequency) * amplitude
		amplitude *= persistence
		frequency *= lacunarity
	}
	return sum
}

// NewDualsDemo creates a demo with a 2x2 grid of chunks using the ChunkManager.
func NewDualsDemo(screenW, screenH int) *DualsDemo {
	cfg := DefaultChunkConfig()

	// Height function: simple noise-like terrain
	heightFunc := FractalNoiseV2

	// Create chunk manager with caching
	manager := NewChunkManager(cfg, heightFunc)

	// Generate a X,Z grid of chunks using the manager (benefits from caching)
	chunks := make([]*core.TerrainChunk, 0, cfg.ChunksZ * cfg.ChunksX)
	for cz := 0; cz < cfg.ChunksZ; cz++ {
		for cx := 0; cx < cfg.ChunksX; cx++ {
			coord := core.ChunkCoord{X: cx, Z: cz}
			chunk, err := manager.GetOrGenerate(coord)
			if err != nil {
				continue
			}
			chunks = append(chunks, chunk)
		}
	}

	// Scale to fit screen
	totalWorldSize := cfg.ChunkSize * 4
	scale := float64(min(screenW, screenH)) / totalWorldSize * 0.9

	demo := &DualsDemo{
		Manager:       manager,
		Chunks:        chunks,
		ScreenW:       screenW,
		ScreenH:       screenH,
		Scale:         scale,
		OffsetX:       float64(screenW) * 0.05,
		OffsetZ:       float64(screenH) * 0.05,
		ShowVoronoi:   false,
		ShowNormals:   false,
		ShowHydrology: true, // Enable hydrology visualization by default
		chunkCfg:      cfg, // Do not use for anything other than regeneration. Use the manager's config instead.
	}

	// Create debug UI with regeneration callback
	demo.DebugUI = NewDebugUI(
		manager.NoiseParams(),
		manager.HydroConfig(),
		demo.regenerate,
	)

	// Pre-compute render data for all chunks
	demo.buildAllRenderData()

	return demo
}

// worldToScreen converts world coordinates to screen coordinates.
func (d *DualsDemo) worldToScreen(x, z float64) (float32, float32) {
	sx := d.OffsetX + x*d.Scale
	sz := d.OffsetZ + z*d.Scale
	return float32(sx), float32(sz)
}

// buildRenderData pre-computes vertices and indices for batched triangle rendering.
// This should be called once when chunks are generated or when view parameters change.
func (d *DualsDemo) buildRenderData(chunk *core.TerrainChunk) {
	mesh := chunk.Mesh
	if mesh == nil {
		chunk.RenderVertices = nil
		chunk.RenderIndices = nil
		return
	}

	// Build a set of core site indices for fast lookup
	coreSet := make(map[core.SiteIndex]struct{}, len(chunk.CoreSiteIndices))
	for _, ci := range chunk.CoreSiteIndices {
		coreSet[ci] = struct{}{}
	}

	// Pre-compute the light direction (same as in drawChunk)
	light := core.Vec3{X: -0.5, Y: 0.8, Z: -0.3}.Normalize()

	// Count core triangles first to pre-allocate
	coreTriCount := 0
	for _, t := range mesh.Tris {
		_, aCore := coreSet[core.SiteIndex(t.A)]
		_, bCore := coreSet[core.SiteIndex(t.B)]
		_, cCore := coreSet[core.SiteIndex(t.C)]
		if aCore || bCore || cCore {
			coreTriCount++
		}
	}

	// Pre-allocate slices (3 vertices and 3 indices per triangle)
	vertices := make([]ebiten.Vertex, 0, coreTriCount*3)
	indices := make([]uint16, 0, coreTriCount*3)

	vertexIdx := uint16(0)

	for ti, t := range mesh.Tris {
		// Only include triangles with at least one core vertex
		_, aCore := coreSet[core.SiteIndex(t.A)]
		_, bCore := coreSet[core.SiteIndex(t.B)]
		_, cCore := coreSet[core.SiteIndex(t.C)]
		if !aCore && !bCore && !cCore {
			continue
		}

		// Get world positions
		a := mesh.Sites[t.A].Pos
		b := mesh.Sites[t.B].Pos
		c := mesh.Sites[t.C].Pos

		// Convert to screen coordinates
		ax, ay := d.worldToScreen(a.X, a.Y)
		bx, by := d.worldToScreen(b.X, b.Y)
		cx, cy := d.worldToScreen(c.X, c.Y)

		// Compute color (same logic as drawChunk)
		var r, g, b_col, alpha float32
		if ti < len(chunk.Mesh.FaceNormals) {
			n := chunk.Mesh.FaceNormals[ti]
			intensity := n.Dot(light)
			if intensity < 0 {
				intensity = 0
			}
			intensity = 0.3 + 0.7*intensity // Ambient + diffuse

			// Base color varies with average height
			avgH := (chunk.Heights[t.A] + chunk.Heights[t.B] + chunk.Heights[t.C])
			baseR := 80.0 + avgH*2
			baseG := 120.0 + avgH*1.5
			baseB := 80.0 + avgH*0.5

			r = float32(clamp(baseR*intensity, 0, 255)) / 255.0
			g = float32(clamp(baseG*intensity, 0, 255)) / 255.0
			b_col = float32(clamp(baseB*intensity, 0, 255)) / 255.0
			alpha = 1.0
		} else {
			r = 100.0 / 255.0
			g = 140.0 / 255.0
			b_col = 100.0 / 255.0
			alpha = 1.0
		}

		// Add three vertices for this triangle
		vertices = append(vertices,
			ebiten.Vertex{DstX: ax, DstY: ay, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b_col, ColorA: alpha},
			ebiten.Vertex{DstX: bx, DstY: by, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b_col, ColorA: alpha},
			ebiten.Vertex{DstX: cx, DstY: cy, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b_col, ColorA: alpha},
		)

		// Add indices for this triangle
		indices = append(indices, vertexIdx, vertexIdx+1, vertexIdx+2)
		vertexIdx += 3
	}

	chunk.RenderVertices = vertices
	chunk.RenderIndices = indices
}

// buildAllRenderData rebuilds render data for all chunks.
func (d *DualsDemo) buildAllRenderData() {
	for _, chunk := range d.Chunks {
		d.buildRenderData(chunk)
		d.buildHydrologyRenderData(chunk)
	}
}

// buildHydrologyRenderData pre-computes vertices for lake and ocean triangles.
func (d *DualsDemo) buildHydrologyRenderData(chunk *core.TerrainChunk) {
	// Clear existing hydrology render data
	chunk.LakeVertices = nil
	chunk.LakeIndices = nil
	chunk.OceanVertices = nil
	chunk.OceanIndices = nil

	if chunk.Hydro == nil || chunk.Mesh == nil {
		return
	}

	// Build lake render data
	d.buildLakeRenderData(chunk)

	// Build ocean render data
	d.buildOceanRenderData(chunk)
}

// buildLakeRenderData pre-computes vertices for all lake triangles in a chunk.
func (d *DualsDemo) buildLakeRenderData(chunk *core.TerrainChunk) {
	if len(chunk.Hydro.Lakes) == 0 {
		return
	}

	// Lake color (same as in drawLakes)
	r := float32(40) / 255.0
	g := float32(140) / 255.0
	b := float32(160) / 255.0
	a := float32(180) / 255.0

	// Collect all lake site indices across all lakes
	allLakeSites := make(map[core.SiteIndex]struct{})
	for _, lake := range chunk.Hydro.Lakes {
		for _, idx := range lake.SiteIndices {
			allLakeSites[idx] = struct{}{}
		}
	}

	// Count lake triangles for pre-allocation
	lakeTriCount := 0
	for _, tri := range chunk.Mesh.Tris {
		_, aLake := allLakeSites[tri.A]
		_, bLake := allLakeSites[tri.B]
		_, cLake := allLakeSites[tri.C]
		if aLake && bLake && cLake {
			lakeTriCount++
		}
	}

	if lakeTriCount == 0 {
		return
	}

	vertices := make([]ebiten.Vertex, 0, lakeTriCount*3)
	indices := make([]uint16, 0, lakeTriCount*3)
	vertexIdx := uint16(0)

	for _, tri := range chunk.Mesh.Tris {
		_, aLake := allLakeSites[tri.A]
		_, bLake := allLakeSites[tri.B]
		_, cLake := allLakeSites[tri.C]
		if !aLake || !bLake || !cLake {
			continue
		}

		posA := chunk.Mesh.Sites[tri.A].Pos
		posB := chunk.Mesh.Sites[tri.B].Pos
		posC := chunk.Mesh.Sites[tri.C].Pos

		ax, ay := d.worldToScreen(posA.X, posA.Y)
		bx, by := d.worldToScreen(posB.X, posB.Y)
		cx, cy := d.worldToScreen(posC.X, posC.Y)

		vertices = append(vertices,
			ebiten.Vertex{DstX: ax, DstY: ay, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
			ebiten.Vertex{DstX: bx, DstY: by, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
			ebiten.Vertex{DstX: cx, DstY: cy, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
		)
		indices = append(indices, vertexIdx, vertexIdx+1, vertexIdx+2)
		vertexIdx += 3
	}

	chunk.LakeVertices = vertices
	chunk.LakeIndices = indices
}

// buildOceanRenderData pre-computes vertices for ocean triangles with depth-based coloring.
func (d *DualsDemo) buildOceanRenderData(chunk *core.TerrainChunk) {
	ocean := chunk.Hydro.Ocean
	if !ocean.IsOcean || len(ocean.SubmergedTris) == 0 {
		return
	}

	// Get sea level for depth calculation
	oceanRegion := d.Manager.HydroManager().GetOceanForChunk(chunk.Coord)
	seaLevel := float64(0)
	if oceanRegion != nil {
		seaLevel = oceanRegion.SeaLevel
	}

	vertices := make([]ebiten.Vertex, 0, len(ocean.SubmergedTris)*3)
	indices := make([]uint16, 0, len(ocean.SubmergedTris)*3)
	vertexIdx := uint16(0)

	for _, triID := range ocean.SubmergedTris {
		if triID < 0 || triID >= len(chunk.Mesh.Tris) {
			continue
		}

		tri := chunk.Mesh.Tris[triID]
		posA := chunk.Mesh.Sites[tri.A].Pos
		posB := chunk.Mesh.Sites[tri.B].Pos
		posC := chunk.Mesh.Sites[tri.C].Pos

		ax, ay := d.worldToScreen(posA.X, posA.Y)
		bx, by := d.worldToScreen(posB.X, posB.Y)
		cx, cy := d.worldToScreen(posC.X, posC.Y)

		// Color based on depth (same logic as drawOcean)
		avgHeight := (chunk.Heights[tri.A] + chunk.Heights[tri.B] + chunk.Heights[tri.C]) / 3.0
		depth := seaLevel - avgHeight
		depthFactor := clamp(depth/20.0, 0, 1)

		r := float32(30-20*depthFactor) / 255.0
		g := float32(80-30*depthFactor) / 255.0
		b := float32(140+40*depthFactor) / 255.0
		a := float32(180+50*depthFactor) / 255.0

		vertices = append(vertices,
			ebiten.Vertex{DstX: ax, DstY: ay, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
			ebiten.Vertex{DstX: bx, DstY: by, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
			ebiten.Vertex{DstX: cx, DstY: cy, SrcX: 1, SrcY: 1, ColorR: r, ColorG: g, ColorB: b, ColorA: a},
		)
		indices = append(indices, vertexIdx, vertexIdx+1, vertexIdx+2)
		vertexIdx += 3
	}

	chunk.OceanVertices = vertices
	chunk.OceanIndices = indices
}

// Update handles input and updates UI state.
func (d *DualsDemo) Update() {
	// Toggle debug UI with D key
	if inpututil.IsKeyJustPressed(ebiten.KeyD) {
		d.DebugUI.Toggle()
	}

	// Update debug UI
	d.DebugUI.Update()
}

// regenerate rebuilds all chunks with new parameters.
func (d *DualsDemo) regenerate(noiseParams core.NoiseParams, hydroConfig hydro.HydroConfig) {

	cfg := DefaultChunkConfig()

	// Update manager configurations
	d.Manager.SetNoiseParams(noiseParams)
	d.Manager.SetHydroConfig(hydroConfig)
	d.Manager.ClearCaches()

	// Regenerate all chunks
	d.Chunks = make([]*core.TerrainChunk, 0, cfg.ChunksZ * cfg.ChunksX)
	for cz := 0; cz < d.Manager.cfg.ChunksZ; cz++ {
		for cx := 0; cx < d.Manager.cfg.ChunksX; cx++ {
			coord := core.ChunkCoord{X: cx, Z: cz}
			chunk, err := d.Manager.GetOrGenerate(coord)
			if err != nil {
				continue
			}
			d.Chunks = append(d.Chunks, chunk)
		}
	}

	// Rebuild render data for all chunks
	d.buildAllRenderData()
}

// Draw renders the terrain chunks.
func (d *DualsDemo) Draw(screen *ebiten.Image) {
	// Clear with dark background
	screen.Fill(color.RGBA{20, 25, 30, 255})

	for _, chunk := range d.Chunks {
		d.drawChunk(screen, chunk)
	}

	// Draw hydrology on top of terrain
	if d.ShowHydrology {
		for _, chunk := range d.Chunks {
			d.drawChunkHydrology(screen, chunk)
		}
	}

	// Draw chunk boundaries
	for _, chunk := range d.Chunks {
		d.drawChunkBoundary(screen, chunk)
	}

	// Draw debug UI on top
	d.DebugUI.Draw(screen)
}

// drawChunkHydrology renders all water features for a chunk.
func (d *DualsDemo) drawChunkHydrology(screen *ebiten.Image, chunk *core.TerrainChunk) {
	if chunk.Hydro == nil {
		return
	}

	// Draw in order: oceans (bottom), lakes (middle), rivers (top)
	d.drawOcean(screen, chunk)
	d.drawLakes(screen, chunk)
	d.drawRivers(screen, chunk)
}

// drawChunk renders a single terrain chunk using batched DrawTriangles.
func (d *DualsDemo) drawChunk(screen *ebiten.Image, chunk *core.TerrainChunk) {
	if len(chunk.RenderVertices) == 0 || len(chunk.RenderIndices) == 0 {
		return
	}

	// Single batched draw call for all triangles in this chunk
	screen.DrawTriangles(chunk.RenderVertices, chunk.RenderIndices, emptyImage(), nil)
}

// drawChunkBoundary draws the chunk boundary.
func (d *DualsDemo) drawChunkBoundary(screen *ebiten.Image, chunk *core.TerrainChunk) {
	boundaryColor := color.RGBA{100, 100, 200, 200}

	x1, y1 := d.worldToScreen(chunk.MinX, chunk.MinZ)
	x2, y2 := d.worldToScreen(chunk.MaxX, chunk.MinZ)
	x3, y3 := d.worldToScreen(chunk.MaxX, chunk.MaxZ)
	x4, y4 := d.worldToScreen(chunk.MinX, chunk.MaxZ)

	vector.StrokeLine(screen, x1, y1, x2, y2, 2, boundaryColor, false)
	vector.StrokeLine(screen, x2, y2, x3, y3, 2, boundaryColor, false)
	vector.StrokeLine(screen, x3, y3, x4, y4, 2, boundaryColor, false)
	vector.StrokeLine(screen, x4, y4, x1, y1, 2, boundaryColor, false)
}

// emptyImage returns a 1x1 white image for solid color fills.
var emptyImageInstance *ebiten.Image

func emptyImage() *ebiten.Image {
	if emptyImageInstance == nil {
		emptyImageInstance = ebiten.NewImage(3, 3)
		emptyImageInstance.Fill(color.White)
	}
	return emptyImageInstance
}

func clamp(v, minV, maxV float64) float64 {
	if v < minV {
		return minV
	}
	if v > maxV {
		return maxV
	}
	return v
}

// -------------------------------------------------------------------
// Hydrology Rendering
// -------------------------------------------------------------------

// drawRivers renders river segments as width-varying blue lines.
func (d *DualsDemo) drawRivers(screen *ebiten.Image, chunk *core.TerrainChunk) {
	if chunk.Hydro == nil {
		return
	}

	for _, segment := range chunk.Hydro.Rivers {
		if len(segment.Vertices) < 2 { // TODO: at this point, this shouldn't be possible
			continue
		}

		// Draw river segments with varying width
		for i := 0; i < len(segment.Vertices)-1; i++ {
			p1 := segment.Vertices[i]
			p2 := segment.Vertices[i+1]

			x1, y1 := d.worldToScreen(p1.X, p1.Y)
			x2, y2 := d.worldToScreen(p2.X, p2.Y)

			// Width varies along the river (grows downstream)
			width := float32(2.0)
			if i < len(segment.Widths) {
				width = float32(segment.Widths[i] * d.Scale * 0.3)
				if width < 1.5 {
					width = 1.5
				}
				if width > 8 {
					width = 8
				}
			}

			// Color: deeper blue for wider sections
			intensity := float64(width) / 8.0
			riverColor := color.RGBA{
				R: uint8(40 + 30*(1-intensity)),
				G: uint8(100 + 50*(1-intensity)),
				B: uint8(180 + 40*intensity),
				A: 220,
			}

			// Draw the river segment
			vector.StrokeLine(screen, x1, y1, x2, y2, width, riverColor, false)

			// Draw lighter edge highlight for wider rivers
			if width > 3 {
				edgeColor := color.RGBA{100, 160, 220, 150}
				vector.StrokeLine(screen, x1, y1, x2, y2, width+1.5, edgeColor, false)
			}
		}

		// Draw source marker (small circle at start)
		if len(segment.Vertices) > 0 {
			src := segment.Vertices[0]
			sx, sy := d.worldToScreen(src.X, src.Y)
			sourceColor := color.RGBA{255, 255, 255, 200}
			vector.FillCircle(screen, sx, sy, 3, sourceColor, false)
		}
	}
}

// drawLakes renders lakes as filled water areas using batched DrawTriangles.
func (d *DualsDemo) drawLakes(screen *ebiten.Image, chunk *core.TerrainChunk) {
	if chunk.Hydro == nil || chunk.Mesh == nil {
		return
	}

	// Draw batched lake triangles
	if len(chunk.LakeVertices) > 0 && len(chunk.LakeIndices) > 0 {
		screen.DrawTriangles(chunk.LakeVertices, chunk.LakeIndices, emptyImage(), nil)
	}

	// Draw lake edges and spillways (these remain unbatched as they're less performance-critical)
	lakeEdgeColor := color.RGBA{60, 180, 200, 200}

	for _, lake := range chunk.Hydro.Lakes {
		if len(lake.SiteIndices) == 0 {
			continue
		}

		// Create a set of lake site indices for fast lookup
		lakeSites := make(map[core.SiteIndex]bool)
		for _, idx := range lake.SiteIndices {
			lakeSites[idx] = true
		}

		// Draw lake outline by finding boundary edges
		for _, tri := range chunk.Mesh.Tris {
			inCount := 0
			if lakeSites[tri.A] {
				inCount++
			}
			if lakeSites[tri.B] {
				inCount++
			}
			if lakeSites[tri.C] {
				inCount++
			}

			// Triangle on lake boundary (some vertices in, some out)
			if inCount > 0 && inCount < 3 {
				a := chunk.Mesh.Sites[tri.A].Pos
				b := chunk.Mesh.Sites[tri.B].Pos
				c := chunk.Mesh.Sites[tri.C].Pos

				ax, ay := d.worldToScreen(a.X, a.Y)
				bx, by := d.worldToScreen(b.X, b.Y)
				cx, cy := d.worldToScreen(c.X, c.Y)

				// Draw edges that are on the boundary
				if lakeSites[tri.A] != lakeSites[tri.B] {
					vector.StrokeLine(screen, ax, ay, bx, by, 2, lakeEdgeColor, false)
				}
				if lakeSites[tri.B] != lakeSites[tri.C] {
					vector.StrokeLine(screen, bx, by, cx, cy, 2, lakeEdgeColor, false)
				}
				if lakeSites[tri.C] != lakeSites[tri.A] {
					vector.StrokeLine(screen, cx, cy, ax, ay, 2, lakeEdgeColor, false)
				}
			}
		}

		// Draw spillway indicator if present
		if lake.Spillway != nil {
			spx, spy := d.worldToScreen(lake.Spillway.X, lake.Spillway.Y)
			spillwayColor := color.RGBA{100, 220, 255, 255}
			vector.FillCircle(screen, spx, spy, 4, spillwayColor, false)
		}
	}
}

// drawOcean renders ocean areas with submerged triangles and coastline using batched DrawTriangles.
func (d *DualsDemo) drawOcean(screen *ebiten.Image, chunk *core.TerrainChunk) {
	if chunk.Hydro == nil || chunk.Mesh == nil {
		return
	}

	ocean := chunk.Hydro.Ocean
	if !ocean.IsOcean {
		return
	}

	// Draw batched ocean triangles
	if len(chunk.OceanVertices) > 0 && len(chunk.OceanIndices) > 0 {
		screen.DrawTriangles(chunk.OceanVertices, chunk.OceanIndices, emptyImage(), nil)
	}

	// Draw coastline points and lines (these remain unbatched as they're less performance-critical)
	if len(ocean.CoastlinePts) > 0 {
		coastColor := color.RGBA{200, 220, 255, 255}
		for _, pt := range ocean.CoastlinePts {
			px, py := d.worldToScreen(pt.X, pt.Y)
			vector.FillCircle(screen, px, py, 2, coastColor, false)
		}

		// Connect coastline points that are close together
		coastLineColor := color.RGBA{180, 200, 240, 200}
		for i := 0; i < len(ocean.CoastlinePts)-1; i++ {
			p1 := ocean.CoastlinePts[i]
			p2 := ocean.CoastlinePts[i+1]

			// Only connect if reasonably close
			dist := p1.Sub(p2).Len()
			if dist < 20 { // Within reasonable distance
				x1, y1 := d.worldToScreen(p1.X, p1.Y)
				x2, y2 := d.worldToScreen(p2.X, p2.Y)
				vector.StrokeLine(screen, x1, y1, x2, y2, 1.5, coastLineColor, false)
			}
		}
	}
}
