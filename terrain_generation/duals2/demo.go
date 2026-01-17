package duals2

import (
	"image/color"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

// DualsDemo is a visualization module for the Delaunay terrain system.
type DualsDemo struct {
	Manager       *ChunkManager   // Chunk manager with caching
	Chunks        []*TerrainChunk // Active chunks for rendering
	ScreenW       int
	ScreenH       int
	Scale         float64 // pixels per world unit
	OffsetX       float64
	OffsetZ       float64
	ShowVoronoi   bool
	ShowNormals   bool
	ShowHydrology bool // Toggle hydrology visualization
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
	cfg := ChunkConfig{
		ChunkSize:    128.0,
		MinPointDist: 8.0,
		HaloWidth:    12.0,
		WorldSeed:    12345,
	}

	// Height function: simple noise-like terrain
	heightFunc := FractalNoiseV2

	// Create chunk manager with caching
	manager := NewChunkManager(cfg, heightFunc)

	// Generate a 3x2 grid of chunks using the manager (benefits from caching)
	chunks := make([]*TerrainChunk, 0, 6)
	for cz := 0; cz < 4; cz++ {
		for cx := 0; cx < 6; cx++ {
			coord := ChunkCoord{X: cx, Z: cz}
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

	return &DualsDemo{
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
	}
}

// worldToScreen converts world coordinates to screen coordinates.
func (d *DualsDemo) worldToScreen(x, z float64) (float32, float32) {
	sx := d.OffsetX + x*d.Scale
	sz := d.OffsetZ + z*d.Scale
	return float32(sx), float32(sz)
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
}

// drawChunkHydrology renders all water features for a chunk.
func (d *DualsDemo) drawChunkHydrology(screen *ebiten.Image, chunk *TerrainChunk) {
	if chunk.Hydro == nil {
		return
	}

	// Draw in order: oceans (bottom), lakes (middle), rivers (top)
	d.drawOcean(screen, chunk)
	d.drawLakes(screen, chunk)
	d.drawRivers(screen, chunk)
}

// drawChunk renders a single terrain chunk.
func (d *DualsDemo) drawChunk(screen *ebiten.Image, chunk *TerrainChunk) {
	mesh := chunk.Mesh
	if mesh == nil {
		return
	}

	// Draw triangles
	for ti, t := range mesh.Tris {
		// Only draw triangles with at least one core vertex
		isCore := false
		for _, ci := range chunk.CoreSiteIndices {
			if t.A == ci || t.B == ci || t.C == ci {
				isCore = true
				break
			}
		}
		if !isCore {
			continue
		}

		a := mesh.Sites[t.A].Pos
		b := mesh.Sites[t.B].Pos
		c := mesh.Sites[t.C].Pos

		ax, ay := d.worldToScreen(a.X, a.Y)
		bx, by := d.worldToScreen(b.X, b.Y)
		cx, cy := d.worldToScreen(c.X, c.Y)

		// Color based on height (using face normal Y component for shading)
		var triColor color.RGBA
		if ti < len(chunk.FaceNormals) {
			n := chunk.FaceNormals[ti]
			// Directional lighting from upper-left
			light := Vec3{-0.5, 0.8, -0.3}.Normalize()
			intensity := n.Dot(light)
			if intensity < 0 {
				intensity = 0
			}
			intensity = 0.3 + 0.7*intensity // Ambient + diffuse

			// Base color varies with average height
			avgH := (chunk.Heights[t.A] + chunk.Heights[t.B] + chunk.Heights[t.C]) / 3.0
			baseR := 80.0 + avgH*2
			baseG := 120.0 + avgH*1.5
			baseB := 80.0 + avgH*0.5

			triColor = color.RGBA{
				R: uint8(clamp(baseR*intensity, 0, 255)),
				G: uint8(clamp(baseG*intensity, 0, 255)),
				B: uint8(clamp(baseB*intensity, 0, 255)),
				A: 255,
			}
		} else {
			triColor = color.RGBA{100, 140, 100, 255}
		}

		// Fill triangle
		drawFilledTriangle(screen, ax, ay, bx, by, cx, cy, triColor)

		// Draw edges
		// edgeColor := color.RGBA{60, 80, 60, 255}
		// vector.StrokeLine(screen, ax, ay, bx, by, 1, edgeColor, false)
		// vector.StrokeLine(screen, bx, by, cx, cy, 1, edgeColor, false)
		// vector.StrokeLine(screen, cx, cy, ax, ay, 1, edgeColor, false)
	}

	// // Draw sites as small dots
	// siteColor := color.RGBA{255, 200, 100, 255}
	// for _, idx := range chunk.CoreSiteIndices {
	// 	site := mesh.Sites[idx]
	// 	sx, sy := d.worldToScreen(site.Pos.X, site.Pos.Y)
	// 	vector.FillCircle(screen, sx, sy, 2, siteColor, false)
	// }
}

// drawChunkBoundary draws the chunk boundary.
func (d *DualsDemo) drawChunkBoundary(screen *ebiten.Image, chunk *TerrainChunk) {
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

// drawFilledTriangle fills a triangle using scanline approach.
func drawFilledTriangle(screen *ebiten.Image, x1, y1, x2, y2, x3, y3 float32, c color.RGBA) {
	// Use Ebiten's vector package to draw a filled triangle as a path
	var path vector.Path
	path.MoveTo(x1, y1)
	path.LineTo(x2, y2)
	path.LineTo(x3, y3)
	path.Close()

	vs, is := path.AppendVerticesAndIndicesForFilling(nil, nil)
	for i := range vs {
		vs[i].SrcX = 1
		vs[i].SrcY = 1
		vs[i].ColorR = float32(c.R) / 255.0
		vs[i].ColorG = float32(c.G) / 255.0
		vs[i].ColorB = float32(c.B) / 255.0
		vs[i].ColorA = float32(c.A) / 255.0
	}

	op := &ebiten.DrawTrianglesOptions{}
	op.FillRule = ebiten.FillRuleNonZero
	screen.DrawTriangles(vs, is, emptyImage(), op)
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
func (d *DualsDemo) drawRivers(screen *ebiten.Image, chunk *TerrainChunk) {
	if chunk.Hydro == nil {
		return
	}

	for _, segment := range chunk.Hydro.Rivers {
		if len(segment.Vertices) < 2 {
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
			sourceColor := color.RGBA{150, 200, 255, 200}
			vector.FillCircle(screen, sx, sy, 3, sourceColor, false)
		}
	}
}

// drawLakes renders lakes as filled water areas.
func (d *DualsDemo) drawLakes(screen *ebiten.Image, chunk *TerrainChunk) {
	if chunk.Hydro == nil || chunk.Mesh == nil {
		return
	}

	lakeColor := color.RGBA{40, 140, 160, 180}
	lakeEdgeColor := color.RGBA{60, 180, 200, 200}

	for _, lake := range chunk.Hydro.Lakes {
		if len(lake.SiteIndices) == 0 {
			continue
		}

		// Create a set of lake site indices for fast lookup
		lakeSites := make(map[int]bool)
		for _, idx := range lake.SiteIndices {
			lakeSites[idx] = true
		}

		// Find and fill triangles where all vertices are lake sites
		for _, tri := range chunk.Mesh.Tris {
			if lakeSites[tri.A] && lakeSites[tri.B] && lakeSites[tri.C] {
				a := chunk.Mesh.Sites[tri.A].Pos
				b := chunk.Mesh.Sites[tri.B].Pos
				c := chunk.Mesh.Sites[tri.C].Pos

				ax, ay := d.worldToScreen(a.X, a.Y)
				bx, by := d.worldToScreen(b.X, b.Y)
				cx, cy := d.worldToScreen(c.X, c.Y)

				drawFilledTriangle(screen, ax, ay, bx, by, cx, cy, lakeColor)
			}
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

// drawOcean renders ocean areas with submerged triangles and coastline.
func (d *DualsDemo) drawOcean(screen *ebiten.Image, chunk *TerrainChunk) {
	if chunk.Hydro == nil || chunk.Mesh == nil {
		return
	}

	ocean := chunk.Hydro.Ocean
	if !ocean.IsOcean {
		return
	}

	// Get the ocean region for sea level
	oceanRegion := d.Manager.HydroManager().GetOceanForChunk(chunk.Coord)
	seaLevel := float64(0)
	if oceanRegion != nil {
		seaLevel = oceanRegion.SeaLevel
	}

	// Draw submerged triangles
	for _, triID := range ocean.SubmergedTris {
		if triID < 0 || triID >= len(chunk.Mesh.Tris) {
			continue
		}

		tri := chunk.Mesh.Tris[triID]
		a := chunk.Mesh.Sites[tri.A].Pos
		b := chunk.Mesh.Sites[tri.B].Pos
		c := chunk.Mesh.Sites[tri.C].Pos

		ax, ay := d.worldToScreen(a.X, a.Y)
		bx, by := d.worldToScreen(b.X, b.Y)
		cx, cy := d.worldToScreen(c.X, c.Y)

		// Color based on depth (darker = deeper)
		avgHeight := (chunk.Heights[tri.A] + chunk.Heights[tri.B] + chunk.Heights[tri.C]) / 3.0
		depth := seaLevel - avgHeight
		depthFactor := clamp(depth/20.0, 0, 1) // Normalize depth to 0-1

		oceanColor := color.RGBA{
			R: uint8(30 - 20*depthFactor),
			G: uint8(80 - 30*depthFactor),
			B: uint8(140 + 40*depthFactor),
			A: uint8(180 + 50*depthFactor),
		}

		drawFilledTriangle(screen, ax, ay, bx, by, cx, cy, oceanColor)
	}

	// Draw coastline points
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
