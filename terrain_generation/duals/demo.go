package duals

import (
	"image/color"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

// DualsDemo is a visualization module for the Delaunay terrain system.
type DualsDemo struct {
	Manager     *ChunkManager   // Chunk manager with caching
	Chunks      []*TerrainChunk // Active chunks for rendering
	ScreenW     int
	ScreenH     int
	Scale       float64 // pixels per world unit
	OffsetX     float64
	OffsetZ     float64
	ShowVoronoi bool
	ShowNormals bool
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
        total += exampleHeightFunc(x*frequency, z*frequency) * amplitude
        maxValue += amplitude
        amplitude *= persistence   // Each octave is quieter
        frequency *= lacunarity    // Each octave is higher frequency
    }
    
    return total / maxValue  // Normalize to [-1, 1]
}

// NewDualsDemo creates a demo with a 2x2 grid of chunks using the ChunkManager.
func NewDualsDemo(screenW, screenH int) *DualsDemo {
	cfg := ChunkConfig{
		ChunkSize:    128.0,
		MinPointDist: 12.0,
		HaloWidth:    12.0,
		WorldSeed:    12345,
	}

	// Height function: simple noise-like terrain
	heightFunc := func(x, z float64) float64 {
		// Multi-octave sine waves for terrain-like appearance
		h := 0.0
		h += 15.0 * math.Sin(x*0.02) * math.Cos(z*0.02)
		h += 8.0 * math.Sin(x*0.05+1.0) * math.Cos(z*0.04+0.5)
		h += 4.0 * math.Sin(x*0.1+2.0) * math.Sin(z*0.08+1.0)
		return h
	}

	// Create chunk manager with caching
	manager := NewChunkManager(cfg, heightFunc)

	// Generate a 3x2 grid of chunks using the manager (benefits from caching)
	chunks := make([]*TerrainChunk, 0, 6)
	for cz := 0; cz < 2; cz++ {
		for cx := 0; cx < 3; cx++ {
			coord := ChunkCoord{X: cx, Z: cz}
			chunk, err := manager.GetOrGenerate(coord)
			if err != nil {
				continue
			}
			chunks = append(chunks, chunk)
		}
	}

	// Scale to fit screen
	totalWorldSize := cfg.ChunkSize * 2
	scale := float64(min(screenW, screenH)) / totalWorldSize * 0.9

	return &DualsDemo{
		Manager:     manager,
		Chunks:      chunks,
		ScreenW:     screenW,
		ScreenH:     screenH,
		Scale:       scale,
		OffsetX:     float64(screenW) * 0.05,
		OffsetZ:     float64(screenH) * 0.05,
		ShowVoronoi: false,
		ShowNormals: false,
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

	// Draw chunk boundaries
	for _, chunk := range d.Chunks {
		d.drawChunkBoundary(screen, chunk)
	}
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
		edgeColor := color.RGBA{60, 80, 60, 255}
		vector.StrokeLine(screen, ax, ay, bx, by, 1, edgeColor, false)
		vector.StrokeLine(screen, bx, by, cx, cy, 1, edgeColor, false)
		vector.StrokeLine(screen, cx, cy, ax, ay, 1, edgeColor, false)
	}

	// Draw sites as small dots
	siteColor := color.RGBA{255, 200, 100, 255}
	for _, idx := range chunk.CoreSiteIndices {
		site := mesh.Sites[idx]
		sx, sy := d.worldToScreen(site.Pos.X, site.Pos.Y)
		vector.FillCircle(screen, sx, sy, 2, siteColor, false)
	}
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
