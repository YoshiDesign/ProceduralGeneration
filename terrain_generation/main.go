package main

import (
	"image/color"
	"log"
	"math/rand"
	"procedural_generation/terrain_generation/hexa"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

var (
	GridW = 48
	GridH = 32
	CellSizePx = 32
)

var gridColor color.Color = color.RGBA{
	R: 255,
	G: 0,
	B: 0,
	A: 255,
}
var clr color.Color
var pixel *ebiten.Image

type App struct {
	rng      *rand.Rand
	//ts 	 	 hexa.TileSet
	autoRun  bool
	accum    float64
	hexGrid  []hexa.Hexagon
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	//tiles := hexa.MakeChunkTiles()

	grid := hexa.MakeHexField(GridW / 5, GridH / 5, float32(CellSizePx))

	return &App{
		//ts: 	  tiles,
		rng:      rng,
		autoRun:  false,
		hexGrid:  grid,
	}

}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return GridW * CellSizePx, GridH * CellSizePx
}

func (a *App) Draw(screen *ebiten.Image) {

	offsetX := float32(GridW*CellSizePx / 3)
	offsetY := float32(GridH*CellSizePx / 3)

	for n, cell := range a.hexGrid {
		corners := cell.CornerWorld // []Vec3, length 6

		if len(corners) < 2 {
			continue
		}

		for i := 0; i < len(corners); i++ {
			j := (i + 1) % len(corners)

			x1, y1 := worldToScreenXZ(corners[i]) 
			x2, y2 := worldToScreenXZ(corners[j])

			vector.StrokeLine(screen, 
				float32(x1) + (float32(CellSizePx)) + offsetX, 
				float32(y1) + (float32(CellSizePx)) + offsetY, 
				float32(x2) + (float32(CellSizePx)) + offsetX, 
				float32(y2) + (float32(CellSizePx)) + offsetY, 1.0, gridColor, false)


			// Draw interior lattices
		}

		for _, point := range cell.LocalLattice {

			world := hexa.Vec3 {
				X: cell.CenterWorld.X + point.Pos.X,
				Y: 0,
				Z: cell.CenterWorld.Z + point.Pos.Z,
			}

			sx, sy := worldToScreenXZ(world)

			var c color.RGBA
			if n % 2 == 1 {
				c = color.RGBA{0, 0, 255, 255}
			} else {
				c = color.RGBA{255, 255, 255, 255}
			}

			vector.FillCircle(screen, sx + float32(CellSizePx) + offsetX, sy + float32(CellSizePx) + offsetY, 1, c, false)

		}
			
	}

}

func worldToScreenXZ(v hexa.Vec3) (float32, float32) {
	// map (world X, world Z) -> (screen X, screen Y)
	return float32(v.X), float32(v.Z)
}

func (a *App) Update() error {

	return nil
}

func main() {
	ebiten.SetWindowTitle("Hexa")
	ebiten.SetWindowSize(GridW*CellSizePx, GridH*CellSizePx)

	app := NewApp()
	if err := ebiten.RunGame(app); err != nil {
		log.Fatal(err)
	}
}