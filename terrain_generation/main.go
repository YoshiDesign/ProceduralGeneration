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
	GridW = 64
	GridH = 48
	CellSizePx = 10
)

var gridColor color.Color = color.RGBA{
	R: 255,
	G: 255,
	B: 255,
	A: 255,
}

type App struct {
	rng      *rand.Rand
	ts 	 	 hexa.TileSet
	autoRun  bool
	accum    float64
	hexGrid  []hexa.Hexagon
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	tiles := hexa.MakeChunkTiles()

	grid := hexa.MakeHexField(GridW, GridH, float32(CellSizePx))

	return &App{
		ts: 	  tiles,
		rng:      rng,
		autoRun:  false,
		hexGrid:  grid,
	}

}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return GridW * CellSizePx, GridH * CellSizePx
}

func (a *App) Draw(screen *ebiten.Image) {

	for _, cell := range a.hexGrid {
		corners := cell.CornerWorld // []Vec3, length 6

		if len(corners) < 2 {
			continue
		}

		for i := 0; i < len(corners); i++ {
			j := (i + 1) % len(corners)

			x1, y1 := worldToScreenXZ(corners[i]) 
			x2, y2 := worldToScreenXZ(corners[j])

			vector.StrokeLine(screen, x1, y1, x2, y2, 1.0, gridColor, false)
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