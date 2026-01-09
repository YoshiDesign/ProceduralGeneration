package main

import (
	"log"
	"math/rand"
	"time"

	"procedural_generation/terrain_generation/hexa"

	"github.com/hajimehoshi/ebiten/v2"
)

var (
	GridW = 10
	GridH = 10
	CellSizePx = 10
)

type App struct {
	rng      *rand.Rand
	ts 	 	 hexa.TileSet
	autoRun  bool
	accum    float64
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	tiles := hexa.MakeChunkTiles() // change this to your own tile set later

	return &App{
		ts: 	  tiles,
		rng:      rng,
		autoRun:  false,
	}

}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return GridW * CellSizePx, GridH * CellSizePx
}

func (a *App) Draw(screen *ebiten.Image) {

}

func (a *App) Update() error {

	return nil
}

func main() {
	ebiten.SetWindowTitle("WFC Prototype (Ebitengine)")
	ebiten.SetWindowSize(GridW*CellSizePx, GridH*CellSizePx)

	app := NewApp()
	if err := ebiten.RunGame(app); err != nil {
		log.Fatal(err)
	}
}