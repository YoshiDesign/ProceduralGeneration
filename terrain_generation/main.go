package main

import (
	"log"
	"math/rand"
	"procedural_generation/terrain_generation/config"
	"procedural_generation/terrain_generation/hexa"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
)

type App struct {
	rng      *rand.Rand
	//ts 	 	 hexa.TileSet
	autoRun  bool
	accum    float64
	hexa	  hexa.HexaModule
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	//tiles := hexa.MakeChunkTiles()

	hexa_ := hexa.HexaModule{}
	hexa_.Init()

	return &App{
		//ts: 	  tiles,
		rng:      rng,
		autoRun:  false,
		hexa:	  hexa_,
	}

}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return config.GridW * config.CellSizePx, config.GridH * config.CellSizePx
}

func (a *App) Draw(screen *ebiten.Image) {

	a.hexa.Run(screen)

}


func (a *App) Update() error {

	return nil
}

func main() {
	ebiten.SetWindowTitle("Hexa")
	ebiten.SetWindowSize(config.GridW*config.CellSizePx, config.GridH*config.CellSizePx)

	app := NewApp()
	if err := ebiten.RunGame(app); err != nil {
		log.Fatal(err)
	}
}