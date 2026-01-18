package main

import (
	"log"
	"math/rand"
	"procedural_generation/terrain_generation/config"
	"procedural_generation/terrain_generation/duals2"
	"procedural_generation/terrain_generation/hexa"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

// ViewMode determines which visualization is active.
type ViewMode int

const (
	ViewHexa ViewMode = iota
	ViewDuals
)

type App struct {
	rng      *rand.Rand
	autoRun  bool
	accum    float64
	hexa     hexa.HexaModule
	duals    *duals2.DualsDemo
	viewMode ViewMode
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))

	screenW := config.GridW * config.CellSizePx
	screenH := config.GridH * config.CellSizePx

	hexa_ := hexa.HexaModule{}
	hexa_.Init()

	dualsDemo := duals2.NewDualsDemo(screenW, screenH)

	return &App{
		rng:      rng,
		autoRun:  false,
		hexa:     hexa_,
		duals:    dualsDemo,
		viewMode: ViewDuals, // Start with Delaunay terrain view
	}
}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return config.GridW * config.CellSizePx, config.GridH * config.CellSizePx
}

func (a *App) Draw(screen *ebiten.Image) {
	switch a.viewMode {
	case ViewHexa:
		a.hexa.Run(screen)
	case ViewDuals:
		a.duals.Draw(screen)
	}
}

func (a *App) Update() error {
	// Toggle view mode with Tab key
	if inpututil.IsKeyJustPressed(ebiten.KeyTab) {
		if a.viewMode == ViewHexa {
			a.viewMode = ViewDuals
		} else {
			a.viewMode = ViewHexa
		}
	}

	// Update the active module
	if a.viewMode == ViewDuals {
		a.duals.Update()
	}

	return nil
}

func main() {
	ebiten.SetWindowTitle("Terrain Generation - Press Tab to toggle views")
	ebiten.SetWindowSize(config.GridW*config.CellSizePx, config.GridH*config.CellSizePx)

	app := NewApp()
	if err := ebiten.RunGame(app); err != nil {
		log.Fatal(err)
	}
}