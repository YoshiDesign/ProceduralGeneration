package main

import (
	"log"
	"math/rand"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/vector"

	"wavegen/wfc"
)

type App struct {
	sim      *wfc.Solver
	rng      *rand.Rand
	autoRun  bool
	stepRate int // steps per second when autoRun
	accum    float64
}

func NewApp() *App {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	tiles := wfc.MakeDemoTiles() // change this to your own tile set later

	s := wfc.NewSolver(wfc.GridW, wfc.GridH, tiles, rng)
	s.Reset()

	return &App{
		sim:      s,
		rng:      rng,
		autoRun:  false,
		stepRate: 30,
	}
}

func (a *App) Update() error {
	// Controls
	if ebiten.IsKeyPressed(ebiten.KeyR) {
		a.sim.Reset()
	}
	if ebiten.IsKeyPressed(ebiten.KeySpace) {
		// single step (edge-trigger would be nicer; fine for skeleton)
		a.sim.Step()
	}
	if ebiten.IsKeyPressed(ebiten.KeyEnter) {
		a.autoRun = true
	}
	if ebiten.IsKeyPressed(ebiten.KeyEscape) {
		a.autoRun = false
	}

	// Optional: left click forces a random tile into a cell (useful for testing propagation)
	if ebiten.IsMouseButtonPressed(ebiten.MouseButtonLeft) {
		x, y := ebiten.CursorPosition()
		cx := x / wfc.CellSizePx
		cy := y / wfc.CellSizePx
		if cx >= 0 && cx < wfc.GridW && cy >= 0 && cy < wfc.GridH {
			a.sim.ForceRandomTile(cx, cy)
		}
	}

	// Auto-run stepping
	if a.autoRun && !a.sim.Done() && !a.sim.Contradiction() {
		// Ebiten update is ~60hz by default
		stepsPerFrame := a.stepRate / 60
		if stepsPerFrame < 1 {
			stepsPerFrame = 1
		}
		for i := 0; i < stepsPerFrame; i++ {
			a.sim.Step()
			if a.sim.Done() || a.sim.Contradiction() {
				break
			}
		}
	}

	return nil
}

func (a *App) Draw(screen *ebiten.Image) {
	// Draw grid
	for y := 0; y < a.sim.H; y++ {
		for x := 0; x < a.sim.W; x++ {
			d := a.sim.DomainAt(x, y)
			c := a.sim.CellColor(d)

			vector.FillRect(screen,
				float32(x*wfc.CellSizePx),
				float32(y*wfc.CellSizePx),
				float32(wfc.CellSizePx-1),
				float32(wfc.CellSizePx-1),
				c,
				false,
			)
		}
	}

	// Draw Roads
	a.sim.RoadLines(screen)

	// HUD
	status := "SPACE=step  ENTER=run  ESC=stop  R=reset  LMB=force\n"
	status += "steps=" + itoa(a.sim.Steps) +
		" collapsed=" + itoa(a.sim.CollapsedCount()) + "/" + itoa(a.sim.W*a.sim.H)

	if a.sim.Contradiction() {
		status += "  [CONTRADICTION]"
	} else if a.sim.Done() {
		status += "  [DONE]"
	} else if a.autoRun {
		status += "  [RUNNING]"
	}

	ebitenutil.DebugPrint(screen, status)
}

func (a *App) Layout(outsideW, outsideH int) (int, int) {
	return wfc.GridW * wfc.CellSizePx, wfc.GridH * wfc.CellSizePx
}

func main() {
	ebiten.SetWindowTitle("WFC Prototype (Ebitengine)")
	ebiten.SetWindowSize(wfc.GridW*wfc.CellSizePx, wfc.GridH*wfc.CellSizePx)

	app := NewApp()
	if err := ebiten.RunGame(app); err != nil {
		log.Fatal(err)
	}
}

func itoa(v int) string {
	// tiny helper to avoid fmt in hot path; fine for demo
	if v == 0 {
		return "0"
	}
	neg := false
	if v < 0 {
		neg = true
		v = -v
	}
	var buf [32]byte
	i := len(buf)
	for v > 0 {
		i--
		buf[i] = byte('0' + (v % 10))
		v /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}
