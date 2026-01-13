package hexa

import (
	"image/color"
	"procedural_generation/terrain_generation/config"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

var	offsetX = float32(config.GridW*config.CellSizePx / 3)
var	offsetY = float32(config.GridH*config.CellSizePx / 3)

type HexaModule struct {
	Grid []Hexagon
}

func worldToScreenXZ(v Vec3) (float32, float32) {
	// map (world X, world Z) -> (screen X, screen Y)
	return float32(v.X), float32(v.Z)
}

func (mod *HexaModule) Init() {
	mod.Grid = MakeHexField(
		config.GridW/5, 
		config.GridH/5, 
		float32(config.CellSizePx))
}

func (mod *HexaModule) Run(screen *ebiten.Image) {

for n, cell := range mod.Grid {
		corners := cell.CornerWorld // []Vec3, length 6

		if len(corners) < 2 {
			continue
		}

		for i := 0; i < len(corners); i++ {
			j := (i + 1) % len(corners)

			x1, y1 := worldToScreenXZ(corners[i]) 
			x2, y2 := worldToScreenXZ(corners[j])

			vector.StrokeLine(screen, 
				float32(x1) + (float32(config.CellSizePx)) + offsetX, 
				float32(y1) + (float32(config.CellSizePx)) + offsetY, 
				float32(x2) + (float32(config.CellSizePx)) + offsetX, 
				float32(y2) + (float32(config.CellSizePx)) + offsetY, 1.0, config.GridColor, false)

			// Draw interior lattices
		}

		for _, point := range cell.LocalLattice {

			world := Vec3 {
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

			vector.FillCircle(screen, 
				sx + float32(config.CellSizePx) + offsetX, 
				sy + float32(config.CellSizePx) + offsetY, 
				1, 
				c, 
				false)

		}
			
	}

}