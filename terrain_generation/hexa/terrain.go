package hexa

import "math"

const N = 8 // Subtile lattice area

func MakeHexField(width, height int, size float32) []Hexagon {

	var field []Hexagon

	for i := range height {
		for j := range width {

			hex := Hexagon{
				X:            j,
				Y:            0,
				Z: 			  i,
				Size_radius:  size,
				inner_radius: size * float32(math.Sqrt(3)) / 2,
				layout:       PointyTop,
				CenterWorld:  Vec3{0,0,0},
				CornerWorld:  	make([]Vec3, CornerCount), // All 6 corners in world space. TODO - are they winding cw or ccw? I think ccw
				LocalLattice:  	make([]FineLatticePoint, 0, N*N),
				LatticeToPoint: make([]int, N*N),
				PointToLattice: make([]FineGridIJ, 0, N*N),
			}

			hex.CenterWorld = hex.center()

			for c := CornerI(0); c < CornerCount; c++ {
				hex.CornerWorld[c] = hex.corner(c)
			}

			// Populate fine lattice representation with invalid indices
			for k := range hex.LatticeToPoint {
				hex.LatticeToPoint[k] = -1
			}

			// Create the hexagon's tile-local lattice
			hex.makeLattice8x8()

			field = append(field, hex)
		}
	}

	return field

}

func MakeWorld(cells []Hexagon) {

}

// func MakeTerrainHeightMap(width, height int) []Vec3 {

// }