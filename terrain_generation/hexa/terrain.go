package hexa

import "math"

func MakeHexField(width, height int, size float32) []Hexagon {

	var field []Hexagon

	for i := range height {
		for j := range width {

			hex := Hexagon{
				x:            j,
				y:            0,
				z: 			  i,
				size_radius:  size,
				inner_radius: size * float32(math.Sqrt(3)) / 2,
				layout:       FlatTop,
				CenterWorld:  Vec3{0,0,0},
				CornerWorld:  make([]Vec3, CornerCount),
			}

			hex.CenterWorld = hex.center()

			for c := CornerI(0); c < CornerCount; c++ {
				hex.CornerWorld[c] = hex.corner(c)
			}

			field = append(field, hex)
		}
	}

	return field

}

func MakeWorld(cells []Hexagon) {

}

// func MakeTerrainHeightMap(width, height int) []Vec3 {

// }