package hexa

import (
	"log"
	"math"
)

/*
	OPTIMIZATION:

	1. The entire Fine Lattice topology data (PointToLattice, LatticeToPoint and LocalLattice) can be made global
	   Right now each hexagon carries their own data.

	   For LatticeToPoint, that's simply a mask. Computing it for every hexagon is a very redundant thing to do.
	   Do it once, make it read only; prosper.

*/

/*
Mental Model:
	World
	└── Hex tiles (axial coords)
		└── Fine lattice (gridToPoint, points)
			└── Hydrology / slope / erosion
*/

/*
	In the flat-top orientation, the horizontal distance between adjacent hexagons centers is:
 	horiz = 3/4 * width = 3/2 * Size_radius.
	The vertical distance is vert = height = sqrt(3) * Size_radius = 2 * inradius.

	In the pointy-top orientation, the horizontal distance between adjacent hexagon centers is:
 	horiz = width == sqrt(3) * Size_radius == 2 * inradius.
	The vertical distance is vert = 3/4 * height == 3/2 * Size_radius.
*/




func makeHexagon (
	x_, y_, z_ int, 
	circumRad float32, 
	layout_ HexType) Hexagon {

	hex := Hexagon{
		X:            x_,
		Y:            y_,
		Z: 			  z_,
		Size_radius:  circumRad,
		inner_radius: circumRad * float32(math.Sqrt(3)) / 2,
		layout:       FlatTop,
		CenterWorld:  Vec3{0,0,0},
		CornerWorld:  make([]Vec3, CornerCount),
	}

	hex.CenterWorld = hex.center()

	for c := CornerI(0); c < CornerCount; c++ {
		hex.CornerWorld[c] = hex.corner(c)
	}

	return hex
}

func (t HexType) String() string {

	switch t {
	case PointyTop: 
		return "pointy"
	case FlatTop  : 
		return "flat"
	}

	return "unknown"
}

func (t HexType) IsValid() bool {

	switch t {
	case PointyTop, FlatTop:
		return true
	default:
		return false
	}

}

func (hex *Hexagon) outerRadius() float32 {
	return hex.Size_radius
}

func (hex *Hexagon) innerRadius() float32 {
	return hex.inner_radius
}

/* Using the cartesian coordinates, define the center point in euclidean space */
func (h Hexagon) center() Vec3 {
	R := h.Size_radius
	q := float32(h.X) // Initialized from Axial space
	r := float32(h.Z) // Initialized from Axial space

	sqrt3 := float32(math.Sqrt(3))

	switch h.layout {
	case PointyTop:
		// Transcribe to euclidean space
		cx := R * sqrt3 * (q + 0.5*float32(h.Z&1))
		//cx := R * sqrt3 * (q + r /2)
		cz := R * 1.5 * r
		return Vec3{cx, float32(h.Y), cz}
	case FlatTop:
		cx := R * 1.5 * q
		//cz := R * sqrt3 * (q + r /2)
		cz := R * sqrt3 * (r + 0.5*float32(h.X&1))
		return Vec3{cx, float32(h.Y), cz}
	default:
		return Vec3{0, 0, 0}
	}
}

/* Corner coordinates are in euclidean space simply because h.center() derives the center coordinate in euclidean space */
func (h Hexagon) corner(i CornerI) Vec3 {

	if i >= CornerCount {
		log.Fatal("too many angle indices")
	}

	center := h.center()
	var angle_deg float32

	switch h.layout {
	case PointyTop:
		angle_deg = 60.0 * float32(i) - 30.0
	case FlatTop:
		angle_deg = 60.0 * float32(i) 
	default:
		log.Fatal("invalid layout")
	}

	var angle_rad = math.Pi / 180.0 * angle_deg

	return Vec3{center.X + h.Size_radius * float32(math.Cos(float64(angle_rad))),
			float32(h.Y),
			center.Z + h.Size_radius * float32(math.Sin(float64(angle_rad)))}

}

func (h Hexagon) slopeTo(h2 Hexagon) float32 {
	
	if !h.congruent(h2) {
		log.Fatal("incongruent hexagons")
	}

	hpos := h.center()
	h2pos := h2.center()

	dx := h2pos.X - hpos.X
	dz := h2pos.Z - hpos.Z

	if dx == 0 {
		return float32(math.Inf(1))
	}
	return dz / dx
}

func (hex *Hexagon) congruent(hex2 Hexagon) bool {
	return hex.Size_radius == hex2.Size_radius && hex.inner_radius == hex2.inner_radius
}

/* `pos` are the axial coordinates */
func ToWorldSpace(pos Vec3, cRad float32, t HexType) Vec3 {
	R := cRad
	q := float32(pos.X)	
	r := float32(pos.Z)

	sqrt3 := float32(math.Sqrt(3))

	switch t {
	case PointyTop:
		// Transcribe to euclidean space
		cx := R * sqrt3 * (q + r/2.0)
		cz := R * 1.5 * r
		return Vec3{cx, pos.Y, cz}
	case FlatTop:
		cx := R * 1.5 * q
		cz := R * sqrt3 * (r + q/2.0)
		return Vec3{cx, pos.Y, cz}
	default:
		return Vec3{0,0,0}
	}

}

/* 
	We could be pre-computing certain weights during this method based on external parameters (tile's world position)
	That acts to make this less generic, more intent specific, but saves us on compute
*/
func (h *Hexagon) makeLattice8x8() {

	lattice := make([]FineLatticePoint, 0, N*N)
	const N = 8

	switch h.layout {

	case PointyTop:
		for j := range N {
			for i := range N {
				// Points in continuous tile-local space. Row major order				
				u := lerpFloat32(-h.Size_radius,  h.Size_radius,  (float32(i)+0.5)/N ) // cols
				v := lerpFloat32(-h.inner_radius, h.inner_radius, (float32(j)+0.5)/N ) // rows
				y := float32(0.0)

				if pointInsideHexPointy(u, v, h.Size_radius) {

					idx := len(lattice)

					// Note: u = x, v = z to match our coordinate system (right handed, +z forward)
					lattice = append(lattice, FineLatticePoint{
						Pos: Vec3{u, y, v},
						weights: WeightField{},
					})

					gridIdx := j*N + i
					h.LatticeToPoint[gridIdx] = idx
					h.PointToLattice = append(h.PointToLattice, FineGridIJ{I: i, J: j}) 
				}
			}
		}

	case FlatTop:
		for j := range N {
			for i := range N {
				// Points in continuous tile-local space. Row major order					
				u := lerpFloat32(-h.inner_radius, h.inner_radius, (float32(j)+0.5)/N ) // cols
				v := lerpFloat32(-h.Size_radius,  h.Size_radius,  (float32(i)+0.5)/N ) // rows
				y := float32(0.0)

				if pointInsideHexFlat(u, v, h.Size_radius) {
					idx := len(lattice)

					// Note: u = x, v = z to match our coordinate system (right handed, +z forward)
					lattice = append(lattice, FineLatticePoint{
						Pos: Vec3{u, y, v},
						weights: WeightField{},
					})

					gridIdx := j*N + i
					h.LatticeToPoint[gridIdx] = idx
					h.PointToLattice = append(h.PointToLattice, FineGridIJ{I: i, J: j}) 
				}

			}
		}

	default:
		log.Fatal("invalid hexagon layout")
	}

	h.LocalLattice = lattice

}

func pointInsideHexPointy(x, z, R float32) bool {
	const sqrt3 = float32(1.7320508075688772)

	ax := absf(x)
	az := absf(z)

	if ax > R {
		return false
	}
	if az > (sqrt3/2)*R {
		return false
	}
	if ax/2 + az/sqrt3 > R {
		return false
	}
	return true
}


func pointInsideHexFlat(x, z, R float32) bool {
	const sqrt3 = float32(1.7320508075688772)

	ax := absf(x)
	az := absf(z)

	if az > R {
		return false
	}
	if ax > (sqrt3/2)*R {
		return false
	}
	if az/2 + ax/sqrt3 > R {
		return false
	}
	return true
}
