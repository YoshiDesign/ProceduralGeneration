package hexa

import (
	"log"
	"math"
)

/*
	In the flat-top orientation, the horizontal distance between adjacent hexagons centers is:
 	horiz = 3/4 * width = 3/2 * size_radius.
	The vertical distance is vert = height = sqrt(3) * size_radius = 2 * inradius.

	In the pointy-top orientation, the horizontal distance between adjacent hexagon centers is:
 	horiz = width == sqrt(3) * size_radius == 2 * inradius.
	The vertical distance is vert = 3/4 * height == 3/2 * size_radius.
*/

type HexType int
type CornerI int

const (
	PointyTop HexType = iota
	FlatTop
)

const (
	First CornerI = iota 
	Second
	Third
	Fourth
	Fifth
	Sixth
	CornerCount // Sentinel
)


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

type Hexagon struct {
	inner_radius float32
	size_radius float32
	x, y, z int // basic cartesian positions
	layout HexType
	CenterWorld Vec3
	CornerWorld []Vec3
}

func (hex *Hexagon) outerRadius() float32 {
	return hex.size_radius
}

func (hex *Hexagon) innerRadius() float32 {
	return hex.inner_radius
}

func (h Hexagon) center() Vec3 {
	R := h.size_radius
	q := float32(h.x)	
	r := float32(h.z)

	sqrt3 := float32(math.Sqrt(3))

	switch h.layout {
	case PointyTop:
		// Transcribe to euclidean space
		cx := R * sqrt3 * (q + 0.5*float32(h.z&1))
		cz := R * 1.5 * r
		return Vec3{cx, float32(h.y), cz}
	case FlatTop:
		cx := R * 1.5 * q
		cz := R * sqrt3 * (r + 0.5*float32(h.x&1))
		return Vec3{cx, float32(h.y), cz}
	default:
		return Vec3{0, 0, 0}
	}
}

/* Corner coordinates are in world space simply because h.center() derives the center coordinate in world space */
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

	return Vec3{center.X + h.size_radius * float32(math.Cos(float64(angle_rad))),
			float32(h.y),
			center.Z + h.size_radius * float32(math.Sin(float64(angle_rad)))}

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
	return hex.size_radius == hex2.size_radius && hex.inner_radius == hex2.inner_radius
}

/* `pos` are the axial coordinates */
func (h Hexagon) toWorldSpace(pos Vec3) Vec3 {
	R := h.size_radius
	q := float32(pos.X)	
	r := float32(pos.Z)

	sqrt3 := float32(math.Sqrt(3))

	switch h.layout {
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