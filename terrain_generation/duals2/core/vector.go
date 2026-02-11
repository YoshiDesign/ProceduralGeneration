package core

import "math"

type Vec2 struct{ X, Y float64 }

func (a Vec2) Add(b Vec2) Vec2    { return Vec2{a.X + b.X, a.Y + b.Y} }
func (a Vec2) Sub(b Vec2) Vec2    { return Vec2{a.X - b.X, a.Y - b.Y} }
func (a Vec2) Mul(s float64) Vec2 { return Vec2{a.X * s, a.Y * s} }
func (a Vec2) Dot(b Vec2) float64 { return a.X*b.X + a.Y*b.Y }
func (a Vec2) Len2() float64      { return a.Dot(a) }
func (a Vec2) Len() float64       { return math.Sqrt(a.Len2()) }
func (a Vec2) Eq(b Vec2) bool {
	return a.X == b.X && a.Y == b.Y
}

// Normalize returns a unit vector (or zero vector if length is zero).
func (v Vec2) Normalize() Vec2 {
	l := v.Len()
	if l < 1e-12 {
		return Vec2{}
	}
	return v.Mul(1.0 / l)
}

// Vec3 represents a 3D point/vector (X = east, Y = up, Z = north).
type Vec3 struct{ X, Y, Z float64 }

func (a Vec3) Add(b Vec3) Vec3    { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }
func (a Vec3) Sub(b Vec3) Vec3    { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }
func (a Vec3) Mul(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }
func (a Vec3) Dot(b Vec3) float64 { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }
func (a Vec3) Len2() float64      { return a.Dot(a) }
func (a Vec3) Len() float64       { return math.Sqrt(a.Len2()) }
func (a Vec3) Eq(b Vec3) bool {
	return a.X == b.X && a.Y == b.Y && a.Z == b.Z
}

func (a Vec3) Cross(b Vec3) Vec3 {
	return Vec3{
		a.Y*b.Z - a.Z*b.Y,
		a.Z*b.X - a.X*b.Z,
		a.X*b.Y - a.Y*b.X,
	}
}

/* TODO: Midnight's default up is -Y */
func (a Vec3) Normalize() Vec3 {
	l := a.Len()
	if l < 1e-12 {
		return Vec3{0, 1, 0} // default up
	}
	return a.Mul(1.0 / l)
}

// XZ returns the horizontal (X, Z) components as a Vec2 (for 2D operations).
func (a Vec3) XZ() Vec2 { return Vec2{a.X, a.Z} }
