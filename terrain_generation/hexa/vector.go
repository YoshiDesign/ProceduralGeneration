package hexa

type Vec3 struct {
	X float32
	Y float32
	Z float32
}

func (v *Vec3) val() (x, y, z float32) {
	return v.X, v.Y, v.Z
}

func (v *Vec3) dot(u Vec3) float32 {

	return (v.X * u.X) + (v.Y * u.Y) + (v.Z * u.Z)

}

// Flat x,z plane
type Vec2 struct {
	x float32
	z float32
}