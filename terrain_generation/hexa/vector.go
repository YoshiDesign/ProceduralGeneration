package hexa

func (v *Vec3) val() (x, y, z float32) {
	return v.X, v.Y, v.Z
}

func (v *Vec3) dot(u Vec3) float32 {

	return (v.X * u.X) + (v.Y * u.Y) + (v.Z * u.Z)

}
