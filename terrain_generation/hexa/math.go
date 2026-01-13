package hexa

func lerpFloat32(a, b, t float32) float32 {
	return a + (b-a)*t
}

func lerpInt(a, b, t int) int {
	return a + (b-a)*t
}

func absf(x float32) float32 {
	if x < 0 {
		return -x
	}
	return x
}

func (v *Vec3) val() (x, y, z float32) {
	return v.X, v.Y, v.Z
}

func (v *Vec3) dot(u Vec3) float32 {

	return (v.X * u.X) + (v.Y * u.Y) + (v.Z * u.Z)

}