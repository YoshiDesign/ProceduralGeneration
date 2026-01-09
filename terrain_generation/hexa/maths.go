package hexa

type Vec3 struct {
	x float32
	y float32
	z float32
}

func (v *Vec3) val() (x, y, z float32) {
	return v.x, v.y, v.z
}

func (v *Vec3) dot(u Vec3) float32 {

	return (v.x * u.x) + (v.y * u.y) + (v.z * u.z)

}