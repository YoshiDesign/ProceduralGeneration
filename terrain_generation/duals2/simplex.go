package duals2

import "math"

// 2D Simplex Noise implementation based on Ken Perlin's algorithm.
// Returns values in the range [-1, 1].

// Skewing factors for 2D simplex grid transformation
var (
	// F2 skews from square to simplex grid: (sqrt(3) - 1) / 2
	f2 = 0.5 * (math.Sqrt(3.0) - 1.0)
	// G2 unskews from simplex back to square: (3 - sqrt(3)) / 6
	g2 = (3.0 - math.Sqrt(3.0)) / 6.0
)

// Gradient vectors for 2D noise (12 directions)
var grad2 = [12][2]float64{
	{1, 1}, {-1, 1}, {1, -1}, {-1, -1},
	{1, 0}, {-1, 0}, {0, 1}, {0, -1},
	{1, 1}, {-1, 1}, {1, -1}, {-1, -1},
}

// Permutation table (doubled to avoid overflow handling)
var perm = [512]int{
	151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
	140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148,
	247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
	57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
	74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122,
	60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
	65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169,
	200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64,
	52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
	207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213,
	119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104,
	218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,
	81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157,
	184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
	222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
	// Duplicate for overflow wrapping
	151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
	140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148,
	247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
	57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175,
	74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122,
	60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
	65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169,
	200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64,
	52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
	207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213,
	119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104,
	218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,
	81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157,
	184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
	222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
}

// dot2 computes the dot product of a gradient vector and a distance vector
func dot2(g [2]float64, x, y float64) float64 {
	return g[0]*x + g[1]*y
}

// fastFloor returns the largest integer less than or equal to x
func fastFloor(x float64) int {
	if x > 0 {
		return int(x)
	}
	return int(x) - 1
}

// Simplex2D computes 2D simplex noise at coordinates (x, z).
// Returns a value in the range [-1, 1].
func Simplex2D(x, z float64) float64 {
	var n0, n1, n2 float64 // Noise contributions from three corners

	// Skew the input space to determine which simplex cell we're in
	s := (x + z) * f2
	i := fastFloor(x + s)
	j := fastFloor(z + s)

	// Unskew the cell origin back to (x, z) space
	t := float64(i+j) * g2
	X0 := float64(i) - t
	Y0 := float64(j) - t
	x0 := x - X0 // Distance from cell origin in (x, z) space
	y0 := z - Y0

	// Determine which simplex we're in (upper or lower triangle)
	var i1, j1 int
	if x0 > y0 {
		i1, j1 = 1, 0 // Lower triangle, XY order: (0,0)->(1,0)->(1,1)
	} else {
		i1, j1 = 0, 1 // Upper triangle, YX order: (0,0)->(0,1)->(1,1)
	}

	// Offsets for middle corner in (x, z) unskewed coords
	x1 := x0 - float64(i1) + g2
	y1 := y0 - float64(j1) + g2

	// Offsets for last corner in (x, z) unskewed coords
	x2 := x0 - 1.0 + 2.0*g2
	y2 := y0 - 1.0 + 2.0*g2

	// Work out the hashed gradient indices of the three simplex corners
	ii := i & 255
	jj := j & 255
	gi0 := perm[ii+perm[jj]] % 12
	gi1 := perm[ii+i1+perm[jj+j1]] % 12
	gi2 := perm[ii+1+perm[jj+1]] % 12

	// Calculate the contribution from the three corners
	t0 := 0.5 - x0*x0 - y0*y0
	if t0 < 0 {
		n0 = 0.0
	} else {
		t0 *= t0
		n0 = t0 * t0 * dot2(grad2[gi0], x0, y0)
	}

	t1 := 0.5 - x1*x1 - y1*y1
	if t1 < 0 {
		n1 = 0.0
	} else {
		t1 *= t1
		n1 = t1 * t1 * dot2(grad2[gi1], x1, y1)
	}

	t2 := 0.5 - x2*x2 - y2*y2
	if t2 < 0 {
		n2 = 0.0
	} else {
		t2 *= t2
		n2 = t2 * t2 * dot2(grad2[gi2], x2, y2)
	}

	// Scale to return values in the range [-1, 1]
	return 70.0 * (n0 + n1 + n2)
}
