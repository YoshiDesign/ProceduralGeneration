package duals

import "math"

// Barycentric returns barycentric weights (wa, wb, wc) for point p in triangle triID.
// Triangle vertices are the *sites* (A,B,C). The weights satisfy:
//
//   p = wa*A + wb*B + wc*C
//   wa + wb + wc = 1
//
// ok=false if the triangle is degenerate (area ~ 0).
// If you want "inside triangle" test: wa>=0 && wb>=0 && wc>=0 (with epsilon).
func (m *DelaunayMesh) Barycentric(triID int, p Vec2) (wa, wb, wc float64, ok bool) {
	if triID < 0 || triID >= len(m.Tris) {
		return 0, 0, 0, false
	}

	t := m.Tris[triID]
	a := m.Sites[t.A].Pos
	b := m.Sites[t.B].Pos
	c := m.Sites[t.C].Pos

	// Using an efficient 2D barycentric formula.
	// denom = cross(b-a, c-a)
	v0 := b.Sub(a)
	v1 := c.Sub(a)
	v2 := p.Sub(a)

	denom := cross2(v0, v1) // 2*area signed
	const eps = 1e-12
	if math.Abs(denom) < eps {
		return 0, 0, 0, false
	}

	// wb = cross(v2, v1) / denom
	// wc = cross(v0, v2) / denom
	// wa = 1 - wb - wc
	wb = cross2(v2, v1) / denom
	wc = cross2(v0, v2) / denom
	wa = 1.0 - wb - wc
	return wa, wb, wc, true
}

// SampleScalar linearly interpolates a scalar field defined per site (vertex)
// across a triangle using barycentric weights.
//
// valuesAtSites must have length >= len(m.Sites).
// Returns (value, ok). ok=false if degenerate tri or invalid indexing.
func (m *DelaunayMesh) SampleScalar(triID int, p Vec2, valuesAtSites []float64) (float64, bool) {
	if len(valuesAtSites) < len(m.Sites) {
		return 0, false
	}
	if triID < 0 || triID >= len(m.Tris) {
		return 0, false
	}

	wa, wb, wc, ok := m.Barycentric(triID, p)
	if !ok {
		return 0, false
	}

	t := m.Tris[triID]
	va := valuesAtSites[t.A]
	vb := valuesAtSites[t.B]
	vc := valuesAtSites[t.C]

	return wa*va + wb*vb + wc*vc, true
}

// TriangleGradient returns the constant gradient of a linearly interpolated scalar field
// over the triangle (triID), assuming the scalar values are defined at the triangle's sites.
//
// For terrain height h(x,z), interpret:
//   dhdx = ∂h/∂x
//   dhdy = ∂h/∂y  (in your terrain: this is ∂h/∂z)
//
// This is extremely useful for slope computation because the gradient is constant per triangle.
//
// Returns (dhdx, dhdy, ok). ok=false if degenerate triangle or bad inputs.
func (m *DelaunayMesh) TriangleGradient(triID int, valuesAtSites []float64) (dhdx, dhdy float64, ok bool) {
	if len(valuesAtSites) < len(m.Sites) {
		return 0, 0, false
	}
	if triID < 0 || triID >= len(m.Tris) {
		return 0, 0, false
	}

	t := m.Tris[triID]
	a := m.Sites[t.A].Pos
	b := m.Sites[t.B].Pos
	c := m.Sites[t.C].Pos

	ha := valuesAtSites[t.A]
	hb := valuesAtSites[t.B]
	hc := valuesAtSites[t.C]

	// We solve for plane coefficients over 2D:
	// h(x,y) = px*x + py*y + k
	//
	// Using differences relative to a:
	// hb - ha = px*(bx-ax) + py*(by-ay)
	// hc - ha = px*(cx-ax) + py*(cy-ay)
	//
	// In matrix form:
	// [ (bx-ax) (by-ay) ] [px] = [hb-ha]
	// [ (cx-ax) (cy-ay) ] [py]   [hc-ha]
	//
	// Solve via inverse of 2x2. Determinant is cross(b-a, c-a).
	ab := b.Sub(a)
	ac := c.Sub(a)

	det := cross2(ab, ac)
	const eps = 1e-12
	if math.Abs(det) < eps {
		return 0, 0, false
	}

	rhs1 := hb - ha
	rhs2 := hc - ha

	// Inverse of [[ab.x ab.y],[ac.x ac.y]] is (1/det) * [[ ac.y, -ab.y],[-ac.x, ab.x]]
	dhdx = (rhs1*ac.Y - rhs2*ab.Y) / det
	dhdy = (-rhs1*ac.X + rhs2*ab.X) / det
	return dhdx, dhdy, true
}

// cross2 returns the 2D scalar cross product (a.x*b.y - a.y*b.x).
// This equals the signed area scale used by many barycentric/tri tests.
func cross2(a, b Vec2) float64 {
	return a.X*b.Y - a.Y*b.X
}
