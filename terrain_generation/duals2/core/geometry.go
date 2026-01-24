package core

import (
	"math"
	"sort"
)

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

// TriangleNormal computes the face normal for a triangle given heights at each site.
// The normal is computed using the cross product of two edges in 3D space.
// Returns (normal, ok) where ok=false if the triangle is degenerate.
//
// For terrain: X = east, Y = up (height), Z = north.
// CCW winding produces an upward-facing normal.
func (m *DelaunayMesh) TriangleNormal(triID int, heights []float64) (Vec3, bool) {
	if triID < 0 || triID >= len(m.Tris) {
		return Vec3{}, false
	}
	if len(heights) < len(m.Sites) {
		return Vec3{}, false
	}

	t := m.Tris[triID]

	// Build 3D positions: (Pos.X, height, Pos.Y) where Pos.Y is Z in world space
	a := Vec3{m.Sites[t.A].Pos.X, heights[t.A], m.Sites[t.A].Pos.Y}
	b := Vec3{m.Sites[t.B].Pos.X, heights[t.B], m.Sites[t.B].Pos.Y}
	c := Vec3{m.Sites[t.C].Pos.X, heights[t.C], m.Sites[t.C].Pos.Y}

	// Two edges from vertex A
	ab := b.Sub(a)
	ac := c.Sub(a)

	// Cross product gives the normal (CCW winding means this points "up")
	n := ab.Cross(ac)
	length := n.Len()
	if length < 1e-12 {
		return Vec3{0, 1, 0}, false // Degenerate, return default up
	}

	return n.Mul(1.0 / length), true
}

// TriangleNormalFromGradient computes the face normal from the height gradient.
// This is an alternative method using the pre-computed gradient (dhdx, dhdz).
//
// The gradient represents the slope in X and Z directions. The normal is:
//   n = normalize(-dhdx, 1, -dhdz)
//
// This method is efficient when you already have the gradient.
func TriangleNormalFromGradient(dhdx, dhdz float64) Vec3 {
	n := Vec3{-dhdx, 1.0, -dhdz}
	return n.Normalize()
}

// AllFaceNormals computes face normals for all triangles in the mesh.
// Returns a slice parallel to m.Tris containing the normalized face normal for each triangle.
func (m *DelaunayMesh) AllFaceNormals(heights []float64) []Vec3 {
	normals := make([]Vec3, len(m.Tris))
	for i := range m.Tris {
		n, ok := m.TriangleNormal(i, heights)
		if !ok {
			normals[i] = Vec3{0, 1, 0} // Default up for degenerate triangles
		} else {
			normals[i] = n
		}
	}
	return normals
}

// SlopeAngle returns the slope angle in radians for a triangle.
// 0 = flat, π/2 = vertical cliff.
func (m *DelaunayMesh) SlopeAngle(triID int, heights []float64) (float64, bool) {
	n, ok := m.TriangleNormal(triID, heights)
	if !ok {
		return 0, false
	}

	// Slope angle is the angle between the normal and the vertical (0, 1, 0)
	// cos(angle) = n · (0,1,0) = n.Y
	cosAngle := n.Y
	if cosAngle > 1.0 {
		cosAngle = 1.0
	}
	if cosAngle < -1.0 {
		cosAngle = -1.0
	}
	return math.Acos(cosAngle), true
}

// SlopePercent returns the slope as a percentage (rise/run * 100).
// A 45° slope returns 100%.
func (m *DelaunayMesh) SlopePercent(triID int, heights []float64) (float64, bool) {
	dhdx, dhdz, ok := m.TriangleGradient(triID, heights)
	if !ok {
		return 0, false
	}
	// Slope magnitude = sqrt(dhdx² + dhdz²)
	slope := math.Sqrt(dhdx*dhdx + dhdz*dhdz)
	return slope * 100.0, true
}

func (m *DelaunayMesh) findAnyOutgoing(site SiteIndex) int {
	for i := range m.HalfEdges {
		if m.HalfEdges[i].Origin == site {
			return i
		}
	}
	return -1
}

func angleSortAround(center Vec2, tris []int, verts []Vec2) {
	type item struct {
		tri int
		v   Vec2
		a   float64
	}
	items := make([]item, len(verts))
	for i := range verts {
		d := verts[i].Sub(center)
		items[i] = item{
			tri: tris[i],
			v:   verts[i],
			a:   math.Atan2(d.Y, d.X),
		}
	}
	sort.Slice(items, func(i, j int) bool { return items[i].a < items[j].a })
	for i := range items {
		tris[i] = items[i].tri
		verts[i] = items[i].v
	}
}

// VoronoiCellForSite walks the 1-ring around a site and returns
// the ordered circumcenters of adjacent triangles (the polygon vertices).
//
// If the site touches the boundary of the triangulation, there will be a missing twin,
// and the cell is “open” (infinite in true Voronoi). We return Closed=false.
func (m *DelaunayMesh) VoronoiCellForSite(site SiteIndex) VoronoiCell {
	start := -1
	if site >= 0 && int(site) < len(m.SiteEdge) {
		start = m.SiteEdge[site]
	}
	if start == -1 {
		return VoronoiCell{Site: site, Closed: false}
	}

	// We need a half-edge whose Origin == site
	// (SiteEdge should satisfy this, but be defensive)
	if m.HalfEdges[start].Origin != site {
		start = m.findAnyOutgoing(site)
		if start == -1 {
			return VoronoiCell{Site: site, Closed: false}
		}
	}

	// If site is on boundary, the ring walk using twin.next will stop.
	tris := make([]int, 0, 8)
	verts := make([]Vec2, 0, 8)

	visited := make(map[int]struct{}, 16)
	e := start
	closed := true

	for {
		if _, ok := visited[e]; ok {
			break
		}
		visited[e] = struct{}{}

		ti := m.HalfEdges[e].Tri
		tris = append(tris, ti)
		verts = append(verts, m.TriCircumcenter[ti])

		tw := m.HalfEdges[e].Twin
		if tw == -1 {
			closed = false
			break
		}
		// Move to next edge around the site: twin.next keeps Origin at same site.
		e = m.HalfEdges[tw].Next

		// Safety: ensure we’re still around the same site; if not, bail.
		if m.HalfEdges[e].Origin != site {
			closed = false
			break
		}
	}

	// Optional: enforce CCW order around site center by sorting by angle.
	// Ring-walk usually already provides correct order if triangle winding is consistent,
	// but sorting gives stability during early prototyping.
	angleSortAround(m.Sites[site].Pos, tris, verts)

	return VoronoiCell{
		Site:     site,
		Tris:     tris,
		Vertices: verts,
		Closed:   closed,
	}
}
