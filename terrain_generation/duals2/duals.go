package duals

import (
	"fmt"
	"math"
	"sort"
)

type Vec2 struct{ X, Y float64 }

func (a Vec2) Add(b Vec2) Vec2   { return Vec2{a.X + b.X, a.Y + b.Y} }
func (a Vec2) Sub(b Vec2) Vec2   { return Vec2{a.X - b.X, a.Y - b.Y} }
func (a Vec2) Mul(s float64) Vec2 { return Vec2{a.X * s, a.Y * s} }
func (a Vec2) Dot(b Vec2) float64 { return a.X*b.X + a.Y*b.Y }
func (a Vec2) Len2() float64      { return a.Dot(a) }
func (a Vec2) Len() float64       { return math.Sqrt(a.Len2()) }

// Vec3 represents a 3D point/vector (X = east, Y = up, Z = north).
type Vec3 struct{ X, Y, Z float64 }

func (a Vec3) Add(b Vec3) Vec3    { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }
func (a Vec3) Sub(b Vec3) Vec3    { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }
func (a Vec3) Mul(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }
func (a Vec3) Dot(b Vec3) float64 { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }
func (a Vec3) Len2() float64      { return a.Dot(a) }
func (a Vec3) Len() float64       { return math.Sqrt(a.Len2()) }

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

// Site represents a terrain vertex with position and height.
type Site struct {
	Pos    Vec2    // 2D position (X, Z plane)
	Height float64 // elevation (Y axis)
}

// Pos3D returns the site position as a 3D point with height as Y.
func (s Site) Pos3D() Vec3 {
	return Vec3{s.Pos.X, s.Height, s.Pos.Y}
}

// Triangle references 3 site indices in CCW order (recommended).
type Triangle struct{ A, B, C int }

// HalfEdge is a directed edge: Origin -> Dest,
// where Dest is the Origin of Next.
type HalfEdge struct {
	Origin   int // site index
	Tri      int // triangle index
	Next     int // half-edge index (within triangle cycle)
	Twin     int // opposite direction half-edge index, or -1 if boundary
	Prev     int // (optional) convenient; set during build
	EdgeDest int // cached dest site (optional, set during build)
}

// DelaunayMesh holds topology + geometry.
// It’s “Delaunay” if the triangles were produced by a Delaunay triangulation.
type DelaunayMesh struct {
	Sites     []Site
	Tris      []Triangle
	HalfEdges []HalfEdge

	// For each triangle, index of its first half-edge (3 consecutive edges)
	TriEdge0 []int

	// For each site, one outgoing half-edge from that site (or -1)
	SiteEdge []int

	// Voronoi vertices: one per triangle
	TriCircumcenter []Vec2
	TriIsDegenerate []bool
}

// BuildHalfEdgeMesh constructs connectivity from Sites + Tris.
// Triangles must reference valid site indices.
// If triangles aren’t consistently wound, Voronoi polygons can come out flipped/odd.
func BuildHalfEdgeMesh(sites []Site, tris []Triangle) (*DelaunayMesh, error) {
	m := &DelaunayMesh{
		Sites:            sites,
		Tris:             tris,
		TriEdge0:         make([]int, len(tris)),
		SiteEdge:         make([]int, len(sites)),
		TriCircumcenter:  make([]Vec2, len(tris)),
		TriIsDegenerate:  make([]bool, len(tris)),
	}
	for i := range m.SiteEdge {
		m.SiteEdge[i] = -1
	}

	// Allocate 3 half-edges per triangle.
	m.HalfEdges = make([]HalfEdge, 0, len(tris)*3)

	// Edge map for twin linking: key = (min,max,dir?) but we need direction.
	// We'll store directed edge (u->v) and look up opposite (v->u).
	type dirKey struct{ U, V int }
	edgeMap := make(map[dirKey]int, len(tris)*3)

	for ti, t := range tris {
		if t.A < 0 || t.A >= len(sites) || t.B < 0 || t.B >= len(sites) || t.C < 0 || t.C >= len(sites) {
			return nil, fmt.Errorf("triangle %d has out-of-range indices: %+v", ti, t)
		}

		e0 := len(m.HalfEdges)
		m.TriEdge0[ti] = e0

		// Create edges in cycle A->B, B->C, C->A
		m.HalfEdges = append(m.HalfEdges,
			HalfEdge{Origin: t.A, Tri: ti, Next: e0 + 1, Twin: -1},
			HalfEdge{Origin: t.B, Tri: ti, Next: e0 + 2, Twin: -1},
			HalfEdge{Origin: t.C, Tri: ti, Next: e0 + 0, Twin: -1},
		)

		// Fill Prev and EdgeDest cache
		m.HalfEdges[e0+0].Prev = e0 + 2
		m.HalfEdges[e0+1].Prev = e0 + 0
		m.HalfEdges[e0+2].Prev = e0 + 1

		m.HalfEdges[e0+0].EdgeDest = m.HalfEdges[m.HalfEdges[e0+0].Next].Origin
		m.HalfEdges[e0+1].EdgeDest = m.HalfEdges[m.HalfEdges[e0+1].Next].Origin
		m.HalfEdges[e0+2].EdgeDest = m.HalfEdges[m.HalfEdges[e0+2].Next].Origin

		// Remember one outgoing edge per site (arbitrary).
		if m.SiteEdge[t.A] == -1 {
			m.SiteEdge[t.A] = e0 + 0
		}
		if m.SiteEdge[t.B] == -1 {
			m.SiteEdge[t.B] = e0 + 1
		}
		if m.SiteEdge[t.C] == -1 {
			m.SiteEdge[t.C] = e0 + 2
		}

		// Register edges for twin linking
		for _, ei := range []int{e0, e0 + 1, e0 + 2} {
			u := m.HalfEdges[ei].Origin
			v := m.HalfEdges[ei].EdgeDest
			k := dirKey{U: u, V: v}
			if _, exists := edgeMap[k]; exists {
				// This usually indicates duplicate triangles or inconsistent input.
				// Not always fatal, but it will break traversal.
				return nil, fmt.Errorf("duplicate directed edge %d->%d detected at halfedge %d", u, v, ei)
			}
			edgeMap[k] = ei
		}
	}

	// Link twins: for each directed edge u->v, find v->u.
	for k, ei := range edgeMap {
		opp := dirKey{U: k.V, V: k.U}
		if tj, ok := edgeMap[opp]; ok {
			m.HalfEdges[ei].Twin = tj
		}
	}

	// Compute circumcenters (Voronoi vertices)
	for ti, t := range tris {
		a := sites[t.A].Pos
		b := sites[t.B].Pos
		c := sites[t.C].Pos
		cc, ok := circumcenter(a, b, c)
		if !ok {
			m.TriIsDegenerate[ti] = true
			// fallback: centroid to keep things stable
			m.TriCircumcenter[ti] = Vec2{(a.X + b.X + c.X) / 3.0, (a.Y + b.Y + c.Y) / 3.0}
		} else {
			m.TriCircumcenter[ti] = cc
		}
	}

	return m, nil
}

// circumcenter returns circumcenter of triangle (a,b,c).
// ok=false if triangle is degenerate (area near zero).
func circumcenter(a, b, c Vec2) (cc Vec2, ok bool) {
	// Based on perpendicular bisector intersection formula.
	d := 2 * (a.X*(b.Y-c.Y) + b.X*(c.Y-a.Y) + c.X*(a.Y-b.Y))
	const eps = 1e-12
	if math.Abs(d) < eps {
		return Vec2{}, false
	}

	a2 := a.X*a.X + a.Y*a.Y
	b2 := b.X*b.X + b.Y*b.Y
	c2 := c.X*c.X + c.Y*c.Y

	ux := (a2*(b.Y-c.Y) + b2*(c.Y-a.Y) + c2*(a.Y-b.Y)) / d
	uy := (a2*(c.X-b.X) + b2*(a.X-c.X) + c2*(b.X-a.X)) / d
	return Vec2{ux, uy}, true
}

type VoronoiCell struct {
	Site     int   // site index
	Tris     []int // triangles around the site, in order
	Vertices []Vec2
	Closed   bool // false if boundary / infinite cell encountered
}

// VoronoiCellForSite walks the 1-ring around a site and returns
// the ordered circumcenters of adjacent triangles (the polygon vertices).
//
// If the site touches the boundary of the triangulation, there will be a missing twin,
// and the cell is “open” (infinite in true Voronoi). We return Closed=false.
func (m *DelaunayMesh) VoronoiCellForSite(site int) VoronoiCell {
	start := -1
	if site >= 0 && site < len(m.SiteEdge) {
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

func (m *DelaunayMesh) findAnyOutgoing(site int) int {
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
