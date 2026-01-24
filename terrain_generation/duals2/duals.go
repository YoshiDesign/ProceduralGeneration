package duals2

import (
	"fmt"
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// BuildHalfEdgeMesh constructs connectivity from Sites + Tris.
// Triangles must reference valid site indices.
// If triangles aren’t consistently wound, Voronoi polygons can come out flipped/odd.
func BuildHalfEdgeMesh(sites []core.Site, tris []core.Triangle) (*core.DelaunayMesh, error) {
	m := &core.DelaunayMesh{
		Sites:            sites,
		Tris:             tris,
		TriEdge0:         make([]int, len(tris)),
		SiteEdge:         make([]int, len(sites)),
		TriCircumcenter:  make([]core.Vec2, len(tris)),
		TriIsDegenerate:  make([]bool, len(tris)),
	}
	for i := range m.SiteEdge {
		m.SiteEdge[i] = -1
	}

	// Allocate 3 half-edges per triangle.
	m.HalfEdges = make([]core.HalfEdge, 0, len(tris)*3)

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
			core.HalfEdge{Origin: core.SiteIndex(t.A), Tri: ti, Next: e0 + 1, Twin: -1},
			core.HalfEdge{Origin: core.SiteIndex(t.B), Tri: ti, Next: e0 + 2, Twin: -1},
			core.HalfEdge{Origin: core.SiteIndex(t.C), Tri: ti, Next: e0 + 0, Twin: -1},
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
			k := dirKey{U: int(u), V: int(v)}
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
			m.TriCircumcenter[ti] = core.Vec2{X: (a.X + b.X + c.X) / 3.0, Y: (a.Y + b.Y + c.Y) / 3.0}
		} else {
			m.TriCircumcenter[ti] = cc
		}
	}

	return m, nil
}

// circumcenter returns circumcenter of triangle (a,b,c).
// ok=false if triangle is degenerate (area near zero).
func circumcenter(a, b, c core.Vec2) (cc core.Vec2, ok bool) {
	// Based on perpendicular bisector intersection formula.
	d := 2 * (a.X*(b.Y-c.Y) + b.X*(c.Y-a.Y) + c.X*(a.Y-b.Y))
	const eps = 1e-12
	if math.Abs(d) < eps {
		return core.Vec2{}, false
	}

	a2 := a.X*a.X + a.Y*a.Y
	b2 := b.X*b.X + b.Y*b.Y
	c2 := c.X*c.X + c.Y*c.Y

	ux := (a2*(b.Y-c.Y) + b2*(c.Y-a.Y) + c2*(a.Y-b.Y)) / d
	uy := (a2*(c.X-b.X) + b2*(a.X-c.X) + c2*(b.X-a.X)) / d
	return core.Vec2{X: ux, Y: uy}, true
}


