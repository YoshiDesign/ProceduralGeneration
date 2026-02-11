package core

import (
	"fmt"
	"math"

	"github.com/hajimehoshi/ebiten/v2"
)

// TerrainChunk represents a generated terrain chunk with mesh data.
type TerrainChunk struct {
	Coord ChunkCoord
	Cfg   ChunkConfig
	Seed  int64

	// Core bounds (what this chunk "owns")
	MinX, MinZ float64
	MaxX, MaxZ float64

	// The Delaunay mesh (includes halo points for boundary continuity)
	Mesh *DelaunayMesh

	// Height values per site (parallel to Mesh.Sites)
	Heights []float64

	// Watershed divide weights per site (parallel to Mesh.Sites)
	// 0.0 = normal terrain, 1.0 = full watershed divide
	// Used by vertex shader to elevate terrain at lake boundaries
	WatershedWeights []float64
	
	// Spatial index for fast point location
	Spatial *SpatialGrid

	// Which sites are in the core region (not halo)
	CoreSiteIndices []SiteIndex

	// Hydrology data
	Hydro *ChunkHydroData

	// Pre-computed render data for batched drawing
	RenderVertices []ebiten.Vertex
	RenderIndices  []uint16

	// Pre-computed hydrology render data
	LakeVertices  []ebiten.Vertex
	LakeIndices   []uint16
	OceanVertices []ebiten.Vertex
	OceanIndices  []uint16

	// Erosion
	ErosionHeightDeltas []float64
}

// ChunkConfig holds parameters for terrain chunk generation.
type ChunkConfig struct {
	ChunkSize    float64 // World units per chunk side (e.g., 256.0)
	MinPointDist float64 // Minimum distance between blue noise points
	HaloWidth    float64 // Boundary overlap region width (typically = MinPointDist)
	WorldSeed    int64   // Global world seed
	ChunksX      int     // Number of chunks along the X axis
	ChunksZ      int     // Number of chunks along the Z axis
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
		if int(t.A) < 0 || int(t.A) >= len(sites) || int(t.B) < 0 || int(t.B) >= len(sites) || int(t.C) < 0 || int(t.C) >= len(sites) {
			return nil, fmt.Errorf("triangle %d has out-of-range indices: %+v", ti, t)
		}

		e0 := len(m.HalfEdges)
		m.TriEdge0[ti] = e0

		// Create edges in cycle A->B, B->C, C->A
		m.HalfEdges = append(m.HalfEdges,
			HalfEdge{Origin: SiteIndex(t.A), Tri: ti, Next: e0 + 1, Twin: -1},
			HalfEdge{Origin: SiteIndex(t.B), Tri: ti, Next: e0 + 2, Twin: -1},
			HalfEdge{Origin: SiteIndex(t.C), Tri: ti, Next: e0 + 0, Twin: -1},
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
			m.TriCircumcenter[ti] = Vec2{X: (a.X + b.X + c.X) / 3.0, Y: (a.Y + b.Y + c.Y) / 3.0}
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
	return Vec2{X: ux, Y: uy}, true
}


// adjTri is an adjacency-aware triangle for the optimized Bowyer-Watson algorithm.
// Vertices a, b, c are in CCW order.
// n0 is the neighbor opposite vertex a (sharing edge b-c), etc.
// -1 means no neighbor (boundary edge).
type adjTri struct {
	a, b, c    int  // vertex indices (CCW)
	n0, n1, n2 int  // neighbor opposite to vertex a, b, c
	alive      bool // false = deleted (lazy deletion)
}

// edgeKey is a canonical representation of an edge (smaller index first).
type edgeKey struct{ a, b int }

// makeEdgeKey creates a canonical edge key.
func makeEdgeKey(a, b int) edgeKey {
	if a > b {
		a, b = b, a
	}
	return edgeKey{a, b}
}

// triEdgeRef stores a reference to a triangle and which edge slot.
type triEdgeRef struct {
	triIdx   int
	edgeSlot int // 0 = edge b-c (opposite a), 1 = edge c-a (opposite b), 2 = edge a-b (opposite c)
}

// orientation2D returns positive if p is left of line a->b, negative if right, zero if collinear.
func orientation2D(a, b, p Vec2) float64 {
	return (b.X-a.X)*(p.Y-a.Y) - (b.Y-a.Y)*(p.X-a.X)
}

// pointInTriangle checks if point p is inside triangle (a, b, c) using orientation tests.
// Returns true if inside or on boundary.
func pointInTriangle(a, b, c, p Vec2) bool {
	o1 := orientation2D(a, b, p)
	o2 := orientation2D(b, c, p)
	o3 := orientation2D(c, a, p)
	// All same sign (or zero) means inside
	return (o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0)
}

// WalkStats tracks walking performance metrics for debugging.
type WalkStats struct {
	TotalWalks     int
	TotalSteps     int
	MaxSteps       int
	FallbackCount  int
	TotalBadTris   int
}

// Global stats for debugging (reset before each Triangulate call)
var walkStats WalkStats

// walkToPoint finds a triangle containing point p using adjacency walking.
// Returns the index of a triangle whose circumcircle contains p (a "bad" triangle).
// startTri should be a valid, alive triangle index.
func walkToPoint(triangles []adjTri, pts []Vec2, startTri int, p Vec2) int {
	current := startTri
	maxSteps := len(triangles) + 100 // Safety limit
	stepsTaken := 0

	for step := 0; step < maxSteps; step++ {
		stepsTaken = step
		t := triangles[current]
		if !t.alive {
			// Find next alive triangle
			for i := 0; i < len(triangles); i++ {
				if triangles[i].alive {
					current = i
					break
				}
			}
			continue
		}

		a, b, c := pts[t.a], pts[t.b], pts[t.c]

		// Check if point is inside this triangle
		if pointInTriangle(a, b, c, p) {
			walkStats.TotalWalks++
			walkStats.TotalSteps += stepsTaken
			if stepsTaken > walkStats.MaxSteps {
				walkStats.MaxSteps = stepsTaken
			}
			return current
		}

		// Find which edge to cross - walk toward p
		// Check each edge and cross if p is on the other side
		o0 := orientation2D(b, c, p) // Edge opposite to a
		o1 := orientation2D(c, a, p) // Edge opposite to b
		o2 := orientation2D(a, b, p) // Edge opposite to c

		// Determine triangle winding to know which direction is "outside"
		triOrient := orientation2D(a, b, c)

		// Cross the edge where p is on the opposite side
		if triOrient > 0 {
			// CCW triangle: negative orientation means p is outside that edge
			if o0 < 0 && t.n0 >= 0 && triangles[t.n0].alive {
				current = t.n0
				continue
			}
			if o1 < 0 && t.n1 >= 0 && triangles[t.n1].alive {
				current = t.n1
				continue
			}
			if o2 < 0 && t.n2 >= 0 && triangles[t.n2].alive {
				current = t.n2
				continue
			}
		} else {
			// CW triangle: positive orientation means p is outside that edge
			if o0 > 0 && t.n0 >= 0 && triangles[t.n0].alive {
				current = t.n0
				continue
			}
			if o1 > 0 && t.n1 >= 0 && triangles[t.n1].alive {
				current = t.n1
				continue
			}
			if o2 > 0 && t.n2 >= 0 && triangles[t.n2].alive {
				current = t.n2
				continue
			}
		}

		// No valid neighbor to cross to - this triangle is our best bet
		walkStats.TotalWalks++
		walkStats.TotalSteps += stepsTaken
		if stepsTaken > walkStats.MaxSteps {
			walkStats.MaxSteps = stepsTaken
		}
		return current
	}

	// Hit max steps - should rarely happen
	walkStats.TotalWalks++
	walkStats.TotalSteps += maxSteps
	if maxSteps > walkStats.MaxSteps {
		walkStats.MaxSteps = maxSteps
	}
	return current
}

// floodFillBadTris finds all triangles whose circumcircles contain point p,
// starting from a known bad triangle and flooding via adjacency.
func floodFillBadTris(triangles []adjTri, pts []Vec2, startTri int, p Vec2) []int {
	badTris := make([]int, 0, 8)
	visited := make(map[int]bool, 16)

	var flood func(ti int)
	flood = func(ti int) {
		if ti < 0 || visited[ti] {
			return
		}
		visited[ti] = true

		t := triangles[ti]
		if !t.alive {
			return
		}

		if inCircumcircle(pts[t.a], pts[t.b], pts[t.c], p) {
			badTris = append(badTris, ti)
			// Flood to neighbors
			flood(t.n0)
			flood(t.n1)
			flood(t.n2)
		}
	}

	flood(startTri)
	return badTris
}

// getTriEdges returns the three edges of a triangle with their edge slots.
// Edge slot 0 = b-c (opposite a), 1 = c-a (opposite b), 2 = a-b (opposite c)
func getTriEdges(t adjTri) [3]struct {
	key  edgeKey
	slot int
} {
	return [3]struct {
		key  edgeKey
		slot int
	}{
		{makeEdgeKey(t.b, t.c), 0},
		{makeEdgeKey(t.c, t.a), 1},
		{makeEdgeKey(t.a, t.b), 2},
	}
}

// getNeighborSlot returns a pointer to the neighbor slot for a given edge slot.
func getNeighborSlot(t *adjTri, slot int) *int {
	switch slot {
	case 0:
		return &t.n0
	case 1:
		return &t.n1
	case 2:
		return &t.n2
	}
	return nil
}

// Triangulate performs Delaunay triangulation on a set of 2D points.
// Uses the Bowyer-Watson algorithm with walking point location for O(n√n) complexity.
func Triangulate(points []Vec2) []Triangle {
	if len(points) < 3 {
		return nil
	}

	// Reset walk stats for this triangulation
	walkStats = WalkStats{}

	// Find bounding box
	minX, maxX := points[0].X, points[0].X
	minZ, maxZ := points[0].Y, points[0].Y
	for _, p := range points {
		if p.X < minX {
			minX = p.X
		}
		if p.X > maxX {
			maxX = p.X
		}
		if p.Y < minZ {
			minZ = p.Y
		}
		if p.Y > maxZ {
			maxZ = p.Y
		}
	}

	// fmt.Println("MinX:", minX, "MaxX:", maxX)
	// fmt.Println("MinZ:", minZ, "MaZZ:", maxZ)

	// Create super-triangle that encompasses all points
	dx := maxX - minX
	dz := maxZ - minZ
	deltaMax := dx
	if dz > dx {
		deltaMax = dz
	}
	midX := (minX + maxX) / 2
	midZ := (minZ + maxZ) / 2

	// Super-triangle vertices (large enough to contain all points)
	superA := Vec2{X: midX - 20*deltaMax, Y: midZ - deltaMax}
	superB := Vec2{X: midX, Y: midZ + 20*deltaMax}
	superC := Vec2{X: midX + 20*deltaMax, Y: midZ - deltaMax}

	// Working list of points including super-triangle
	allPts := make([]Vec2, len(points)+3)
	copy(allPts, points)
	allPts[len(points)] = superA
	allPts[len(points)+1] = superB
	allPts[len(points)+2] = superC

	superI := [3]int{len(points), len(points) + 1, len(points) + 2}

	// Initialize triangles with adjacency
	triangles := make([]adjTri, 1, len(points)*2+1)
	triangles[0] = adjTri{
		a: superI[0], b: superI[1], c: superI[2],
		n0: -1, n1: -1, n2: -1,
		alive: true,
	}

	// Edge-to-triangle map for O(1) adjacency lookups
	edgeToTri := make(map[edgeKey]triEdgeRef, len(points)*3)

	// Register initial triangle's edges
	for _, e := range getTriEdges(triangles[0]) {
		edgeToTri[e.key] = triEdgeRef{triIdx: 0, edgeSlot: e.slot}
	}

	// Track a known alive triangle for walking start point
	lastInsertedTri := 0

	// Insert points one at a time
	for pi := 0; pi < len(points); pi++ {
		p := allPts[pi]

		// Walk to find a triangle containing p
		startTri := walkToPoint(triangles, allPts, lastInsertedTri, p)

		// Flood-fill to find all bad triangles
		badTris := floodFillBadTris(triangles, allPts, startTri, p)

		if len(badTris) == 0 {
			// Fallback: linear scan (shouldn't happen normally)
			walkStats.FallbackCount++
			for ti := range triangles {
				t := triangles[ti]
				if t.alive && inCircumcircle(allPts[t.a], allPts[t.b], allPts[t.c], p) {
					badTris = append(badTris, ti)
				}
			}
		}

		walkStats.TotalBadTris += len(badTris)

		// Find boundary polygon edges (edges that appear exactly once)
		type polyEdge struct {
			a, b     int // Original edge direction for new triangle winding
			neighbor int // The triangle on the other side (outside the cavity)
		}
		edgeCount := make(map[edgeKey]int, len(badTris)*3)
		edgeInfo := make(map[edgeKey]polyEdge, len(badTris)*3)

		for _, ti := range badTris {
			t := triangles[ti]
			// Edge 0: b-c, neighbor n0
			// Edge 1: c-a, neighbor n1
			// Edge 2: a-b, neighbor n2
			edges := []struct {
				a, b     int
				neighbor int
			}{
				{t.b, t.c, t.n0},
				{t.c, t.a, t.n1},
				{t.a, t.b, t.n2},
			}

			for _, e := range edges {
				key := makeEdgeKey(e.a, e.b)
				edgeCount[key]++
				edgeInfo[key] = polyEdge{a: e.a, b: e.b, neighbor: e.neighbor}
			}
		}

		// Remove edges from edgeToTri for bad triangles
		for _, ti := range badTris {
			t := triangles[ti]
			for _, e := range getTriEdges(t) {
				delete(edgeToTri, e.key)
			}
			// Mark as dead
			triangles[ti].alive = false
		}

		// Collect boundary edges (those appearing exactly once)
		polygon := make([]polyEdge, 0, len(badTris)+2)
		for key, count := range edgeCount {
			if count == 1 {
				polygon = append(polygon, edgeInfo[key])
			}
		}

		// Create new triangles from polygon edges to the new point
		newTriIndices := make([]int, 0, len(polygon))
		for _, e := range polygon {
			newTri := adjTri{
				a: e.a, b: e.b, c: pi,
				n0: -1, n1: -1, n2: e.neighbor, // n2 is opposite vertex c, which is edge a-b
				alive: true,
			}

			triIdx := len(triangles)
			triangles = append(triangles, newTri)
			newTriIndices = append(newTriIndices, triIdx)

			// Update the neighbor's adjacency to point back to us
			if e.neighbor >= 0 && triangles[e.neighbor].alive {
				neighborTri := &triangles[e.neighbor]
				edgeKeyAB := makeEdgeKey(e.a, e.b)
				// Find which edge slot in neighbor matches this edge
				for _, ne := range getTriEdges(*neighborTri) {
					if ne.key == edgeKeyAB {
						*getNeighborSlot(neighborTri, ne.slot) = triIdx
						break
					}
				}
			}

			// Register edges in edgeToTri
			for _, edge := range getTriEdges(newTri) {
				edgeToTri[edge.key] = triEdgeRef{triIdx: triIdx, edgeSlot: edge.slot}
			}
		}

		// Link new triangles to each other via shared edges
		for i, ti := range newTriIndices {
			t := &triangles[ti]
			// Edge 0: b-c (opposite a) -> b is polygon edge's a, c is pi
			// Edge 1: c-a (opposite b) -> c is pi, a is polygon edge's a
			// These edges are shared with other new triangles

			// Edge with pi as one vertex - find sibling triangle sharing this edge
			// Edge 0: vertices b, c (where c = pi)
			key0 := makeEdgeKey(t.b, t.c)
			// Edge 1: vertices c, a (where c = pi)
			key1 := makeEdgeKey(t.c, t.a)

			for j, tj := range newTriIndices {
				if i == j {
					continue
				}
				ot := &triangles[tj]
				// Check if they share edge 0
				if t.n0 < 0 {
					for _, oe := range getTriEdges(*ot) {
						if oe.key == key0 {
							t.n0 = tj
							*getNeighborSlot(ot, oe.slot) = ti
							break
						}
					}
				}
				// Check if they share edge 1
				if t.n1 < 0 {
					for _, oe := range getTriEdges(*ot) {
						if oe.key == key1 {
							t.n1 = tj
							*getNeighborSlot(ot, oe.slot) = ti
							break
						}
					}
				}
			}
		}

		// Update walking start point
		if len(newTriIndices) > 0 {
			lastInsertedTri = newTriIndices[0]
		}
	}

	// Collect final triangles (excluding super-triangle vertices)
	result := make([]Triangle, 0, len(triangles))
	for _, t := range triangles {
		if !t.alive {
			continue
		}
		usesSuperVertex := false
		for _, si := range superI {
			if t.a == si || t.b == si || t.c == si {
				usesSuperVertex = true
				break
			}
		}
		if !usesSuperVertex {
			// Ensure CCW winding
			result = append(result, ensureCCW(allPts, SiteIndex(t.a), SiteIndex(t.b), SiteIndex(t.c)))
		}
	}

	// // Log walking stats for performance validation
	// if walkStats.TotalWalks > 0 {
	// 	avgSteps := float64(walkStats.TotalSteps) / float64(walkStats.TotalWalks)
	// 	avgBadTris := float64(walkStats.TotalBadTris) / float64(walkStats.TotalWalks)
	// 	fmt.Printf("[Triangulate] n=%d: walks=%d, avgSteps=%.2f, maxSteps=%d, fallbacks=%d, avgBadTris=%.2f\n",
	// 		len(points), walkStats.TotalWalks, avgSteps, walkStats.MaxSteps, walkStats.FallbackCount, avgBadTris)
	// }

	return result
}

// inCircumcircle returns true if point p is inside the circumcircle of triangle (a, b, c).
func inCircumcircle(a, b, c, p Vec2) bool {
	// Using the determinant method
	ax, ay := a.X-p.X, a.Y-p.Y
	bx, by := b.X-p.X, b.Y-p.Y
	cx, cy := c.X-p.X, c.Y-p.Y

	det := (ax*ax+ay*ay)*(bx*cy-cx*by) -
		(bx*bx+by*by)*(ax*cy-cx*ay) +
		(cx*cx+cy*cy)*(ax*by-bx*ay)

	// For CCW triangles, det > 0 means p is inside
	// Handle both windings by checking the sign of the triangle area
	area := (b.X-a.X)*(c.Y-a.Y) - (c.X-a.X)*(b.Y-a.Y)
	if area < 0 {
		return det < 0
	}
	return det > 0
}

// ensureCCW returns a Triangle with vertices in counter-clockwise order.
func ensureCCW(pts []Vec2, a, b, c SiteIndex) Triangle {
	// Cross product of (b-a) × (c-a)
	cross := (pts[b].X-pts[a].X)*(pts[c].Y-pts[a].Y) - (pts[c].X-pts[a].X)*(pts[b].Y-pts[a].Y)
	if cross < 0 {
		return Triangle{A: a, C: b, B: c} // Swap b and c
	}
	return Triangle{A: a, B: b, C: c}
}

