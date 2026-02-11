package core

type SiteIndex int
type SideIndex int //

// ChunkCoord identifies a chunk by integer grid coordinates.
type ChunkCoord struct {
	X, Z int
}

const (
	NORTH_EDGE SideIndex = iota
	SOUTH_EDGE
	EAST_EDGE
	WEST_EDGE
)

// Site represents a terrain vertex with position and height.
type Site struct {
	Pos    Vec2    // 2D position (X, Z plane)
	Height float64 // elevation (Y axis)
}

// Triangle references 3 site indices in CCW order (recommended).
type Triangle struct{ A, B, C SiteIndex }

// HalfEdge is a directed edge: Origin -> Dest,
// where Dest is the Origin of Next.
type HalfEdge struct {
	Origin   SiteIndex // site index
	Tri      int       // triangle index
	Next     int       // half-edge index (within triangle cycle)
	Twin     int       // opposite direction half-edge index, or -1 if boundary
	Prev     int       // (optional) convenient; set during build. Otherwise unused (currently)
	EdgeDest SiteIndex // cached dest site (optional, set during build)
}

// DelaunayMesh holds topology + geometry.
// It’s “Delaunay” if the triangles were produced by a Delaunay triangulation.
type DelaunayMesh struct {
	Sites     []Site     // Vertices
	Tris      []Triangle // Faces
	HalfEdges []HalfEdge // Edges

	// Face normals per triangle (parallel to Tris)
	FaceNormals []Vec3

	// For each triangle, index of its first half-edge (3 consecutive edges)
	TriEdge0 []int

	// For each site, one outgoing half-edge from that site (or -1)
	SiteEdge []int

	// Voronoi vertices: one per triangle
	TriCircumcenter []Vec2
	TriIsDegenerate []bool
}

type VoronoiCell struct {
	Site     SiteIndex // site index
	Tris     []int     // triangles around the site, in order
	Vertices []Vec2
	Closed   bool // false if boundary / infinite cell encountered
}
