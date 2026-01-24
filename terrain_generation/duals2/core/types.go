package core

import (
	"math"
)

type SiteIndex int
type SideIndex int

const (
    NORTH_EDGE SideIndex = iota
    SOUTH_EDGE
    EAST_EDGE
    WEST_EDGE
)

// NoiseParams holds configurable parameters for fractal noise terrain generation.
type NoiseParams struct {
	Octaves     int
	Frequency   float64
	Amplitude   float64
	Persistence float64
	Lacunarity  float64
}

type Vec2 struct{ X, Y float64 }

func (a Vec2) Add(b Vec2) Vec2    { return Vec2{a.X + b.X, a.Y + b.Y} }
func (a Vec2) Sub(b Vec2) Vec2    { return Vec2{a.X - b.X, a.Y - b.Y} }
func (a Vec2) Mul(s float64) Vec2 { return Vec2{a.X * s, a.Y * s} }
func (a Vec2) Dot(b Vec2) float64 { return a.X*b.X + a.Y*b.Y }
func (a Vec2) Len2() float64      { return a.Dot(a) }
func (a Vec2) Len() float64       { return math.Sqrt(a.Len2()) }
func (a Vec2) Eq(b Vec2) bool {
	return a.X == b.X && a.Y == b.Y
}

// Normalize returns a unit vector (or zero vector if length is zero).
func (v Vec2) Normalize() Vec2 {
	l := v.Len()
	if l < 1e-12 {
		return Vec2{}
	}
	return v.Mul(1.0 / l)
}

// Vec3 represents a 3D point/vector (X = east, Y = up, Z = north).
type Vec3 struct{ X, Y, Z float64 }

func (a Vec3) Add(b Vec3) Vec3    { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }
func (a Vec3) Sub(b Vec3) Vec3    { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }
func (a Vec3) Mul(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }
func (a Vec3) Dot(b Vec3) float64 { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }
func (a Vec3) Len2() float64      { return a.Dot(a) }
func (a Vec3) Len() float64       { return math.Sqrt(a.Len2()) }
func (a Vec3) Eq(b Vec3) bool {
	return a.X == b.X && a.Y == b.Y && a.Z == b.Z
}

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
	Origin   SiteIndex // site index
	Tri      int // triangle index
	Next     int // half-edge index (within triangle cycle)
	Twin     int // opposite direction half-edge index, or -1 if boundary
	Prev     int // (optional) convenient; set during build. Otherwise unused (currently)
	EdgeDest SiteIndex // cached dest site (optional, set during build)
}

// DelaunayMesh holds topology + geometry.
// It’s “Delaunay” if the triangles were produced by a Delaunay triangulation.
type DelaunayMesh struct {
	Sites     []Site	// Vertices
	Tris      []Triangle // Faces
	HalfEdges []HalfEdge // Edges

	// For each triangle, index of its first half-edge (3 consecutive edges)
	TriEdge0 []int

	// For each site, one outgoing half-edge from that site (or -1)
	SiteEdge []int

	// Voronoi vertices: one per triangle
	TriCircumcenter []Vec2
	TriIsDegenerate []bool
}

// ChunkCoord identifies a chunk by integer grid coordinates.
type ChunkCoord struct {
	X, Z int
}

type VoronoiCell struct {
	Site     SiteIndex   // site index
	Tris     []int // triangles around the site, in order
	Vertices []Vec2
	Closed   bool // false if boundary / infinite cell encountered
}
