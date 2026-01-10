package hexa

type HexType int
type CornerI int

const (
	PointyTop HexType = iota
	FlatTop
)

const (
	First CornerI = iota
	Second
	Third
	Fourth
	Fifth
	Sixth
	CornerCount // Sentinel
)

type Hexagon struct {
	inner_radius   float32 // Apothem
	Size_radius    float32 // Circumradius
	X, Y, Z        int     // basic cartesian positions
	layout         HexType
	CenterWorld    Vec3
	CornerWorld    []Vec3
	LocalLattice   []FineLatticePoint // Points mapped within a hex tile
	PointToLattice []FineGridIJ       // Tells you the (i, j) of the superimposed NxN lattice at a point inside of the hex tile
	LatticeToPoint []int              // Lookup Table: Either -1 (outside of hex tile) or the index of the corresponding FineLatticePoint
}

type FineGridIJ struct{ I, J int }

type FineLatticePoint struct {
	Pos     Vec3 // Tile local coordinates
	weights WeightField
}

type WeightField struct {
	temperature float32
	elevation   float32
}

// Flat x,z plane
type Vec2 struct {
	x float32
	z float32
}

type Vec3 struct {
	X float32
	Y float32
	Z float32
}
