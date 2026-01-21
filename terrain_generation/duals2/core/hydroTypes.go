package core

// SourceType identifies how a river originates.
type SourceType int

const (
	SourceMountain SourceType = iota // High elevation in cold/wet biome
	SourceLake                       // Lake spillway
	SourceSpring                     // Underground spring (future)
)


// ChunkHydroData holds all hydrology information for a single chunk.
type ChunkHydroData struct {
	Rivers []RiverSegment
	Lakes  []Lake
	Ocean  OceanChunkData

	// Cross-chunk coordination for rivers
	RiverEntries []RiverInterPoint // Rivers entering from neighbors
	RiverExits   []RiverInterPoint // Rivers exiting to neighbors

	// Cross-chunk coordination for lakes
	LakeEntries []LakeBoundaryPoint // Lakes entering from neighbors
	LakeExits   []LakeBoundaryPoint // Lakes exiting to neighbors
}

// LakeBoundaryPoint stores information about a lake crossing a chunk boundary.
// This enables cross-chunk lake continuity.
type LakeBoundaryPoint struct {
	LakeID      int        // Global lake ID
	WaterLevel  float64    // Current water level (may be provisional)
	Position    Vec2       // World position at boundary
	SiteIndex   int        // Site index in source chunk's mesh
	EdgeIndex   int        // Which edge: 0=minX, 1=maxX, 2=minZ, 3=maxZ
	SourceChunk ChunkCoord // Chunk this boundary came from
}

// ChunkBounds represents the bounds of a chunk for boundary detection.
type ChunkBounds struct {
	MinX, MinZ, MaxX, MaxZ float64
}

// HydroBodyType identifies the type of water body.
type HydroBodyType int

const (
	HydroRiver HydroBodyType = iota
	HydroLake
	HydroOcean
)

// HydroBodyRef references a water body by type and ID.
type HydroBodyRef struct {
	Type HydroBodyType
	ID   int
}

// Lake represents a filled terrain depression.
type Lake struct {
	ID          int
	CenterPos   Vec2      // Approximate center of the lake
	SiteIndices []int     // Mesh vertex indices within the lake
	WaterLevel  float64   // Surface elevation (Y)
	MinDepth    float64   // Deepest point below water level
	Spillway    *Vec2     // Outlet point (nil if endorheic/closed basin)
	SpillwayDir Vec2      // Direction water flows out of spillway
	Inflows     []int     // River IDs flowing into this lake
	Outflow     *int      // River ID flowing out (nil if endorheic)
}

// OceanRegion represents an ocean spanning multiple chunks.
// Oceans are defined by their bounding chunk region and coastline.
type OceanRegion struct {
	ID        int
	MinChunkX int // Inclusive
	MinChunkZ int // Inclusive
	MaxChunkX int // Inclusive
	MaxChunkZ int // Inclusive
	SeaLevel  float64
}

// OceanChunkData holds per-chunk ocean information.
type OceanChunkData struct {
	IsOcean       bool      // True if this chunk is part of an ocean
	OceanID       int       // Which ocean region this belongs to
	CoastlinePts  []Vec2    // Coastline vertices within this chunk
	SubmergedTris []int     // Triangle indices below sea level
}

// RiverInterPoint stores information about a river when transitioning through a chunk boundary.
// This enables cross-chunk river continuity.
type RiverInterPoint struct {
	Position    Vec2    // World position at chunk boundary
	Width       float64 // River width at exit
	Depth       float64 // River depth at exit
	FlowDir     Vec2    // Flow direction at exit
	RiverID     int     // Source river ID (for linking segments)
	Distance    float64 // Total distance traveled from source
	ExitEdge    int     // Which edge: 0=minX, 1=maxX, 2=minZ, 3=maxZ
}


// RiverSegment represents a contiguous path of river vertices within a chunk.
// Vertices are ordered from upstream to downstream.
type RiverSegment struct {
	Vertices []Vec2    // Path points (source toward destination)
	Widths   []float64 // Width at each vertex (grows downstream)
	Depths   []float64 // Depth at each vertex
	FlowDir  Vec2      // Average normalized flow direction
}

// River represents a complete river with potential tributaries.
type River struct {
	ID          int
	Segments    []RiverSegment
	SourceType  SourceType
	SourcePos   Vec2         // World position where river originates
	FlowsInto   HydroBodyRef // What this river terminates into
	Tributaries []int        // River IDs that merge into this river
}
