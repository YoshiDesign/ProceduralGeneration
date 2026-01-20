package hydro

import (
	"math"
	"math/rand"
	"procedural_generation/terrain_generation/duals2/core"
)

// HydroManager manages hydrology across all chunks.
type HydroManager struct {
	Cfg HydroConfig
	Rng *rand.Rand

	// Global water bodies
	Oceans []core.OceanRegion
	Rivers map[int]*core.River // RiverID -> River
	Lakes  map[int]*core.Lake  // LakeID -> Lake

	// Cross-chunk river propagation
	// Key: destination ChunkCoord, Value: rivers entering that chunk
	RiverPropagation map[core.ChunkCoord][]core.RiverInterPoint

	// Ocean region ownership
	// Key: ChunkCoord, Value: OceanID (-1 if not ocean)
	OceanChunks map[core.ChunkCoord]int

	// Tracks visited sites during river tracing within a chunk.
	// Reset before processing each chunk's hydrology.
	VisitedSites map[int]bool

	// ID counters
	NextRiverID int
	NextLakeID  int
	NextOceanID int
}


// HydroConfig holds parameters for hydrology generation.
type HydroConfig struct {
	// Sea level
	SeaLevel float64 // Global sea level elevation (e.g., 0.0)

	// River parameters
	RiverMinSlope       float64 // Minimum slope to continue river flow
	RiverWidthBase      float64 // Starting width at source
	RiverWidthGrowth    float64 // Width increase per unit distance traveled
	RiverDepthBase      float64 // Starting depth at source
	RiverDepthGrowth    float64 // Depth increase per unit distance traveled
	RiverSourceMinElev  float64 // Minimum elevation for river sources
	RiverSourceMaxElev  float64 // TODO Maximum elevation for river sources
	RiverSourceProbability float64 // Probability a valid source spawns a river (0-1)

	// Lake parameters
	LakeMinDepth     float64 // Minimum depression depth to form a lake
	LakeMinArea      int     // Minimum number of vertices to form a lake
	LakeProbability  float64 // Probability an eligible depression becomes a lake (0-1)

	// Ocean parameters
	OceanMaxChunksX  int     // Maximum X span in chunks
	OceanMaxChunksZ  int     // Maximum Z span in chunks
	OceanProbability float64 // Probability of ocean region spawning (per seed point)

	// Flow direction control
	EquatorZ float64 // Z-coordinate of the equator for flow bias
}

// DefaultHydroConfig returns sensible defaults for hydrology generation.
func DefaultHydroConfig() HydroConfig {
	return HydroConfig{
		SeaLevel: 0.0,

		RiverMinSlope:          0.001,  // Very gentle minimum slope
		RiverWidthBase:         2.0,    // 2 units wide at source
		RiverWidthGrowth:       0.01,   // Grows 1% per unit distance
		RiverDepthBase:         0.5,    // 0.5 units deep at source
		RiverDepthGrowth:       0.005,  // Grows slower than width
		RiverSourceMinElev:     1.0,   // Sources above 10 units elevation (lowered for demo compatibility)
		RiverSourceMaxElev:     100.0, // Sources below 100 units elevation
		RiverSourceProbability: 0.1,    // 30% of valid sources spawn rivers

		LakeMinDepth:    2.0, // Depression must be at least 2 units deep
		LakeMinArea:     3,   // At least 3 vertices
		LakeProbability: 0.0, // 50% of eligible depressions become lakes

		OceanMaxChunksX:  8,   // Ocean can span up to 8 chunks in X
		OceanMaxChunksZ:  8,   // Ocean can span up to 8 chunks in Z
		OceanProbability: 0.0, // 10% base probability

		EquatorZ: 0.0, // Equator at Z=0
	}
}


// NewHydroManager creates a new hydrology manager.
	func NewHydroManager(cfg HydroConfig, seed int64) *HydroManager {
	return &HydroManager{
		Cfg:              cfg,
		Rng:              rand.New(rand.NewSource(seed)),
		Rivers:           make(map[int]*core.River),
		Lakes:            make(map[int]*core.Lake),
		RiverPropagation: make(map[core.ChunkCoord][]core.RiverInterPoint),
		OceanChunks:      make(map[core.ChunkCoord]int),
		VisitedSites:     make(map[int]bool),
	}
}


// FlowBias returns the preferred flow direction based on position relative to equator.
// North of equator (Z > EquatorZ): bias toward south (negative Z).
// South of equator (Z < EquatorZ): bias toward north (positive Z).
func (cfg *HydroConfig) FlowBias(worldZ float64) core.Vec2 {

	if worldZ > cfg.EquatorZ {
		return core.Vec2{X: 0, Y: -1} // Flow south
	}
	return core.Vec2{X: 0, Y: 1} // Flow north
}


// ResetVisitedSites clears the visited sites map for a new chunk's hydrology pass.
func (hm *HydroManager) ResetVisitedSites() {
	hm.VisitedSites = make(map[int]bool)
}

// IsSourceVisited returns true if the given site has already been visited by a river trace.
func (hm *HydroManager) IsSourceVisited(site int) bool {
	return hm.VisitedSites != nil && hm.VisitedSites[site]
}

// -------------------------------------------------------------------
// River Tracing Algorithm
// -------------------------------------------------------------------

// getLowerNeighbors returns site indices of neighbors with lower elevation.
func getLowerNeighbors(mesh *core.DelaunayMesh, cfg HydroConfig, heights []float64, site int) []int {
	neighbors := make([]int, 0, 8)
	currentHeight := heights[site]

	// Use half-edge structure to find neighbors
	startEdge := mesh.SiteEdge[site]
	if startEdge == -1 {
		return neighbors
	}

	edge := startEdge
	visited := make(map[int]bool)

	for {
		if visited[edge] {
			break
		}
		visited[edge] = true

		// Get destination of this edge
		dest := mesh.HalfEdges[edge].EdgeDest
		if heights[dest] < currentHeight-cfg.RiverMinSlope {
			neighbors = append(neighbors, dest)
		}

		// Move to next edge around the vertex
		twin := mesh.HalfEdges[edge].Twin
		if twin == -1 {
			break
		}
		edge = mesh.HalfEdges[twin].Next
		if edge == startEdge {
			break
		}
	}

	return neighbors
}

// selectNextVertex chooses the next vertex for river flow.
// Prefers steeper slopes but applies flow bias for tie-breaking.
// TODO: This wont matter. Vulkan can carve f*cking valleys.
//		 Don't bias based on slope. This will be for the tectonics layer in the near future.
func (hm *HydroManager) selectNextVertex(
	mesh *core.DelaunayMesh,
	heights []float64,
	current int,
	candidates []int,
	currentHeight float64,
	bias core.Vec2,
) int {
	if len(candidates) == 1 {
		return candidates[0]
	}

	currentPos := mesh.Sites[current].Pos
	bestScore := math.Inf(-1)
	bestCandidate := candidates[0]

	for _, c := range candidates {
		pos := mesh.Sites[c].Pos
		h := heights[c]

		// Calculate slope (height drop per unit distance)
		dist := pos.Sub(currentPos).Len()
		if dist < 1e-9 {
			continue
		}
		slope := (currentHeight - h) / dist

		// Calculate alignment with flow bias
		dir := pos.Sub(currentPos).Mul(1.0 / dist)
		biasAlignment := dir.Dot(bias) // -1 to 1

		// Score: primarily slope, with bias as tie-breaker
		// Slope is weighted heavily (10x), bias is secondary
		score := slope*10.0 + biasAlignment*0.5

		// Add small random noise for meandering (±5%)
		score *= 1.0 + (hm.Rng.Float64()-0.5)*0.1

		if score > bestScore {
			bestScore = score
			bestCandidate = c
		}
	}

	return bestCandidate
}

// Prefer the vertex most aligned with the flow bias
// NOTE: This will ignore slope, implying the vertex shader will reduce height to the selected vertex-Y, creating a valley.
func (hm *HydroManager) selectNextVertexV2(
	mesh *core.DelaunayMesh,
	current int,
	candidates []int,
	bias core.Vec2,
) int {

	if len(candidates) == 1 {
		return candidates[0]
	}

	currentPos := mesh.Sites[current].Pos
	bestScore := math.Inf(-1)
	bestCandidate := candidates[0]

	for _, c := range candidates {
		pos := mesh.Sites[c].Pos

		// Too close to constitute a river segment, I guess!
		dist := pos.Sub(currentPos).Len()
		if dist < 1e-9 {
			continue
		}

		// Calculate alignment with flow bias
		dir := pos.Sub(currentPos).Mul(1.0 / dist)
		biasAlignment := dir.Dot(bias) // -1 to 1

		// If we weight by slope we can cut mountains in half... make this a V3.

		score := biasAlignment*0.5

		// Add small random noise for meandering (±5%)
		score *= 1.0 + (hm.Rng.Float64()-0.5)*0.1

		if score > bestScore {
			bestScore = score
			bestCandidate = c
		}
	}

	return bestCandidate

}
