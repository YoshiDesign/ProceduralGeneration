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

	// Cross-chunk lake propagation
	// Key: destination ChunkCoord, Value: lake boundaries entering that chunk
	LakePropagation map[core.ChunkCoord][]core.LakeInterPoint

	// Track which chunks a lake spans (for water level updates)
	LakeChunks map[int][]core.ChunkCoord // LakeID -> list of chunks

	// Ocean region ownership
	// Key: ChunkCoord, Value: OceanID (-1 if not ocean)
	OceanChunks map[core.ChunkCoord]int

	// Tracks visited sites during river tracing within a chunk.
	// Reset before processing each chunk's hydrology.
	VisitedSites map[core.SiteIndex]bool

	// Tracks which sites are already part of a lake (for deduplication)
	// Key is (ChunkCoord, siteIndex) to ensure uniqueness across chunks
	LakeSites map[core.ChunkSiteKey]int

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
	RiverSourceMaxElev  float64 // Maximum elevation for river sources
	RiverSourceProbability float64 // Probability a valid source spawns a river (0-1)

	// Lake parameters
	LakeMinDepth                 float64 // Minimum depression depth to form a lake
	LakeMinArea                  int     // Minimum number of vertices to form a lake
	LakeProbability              float64 // Probability an eligible depression becomes a lake (0-1)
	LakeWaterLevelMergeThreshold float64 // If two lakes' water levels differ by less than this, merge them

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

		LakeMinDepth:                 2.0, // Depression must be at least 2 units deep
		LakeMinArea:                  3,   // At least 3 vertices
		LakeProbability:              0.0, // 50% of eligible depressions become lakes
		LakeWaterLevelMergeThreshold: 2.0, // Merge lakes if within 2 units of elevation

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
		LakePropagation:  make(map[core.ChunkCoord][]core.LakeInterPoint),
		LakeChunks:       make(map[int][]core.ChunkCoord),
		OceanChunks:      make(map[core.ChunkCoord]int),
		VisitedSites:     make(map[core.SiteIndex]bool),
		LakeSites:        make(map[core.ChunkSiteKey]int),
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
	hm.VisitedSites = make(map[core.SiteIndex]bool)
}

// IsSourceVisited returns true if the given site has already been visited by a river trace.
func (hm *HydroManager) IsSourceVisited(site core.SiteIndex) bool {
	return hm.VisitedSites != nil && hm.VisitedSites[site]
}

// -------------------------------------------------------------------
// River Tracing Algorithm
// -------------------------------------------------------------------

// selectNextVertex chooses the next vertex for river flow.
// Prefers steeper slopes but applies flow bias for tie-breaking.
// V1: Bias based on slope
func (hm *HydroManager) selectNextVertexV1(
	mesh *core.DelaunayMesh,
	heights []float64,
	current core.SiteIndex,
	candidates []core.SiteIndex,
	currentHeight float64,
	bias core.Vec2,
) core.SiteIndex {
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
// NOTE: V2 will ignore slope
func (hm *HydroManager) selectNextVertexV2(
	mesh *core.DelaunayMesh,
	current core.SiteIndex,
	candidates []core.SiteIndex,
	bias core.Vec2,
) core.SiteIndex {

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
