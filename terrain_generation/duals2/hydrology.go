package duals2

import (
	"fmt"
	"math"
	"math/rand"
)

// SourceType identifies how a river originates.
type SourceType int

const (
	SourceMountain SourceType = iota // High elevation in cold/wet biome
	SourceLake                       // Lake spillway
	SourceSpring                     // Underground spring (future)
)

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

// FlowBias returns the preferred flow direction based on position relative to equator.
// North of equator (Z > EquatorZ): bias toward south (negative Z).
// South of equator (Z < EquatorZ): bias toward north (positive Z).
func (cfg *HydroConfig) FlowBias(worldZ float64) Vec2 {

	if worldZ > cfg.EquatorZ {
		return Vec2{0, -1} // Flow south
	}
	return Vec2{0, 1} // Flow north
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

// ChunkHydroData holds all hydrology information for a single chunk.
type ChunkHydroData struct {
	Rivers []RiverSegment
	Lakes  []Lake
	Ocean  OceanChunkData

	// Cross-chunk coordination
	RiverEntries []RiverInterPoint // Rivers entering from neighbors
	RiverExits   []RiverInterPoint // Rivers exiting to neighbors
}

// HydroManager manages hydrology across all chunks.
type HydroManager struct {
	cfg HydroConfig
	rng *rand.Rand

	// Global water bodies
	Oceans []OceanRegion
	Rivers map[int]*River // RiverID -> River
	Lakes  map[int]*Lake  // LakeID -> Lake

	// Cross-chunk river propagation
	// Key: destination ChunkCoord, Value: rivers entering that chunk
	riverPropagation map[ChunkCoord][]RiverInterPoint

	// Ocean region ownership
	// Key: ChunkCoord, Value: OceanID (-1 if not ocean)
	oceanChunks map[ChunkCoord]int

	// Tracks visited sites during river tracing within a chunk.
	// Reset before processing each chunk's hydrology.
	visitedSites map[int]bool

	// ID counters
	nextRiverID int
	nextLakeID  int
	nextOceanID int
}

// NewHydroManager creates a new hydrology manager.
func NewHydroManager(cfg HydroConfig, seed int64) *HydroManager {
	return &HydroManager{
		cfg:              cfg,
		rng:              rand.New(rand.NewSource(seed)),
		Rivers:           make(map[int]*River),
		Lakes:            make(map[int]*Lake),
		riverPropagation: make(map[ChunkCoord][]RiverInterPoint),
		oceanChunks:      make(map[ChunkCoord]int),
		visitedSites:     make(map[int]bool),
	}
}

// ResetVisitedSites clears the visited sites map for a new chunk's hydrology pass.
func (hm *HydroManager) ResetVisitedSites() {
	hm.visitedSites = make(map[int]bool)
}

// IsSourceVisited returns true if the given site has already been visited by a river trace.
func (hm *HydroManager) IsSourceVisited(site int) bool {
	return hm.visitedSites != nil && hm.visitedSites[site]
}

// -------------------------------------------------------------------
// River Tracing Algorithm
// -------------------------------------------------------------------

// TraceRiver traces a river path from a starting position using gradient descent.
// Returns the river segment and any exit point if the river leaves the chunk bounds.
func (hm *HydroManager) TraceRiver(
	mesh *DelaunayMesh,
	heights []float64,
	startSite int,	// source - A Site/vertex
	chunkBounds struct{ MinX, MinZ, MaxX, MaxZ float64 },
	initialWidth, initialDepth, distanceTraveled float64,
) (RiverSegment, *RiverInterPoint) {

	segment := RiverSegment{
		Vertices: make([]Vec2, 0, 64),
		Widths:   make([]float64, 0, 64),
		Depths:   make([]float64, 0, 64),
	}

	current := startSite
	visited := make(map[int]bool)
	distance := distanceTraveled
	width := initialWidth
	depth := initialDepth

	var flowDirAccum Vec2
	var exitPoint *RiverInterPoint

	for iterations := 0; iterations < 10000; iterations++ {
		if visited[current] {
			break // Cycle detected, stop
		}
		visited[current] = true
		if hm.visitedSites != nil {
			hm.visitedSites[current] = true
		}

		pos := mesh.Sites[current].Pos
		h := heights[current]

		// Record vertex
		segment.Vertices = append(segment.Vertices, pos)
		segment.Widths = append(segment.Widths, width)
		segment.Depths = append(segment.Depths, depth)

		// Check if we've exited chunk bounds
		if pos.X < chunkBounds.MinX || pos.X >= chunkBounds.MaxX ||
			pos.Y < chunkBounds.MinZ || pos.Y >= chunkBounds.MaxZ {

			// Determine exit edge
			exitEdge := 0
			if pos.X < chunkBounds.MinX {
				exitEdge = 0
			} else if pos.X >= chunkBounds.MaxX {
				exitEdge = 1
			} else if pos.Y < chunkBounds.MinZ {
				exitEdge = 2
			} else {
				exitEdge = 3
			}

			exitPoint = &RiverInterPoint{
				Position: pos,
				Width:    width,
				Depth:    depth,
				FlowDir:  flowDirAccum.Normalize(),
				Distance: distance,
				ExitEdge: exitEdge,
			}
			break
		}

		// Find neighbors with lower elevation
		neighbors := hm.getLowerNeighbors(mesh, heights, current)
		if len(neighbors) == 0 {
			// Local minimum - potential lake site
			break
		}

		// Apply flow bias for tie-breaking
		bias := hm.cfg.FlowBias(pos.Y) // pos.Y is Z in world space
		next := hm.selectNextVertex(mesh, heights, current, neighbors, h, bias)

		// Calculate step
		nextPos := mesh.Sites[next].Pos
		stepVec := nextPos.Sub(pos)
		stepDist := stepVec.Len()

		// Accumulate flow direction
		if stepDist > 1e-9 {
			flowDirAccum = flowDirAccum.Add(stepVec.Mul(1.0 / stepDist))
		}

		// Update river properties
		distance += stepDist
		width = hm.cfg.RiverWidthBase + distance*hm.cfg.RiverWidthGrowth
		depth = hm.cfg.RiverDepthBase + distance*hm.cfg.RiverDepthGrowth

		current = next
	}

	// Compute average flow direction
	if len(segment.Vertices) > 1 {
		segment.FlowDir = flowDirAccum.Normalize()
	}

	return segment, exitPoint
}

// getLowerNeighbors returns site indices of neighbors with lower elevation.
func (hm *HydroManager) getLowerNeighbors(mesh *DelaunayMesh, heights []float64, site int) []int {
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
		if heights[dest] < currentHeight-hm.cfg.RiverMinSlope {
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
	mesh *DelaunayMesh,
	heights []float64,
	current int,
	candidates []int,
	currentHeight float64,
	bias Vec2,
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
		score *= 1.0 + (hm.rng.Float64()-0.5)*0.1

		if score > bestScore {
			bestScore = score
			bestCandidate = c
		}
	}

	return bestCandidate
}

func (hm *HydroManager) TraceRiverV2(
	mesh *DelaunayMesh,
	heights []float64,
	startSite int,	// source - A Site/vertex
	chunkBounds struct{ MinX, MinZ, MaxX, MaxZ float64 },
	initialWidth, initialDepth, distanceTraveled float64,
) (RiverSegment, *RiverInterPoint) {

	segment := RiverSegment{
		Vertices: make([]Vec2, 0, 64),
		Widths:   make([]float64, 0, 64),
		Depths:   make([]float64, 0, 64),
	}

	current := startSite
	visited := make(map[int]bool)
	distance := distanceTraveled
	width := initialWidth
	depth := initialDepth

	var flowDirAccum Vec2
	var exitPoint *RiverInterPoint

	// Note: Iterations is an arbitrary limit
	for iterations := 0; iterations < 10000; iterations++ {
		if visited[current] {
			break // Cycle detected, stop
		}
		visited[current] = true
		if hm.visitedSites != nil {
			hm.visitedSites[current] = true
		}

		pos := mesh.Sites[current].Pos
		// h := heights[current]

		// Note: First iteration records vertex unless already visited somehow
		// Record vertex
		segment.Vertices = append(segment.Vertices, pos)
		segment.Widths = append(segment.Widths, width)
		segment.Depths = append(segment.Depths, depth)

		// Check if we've exited chunk bounds
		if pos.X < chunkBounds.MinX || pos.X >= chunkBounds.MaxX ||
			pos.Y < chunkBounds.MinZ || pos.Y >= chunkBounds.MaxZ {

			// Determine exit edge
			exitEdge := 0
			if pos.X < chunkBounds.MinX {
				exitEdge = 0
			} else if pos.X >= chunkBounds.MaxX {
				exitEdge = 1
			} else if pos.Y < chunkBounds.MinZ {
				exitEdge = 2
			} else {
				exitEdge = 3
			}

			exitPoint = &RiverInterPoint{
				Position: pos,
				Width:    width,
				Depth:    depth,
				FlowDir:  flowDirAccum.Normalize(),
				Distance: distance,
				ExitEdge: exitEdge,
			}
			break
		}

		bias := hm.cfg.FlowBias(pos.Y) // pos.Y is Z in world space

		// Find neighbors with lower elevation
		neighbors := hm.getFlowBiasedNeighbors(mesh, current, bias)
		if len(neighbors) == 0 {
			// Local minimum - potential lake site
			break
		}

		// Select the next segment vertex (Site Index)
		next := hm.selectNextVertexV2(mesh, current, neighbors, bias)

		// Calculate step
		nextPos := mesh.Sites[next].Pos
		stepVec := nextPos.Sub(pos)
		stepDist := stepVec.Len()

		// Accumulate flow direction
		if stepDist > 1e-9 {
			flowDirAccum = flowDirAccum.Add(stepVec.Mul(1.0 / stepDist))
		}

		// Update river properties
		distance += stepDist
		width = hm.cfg.RiverWidthBase + distance*hm.cfg.RiverWidthGrowth
		depth = hm.cfg.RiverDepthBase + distance*hm.cfg.RiverDepthGrowth

		current = next
	}

	// Compute average flow direction
	if len(segment.Vertices) > 1 {
		segment.FlowDir = flowDirAccum.Normalize()
	}

	return segment, exitPoint

}

func (hm *HydroManager) getFlowBiasedNeighbors(mesh *DelaunayMesh, site int, bias Vec2) []int {

	neighbors := make([]int, 0, 8)

	// Use half-edge structure to find neighbors
	startEdge := mesh.SiteEdge[site]
	if startEdge == -1 {
		return neighbors
	}

	edge := startEdge
	visited := make(map[int]bool)

	currentPos := mesh.Sites[site].Pos

	for {
		if visited[edge] {
			break
		}
		visited[edge] = true

		// Get destination of this edge (Site index)
		dest := mesh.HalfEdges[edge].EdgeDest
		
		// Calculate the destination based on similarity to the flow bias
		destPos := mesh.Sites[dest].Pos
		destVec := destPos.Sub(currentPos)
		destDist := destVec.Len()
		if destDist < 1e-9 {
			continue
		}
		destDir := destVec.Mul(1.0 / destDist)
		biasAlignment := destDir.Dot(bias) // -1 to 1

		// Aligned (enough) with our flow bias
		if biasAlignment > 0.1 {
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

// Prefer the vertex most aligned with the flow bias
// NOTE: This will ignore slope, implying the vertex shader will reduce height to the selected vertex-Y, creating a valley.
func (hm *HydroManager) selectNextVertexV2(
	mesh *DelaunayMesh,
	current int,
	candidates []int,
	bias Vec2,
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

		// If we weight by slope we can cut mountains in half... make this a config. TODO

		score := biasAlignment*0.5

		// Add small random noise for meandering (±5%)
		score *= 1.0 + (hm.rng.Float64()-0.5)*0.1

		if score > bestScore {
			bestScore = score
			bestCandidate = c
		}
	}

	return bestCandidate

}

// -------------------------------------------------------------------
// Lake Detection
// -------------------------------------------------------------------

// FindLocalMinima identifies vertices that are lower than all their neighbors.
// These are potential lake sites.
func (hm *HydroManager) FindLocalMinima(mesh *DelaunayMesh, heights []float64) []int {
	minima := make([]int, 0)

	for site := range mesh.Sites {
		if hm.isLocalMinimum(mesh, heights, site) {
			minima = append(minima, site)
		}
	}

	return minima
}

// isLocalMinimum returns true if the site is lower than all its neighbors.
func (hm *HydroManager) isLocalMinimum(mesh *DelaunayMesh, heights []float64, site int) bool {
	currentHeight := heights[site]

	startEdge := mesh.SiteEdge[site]
	if startEdge == -1 {
		return false
	}

	edge := startEdge
	visited := make(map[int]bool)

	for {
		if visited[edge] {
			break
		}
		visited[edge] = true

		dest := mesh.HalfEdges[edge].EdgeDest
		if heights[dest] <= currentHeight {
			return false // Not a minimum
		}

		twin := mesh.HalfEdges[edge].Twin
		if twin == -1 {
			break
		}
		edge = mesh.HalfEdges[twin].Next
		if edge == startEdge {
			break
		}
	}

	return true
}

// FloodFillLake fills a depression from a local minimum to find lake extent.
// Returns the lake if valid, nil otherwise.
func (hm *HydroManager) FloodFillLake(
	mesh *DelaunayMesh,
	heights []float64,
	minimum int,
) *Lake {
	// Check probability - if 0, never create lakes
	if hm.rng.Float64() >= hm.cfg.LakeProbability {
		return nil
	}

	minHeight := heights[minimum]

	// Find all connected vertices below the spillway level
	filled := make(map[int]bool)
	frontier := []int{minimum}
	filled[minimum] = true

	var spillway *int
	spillwayHeight := math.Inf(1)

	for len(frontier) > 0 {
		current := frontier[0]
		frontier = frontier[1:]

		// Check neighbors
		startEdge := mesh.SiteEdge[current]
		if startEdge == -1 {
			continue
		}

		edge := startEdge
		visited := make(map[int]bool)

		for {
			if visited[edge] {
				break
			}
			visited[edge] = true

			dest := mesh.HalfEdges[edge].EdgeDest
			destHeight := heights[dest]

			if !filled[dest] {
				if destHeight < spillwayHeight {
					// This could be part of the lake or the spillway
					if destHeight <= minHeight+hm.cfg.LakeMinDepth*2 {
						// Still in the depression
						filled[dest] = true
						frontier = append(frontier, dest)
					} else if destHeight < spillwayHeight {
						// Potential spillway
						spillwayHeight = destHeight
						spillway = &dest
					}
				}
			}

			twin := mesh.HalfEdges[edge].Twin
			if twin == -1 {
				break
			}
			edge = mesh.HalfEdges[twin].Next
			if edge == startEdge {
				break
			}
		}
	}

	// Check if lake meets minimum requirements
	if len(filled) < hm.cfg.LakeMinArea {
		return nil
	}

	// Calculate lake properties
	siteIndices := make([]int, 0, len(filled))
	var centerX, centerZ float64
	maxDepth := 0.0

	for site := range filled {
		siteIndices = append(siteIndices, site)
		pos := mesh.Sites[site].Pos
		centerX += pos.X
		centerZ += pos.Y
		depth := spillwayHeight - heights[site]
		if depth > maxDepth {
			maxDepth = depth
		}
	}

	centerX /= float64(len(siteIndices))
	centerZ /= float64(len(siteIndices))

	lake := &Lake{
		ID:          hm.nextLakeID,
		CenterPos:   Vec2{centerX, centerZ},
		SiteIndices: siteIndices,
		WaterLevel:  spillwayHeight,
		MinDepth:    maxDepth,
	}

	if spillway != nil {
		sp := mesh.Sites[*spillway].Pos
		lake.Spillway = &sp

		// Calculate spillway direction (away from lake center)
		lake.SpillwayDir = sp.Sub(lake.CenterPos).Normalize()
	}

	hm.nextLakeID++
	return lake
}

// -------------------------------------------------------------------
// Ocean Region System
// -------------------------------------------------------------------

// OceanSeedPoint represents a potential ocean origin.
type OceanSeedPoint struct {
	ChunkX int
	ChunkZ int
	Seed   int64
}

// GenerateOceanRegion creates an ocean region starting from a seed chunk.
// The ocean expands from the seed up to the configured maximum size.
func (hm *HydroManager) GenerateOceanRegion(seedChunk ChunkCoord, worldSeed int64) *OceanRegion {
	// Create deterministic RNG for this ocean
	oceanSeed := combineSeeds(worldSeed, int64(seedChunk.X), int64(seedChunk.Z))
	rng := rand.New(rand.NewSource(oceanSeed))

	// Check probability
	if rng.Float64() > hm.cfg.OceanProbability {
		return nil
	}

	// Determine ocean extent (random within max bounds)
	extentX := rng.Intn(hm.cfg.OceanMaxChunksX) + 1
	extentZ := rng.Intn(hm.cfg.OceanMaxChunksZ) + 1

	// Center the ocean on the seed chunk (with some randomness)
	offsetX := rng.Intn(extentX)
	offsetZ := rng.Intn(extentZ)

	ocean := &OceanRegion{
		ID:        hm.nextOceanID,
		MinChunkX: seedChunk.X - offsetX,
		MinChunkZ: seedChunk.Z - offsetZ,
		MaxChunkX: seedChunk.X - offsetX + extentX - 1,
		MaxChunkZ: seedChunk.Z - offsetZ + extentZ - 1,
		SeaLevel:  hm.cfg.SeaLevel,
	}

	// Register chunks as ocean
	for x := ocean.MinChunkX; x <= ocean.MaxChunkX; x++ {
		for z := ocean.MinChunkZ; z <= ocean.MaxChunkZ; z++ {
			hm.oceanChunks[ChunkCoord{x, z}] = ocean.ID
		}
	}

	hm.Oceans = append(hm.Oceans, *ocean)
	hm.nextOceanID++

	return ocean
}

// IsOceanChunk returns true if the chunk is part of an ocean region.
func (hm *HydroManager) IsOceanChunk(coord ChunkCoord) bool {
	_, ok := hm.oceanChunks[coord]
	return ok
}

// GetOceanForChunk returns the ocean region for a chunk, or nil if not ocean.
func (hm *HydroManager) GetOceanForChunk(coord ChunkCoord) *OceanRegion {
	oceanID, ok := hm.oceanChunks[coord]
	if !ok {
		return nil
	}
	for i := range hm.Oceans {
		if hm.Oceans[i].ID == oceanID {
			return &hm.Oceans[i]
		}
	}
	return nil
}

// ComputeOceanChunkData calculates ocean data for a specific chunk.
func (hm *HydroManager) ComputeOceanChunkData(
	coord ChunkCoord,
	mesh *DelaunayMesh,
	heights []float64,
) OceanChunkData {
	ocean := hm.GetOceanForChunk(coord)
	if ocean == nil {
		return OceanChunkData{IsOcean: false, OceanID: -1}
	}

	data := OceanChunkData{
		IsOcean:       true,
		OceanID:       ocean.ID,
		SubmergedTris: make([]int, 0),
		CoastlinePts:  make([]Vec2, 0),
	}

	// Find submerged triangles and coastline
	for triID, tri := range mesh.Tris {
		hA := heights[tri.A]
		hB := heights[tri.B]
		hC := heights[tri.C]

		belowA := hA <= ocean.SeaLevel
		belowB := hB <= ocean.SeaLevel
		belowC := hC <= ocean.SeaLevel

		if belowA && belowB && belowC {
			// Fully submerged
			data.SubmergedTris = append(data.SubmergedTris, triID)
		} else if belowA || belowB || belowC {
			// Partial - coastline passes through this triangle
			// Record the crossing points for coastline rendering
			data.SubmergedTris = append(data.SubmergedTris, triID)

			// Find edge crossings
			pA := mesh.Sites[tri.A].Pos
			pB := mesh.Sites[tri.B].Pos
			pC := mesh.Sites[tri.C].Pos

			if belowA != belowB {
				t := (ocean.SeaLevel - hA) / (hB - hA)
				crossing := Vec2{pA.X + t*(pB.X-pA.X), pA.Y + t*(pB.Y-pA.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
			if belowB != belowC {
				t := (ocean.SeaLevel - hB) / (hC - hB)
				crossing := Vec2{pB.X + t*(pC.X-pB.X), pB.Y + t*(pC.Y-pB.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
			if belowC != belowA {
				t := (ocean.SeaLevel - hC) / (hA - hC)
				crossing := Vec2{pC.X + t*(pA.X-pC.X), pC.Y + t*(pA.Y-pC.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
		}
	}

	return data
}

// -------------------------------------------------------------------
// Cross-Chunk River Propagation
// -------------------------------------------------------------------

// RegisterRiverExit records a river exiting a chunk for propagation.
func (hm *HydroManager) RegisterRiverExit(fromChunk ChunkCoord, exit RiverInterPoint) {
	// Determine destination chunk based on exit edge
	var destChunk ChunkCoord
	switch exit.ExitEdge {
	case 0: // minX
		destChunk = ChunkCoord{fromChunk.X - 1, fromChunk.Z}
	case 1: // maxX
		destChunk = ChunkCoord{fromChunk.X + 1, fromChunk.Z}
	case 2: // minZ
		destChunk = ChunkCoord{fromChunk.X, fromChunk.Z - 1}
	case 3: // maxZ
		destChunk = ChunkCoord{fromChunk.X, fromChunk.Z + 1}
	}

	hm.riverPropagation[destChunk] = append(hm.riverPropagation[destChunk], exit)
}

// GetRiverEntries returns rivers that should enter this chunk from neighbors.
func (hm *HydroManager) GetRiverEntries(chunk ChunkCoord) []RiverInterPoint {
	return hm.riverPropagation[chunk]
}

// ClearRiverEntries removes processed river entries for a chunk.
func (hm *HydroManager) ClearRiverEntries(chunk ChunkCoord) {
	delete(hm.riverPropagation, chunk)
}

// -------------------------------------------------------------------
// River Source Detection
// -------------------------------------------------------------------

// FindRiverSources identifies potential river sources in a chunk.
// Sources are considered at each Site (a triangle vertex)
// TODO: Consider the circumcenter that we're already plotting for the Voronoi diagram
// Sources are high-elevation vertices that pass the probability check.
func (hm *HydroManager) FindRiverSources(
	mesh *DelaunayMesh,
	heights []float64,
	coreSiteIndices []int,
) []int {
	sources := make([]int, 0)

	candidateCount := 0
	for _, site := range coreSiteIndices {
		h := heights[site]

		// Check elevation threshold
		if h < hm.cfg.RiverSourceMinElev {
			continue
		}

		// TODO - This prevents sources in basins. See docs
		// Check if this is a local maximum or near-maximum
		// (sources shouldn't be in valleys)
		if !hm.isNearLocalMaximum(mesh, heights, site) {
			continue
		}

		candidateCount++
		randVal := hm.rng.Float64()
		// Apply probability
		if randVal > hm.cfg.RiverSourceProbability {
			continue
		}

		sources = append(sources, site)
	}

	return sources
}


/*
	TODO: This is just one of many ways to pick a candidate. 
	Primary traversal mechanism is our half-edge mesh. (Acceleration struct)
*/
// isNearLocalMaximum returns true if the site is higher than most neighbors.
func (hm *HydroManager) isNearLocalMaximum(mesh *DelaunayMesh, heights []float64, site int) bool {
	currentHeight := heights[site]
	higherCount := 0
	totalNeighbors := 0

	startEdge := mesh.SiteEdge[site]
	if startEdge == -1 {
		return false
	}

	edge := startEdge
	visited := make(map[int]bool)

	for {

		// On first repeat visit to an edge. Break
		if visited[edge] {
			break
		}
		visited[edge] = true

		dest := mesh.HalfEdges[edge].EdgeDest
		totalNeighbors++
		if heights[dest] > currentHeight {
			higherCount++
		}

		// -1 is an init index. Might imply a min/max boundary... moot
		twin := mesh.HalfEdges[edge].Twin
		if twin == -1 {
			fmt.Printf("[isNearLocalMaximum] Twin is -1\nSite: %d\nEdge:\t %d\n", site, edge)
			break
		}
		edge = mesh.HalfEdges[twin].Next
		if edge == startEdge {
			break
		}
	}

	// Source should have at most 1/3 of neighbors higher
	return totalNeighbors > 0 && higherCount <= totalNeighbors/3
}

// -------------------------------------------------------------------
// Utility Functions
// -------------------------------------------------------------------

// combineSeeds creates a deterministic seed from multiple values.
func combineSeeds(a, b, c int64) int64 {
	// Simple hash combination
	h := a
	h = h*31 + b
	h = h*31 + c
	return h
}

// Normalize returns a unit vector (or zero vector if length is zero).
func (v Vec2) Normalize() Vec2 {
	l := v.Len()
	if l < 1e-12 {
		return Vec2{}
	}
	return v.Mul(1.0 / l)
}
