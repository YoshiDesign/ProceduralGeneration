package duals2

import (
	"fmt"
	"hash/fnv"
	"math/rand"
	"sync"

	"github.com/hajimehoshi/ebiten/v2"
)

// NoiseParams holds configurable parameters for fractal noise terrain generation.
type NoiseParams struct {
	Octaves     int
	Frequency   float64
	Amplitude   float64
	Persistence float64
	Lacunarity  float64
}

// DefaultNoiseParams returns sensible defaults for terrain noise generation.
func DefaultNoiseParams() NoiseParams {
	return NoiseParams{
		Octaves:     6,
		Frequency:   0.1,
		Amplitude:   3.0,
		Persistence: 0.2,
		Lacunarity:  2.0,
	}
}

// ChunkCoord identifies a chunk by integer grid coordinates.
type ChunkCoord struct {
	X, Z int
}

// ChunkManager manages chunk generation and caching for optimized neighbor lookups.
type ChunkManager struct {
	mu          sync.RWMutex
	cache       map[ChunkCoord]*TerrainChunk
	pointsCache map[ChunkCoord][]Vec2 // Raw blue noise points (before full chunk is built)
	cfg         ChunkConfig
	noiseParams NoiseParams
	heightFunc  func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64

	// Hydrology system
	hydro *HydroManager
}


// ChunkConfig holds parameters for terrain chunk generation.
type ChunkConfig struct {
	ChunkSize   float64 // World units per chunk side (e.g., 256.0)
	MinPointDist float64 // Minimum distance between blue noise points
	HaloWidth   float64 // Boundary overlap region width (typically = MinPointDist)
	WorldSeed   int64   // Global world seed
	ChunksX     int     // Number of chunks along the X axis
	ChunksZ     int     // Number of chunks along the Z axis
}

// DefaultChunkConfig returns sensible defaults for a terrain chunk.
func DefaultChunkConfig() ChunkConfig {
	return ChunkConfig{
		ChunkSize:    128.0,
		MinPointDist: 8.0,
		HaloWidth:    8.0,
		WorldSeed:    42,
		ChunksX: 	  4,
		ChunksZ: 	  3,
	}
}

// TerrainChunk represents a generated terrain chunk with mesh data.
type TerrainChunk struct {
	Coord ChunkCoord
	Cfg   ChunkConfig

	// Core bounds (what this chunk "owns")
	MinX, MinZ float64
	MaxX, MaxZ float64

	// The Delaunay mesh (includes halo points for boundary continuity)
	Mesh *DelaunayMesh

	// Height values per site (parallel to Mesh.Sites)
	Heights []float64

	// Face normals per triangle (parallel to Mesh.Tris)
	FaceNormals []Vec3

	// Spatial index for fast point location
	Spatial *SpatialGrid

	// Which sites are in the core region (not halo)
	CoreSiteIndices []int

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
}

// chunkSeed computes a deterministic seed for a chunk based on world seed and coordinates.
func chunkSeed(worldSeed int64, coord ChunkCoord) int64 {
	h := fnv.New64a()
	// Write world seed
	buf := make([]byte, 8)
	buf[0] = byte(worldSeed)
	buf[1] = byte(worldSeed >> 8)
	buf[2] = byte(worldSeed >> 16)
	buf[3] = byte(worldSeed >> 24)
	buf[4] = byte(worldSeed >> 32)
	buf[5] = byte(worldSeed >> 40)
	buf[6] = byte(worldSeed >> 48)
	buf[7] = byte(worldSeed >> 56)
	h.Write(buf)

	// Write chunk X
	buf[0] = byte(coord.X)
	buf[1] = byte(coord.X >> 8)
	buf[2] = byte(coord.X >> 16)
	buf[3] = byte(coord.X >> 24)
	buf[4] = 0
	buf[5] = 0
	buf[6] = 0
	buf[7] = 0
	h.Write(buf)

	// Write chunk Z
	buf[0] = byte(coord.Z)
	buf[1] = byte(coord.Z >> 8)
	buf[2] = byte(coord.Z >> 16)
	buf[3] = byte(coord.Z >> 24)
	h.Write(buf)

	return int64(h.Sum64())
}

// boundarySeed computes a deterministic seed for the shared boundary between two chunks.
// It uses the minimum of the two chunk coords to ensure both chunks get the same seed.
func boundarySeed(worldSeed int64, coord1, coord2 ChunkCoord) int64 {
	// Use lexicographically smaller coord first
	var first, second ChunkCoord
	if coord1.X < coord2.X || (coord1.X == coord2.X && coord1.Z < coord2.Z) {
		first, second = coord1, coord2
	} else {
		first, second = coord2, coord1
	}

	h := fnv.New64a()
	buf := make([]byte, 8)

	// World seed
	buf[0] = byte(worldSeed)
	buf[1] = byte(worldSeed >> 8)
	buf[2] = byte(worldSeed >> 16)
	buf[3] = byte(worldSeed >> 24)
	buf[4] = byte(worldSeed >> 32)
	buf[5] = byte(worldSeed >> 40)
	buf[6] = byte(worldSeed >> 48)
	buf[7] = byte(worldSeed >> 56)
	h.Write(buf)

	// First chunk
	buf[0] = byte(first.X)
	buf[1] = byte(first.X >> 8)
	buf[2] = byte(first.Z)
	buf[3] = byte(first.Z >> 8)
	// Second chunk
	buf[4] = byte(second.X)
	buf[5] = byte(second.X >> 8)
	buf[6] = byte(second.Z)
	buf[7] = byte(second.Z >> 8)
	h.Write(buf)

	return int64(h.Sum64())
}

// NewChunkManager creates a new chunk manager with the given configuration.
func NewChunkManager(cfg ChunkConfig, heightFunc func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64) *ChunkManager {
	return &ChunkManager{
		cache:       make(map[ChunkCoord]*TerrainChunk),
		pointsCache: make(map[ChunkCoord][]Vec2),
		cfg:         cfg,
		noiseParams: DefaultNoiseParams(),
		heightFunc:  heightFunc,
		hydro:       NewHydroManager(DefaultHydroConfig(), cfg.WorldSeed),
	}
}

// NewChunkManagerWithHydro creates a chunk manager with custom hydrology configuration.
func NewChunkManagerWithHydro(cfg ChunkConfig, hydroCfg HydroConfig, heightFunc func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64) *ChunkManager {
	return &ChunkManager{
		cache:       make(map[ChunkCoord]*TerrainChunk),
		pointsCache: make(map[ChunkCoord][]Vec2),
		cfg:         cfg,
		noiseParams: DefaultNoiseParams(),
		heightFunc:  heightFunc,
		hydro:       NewHydroManager(hydroCfg, cfg.WorldSeed),
	}
}

// NewChunkManagerFull creates a chunk manager with all custom configurations.
func NewChunkManagerFull(cfg ChunkConfig, noiseParams NoiseParams, hydroCfg HydroConfig, heightFunc func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64) *ChunkManager {
	return &ChunkManager{
		cache:       make(map[ChunkCoord]*TerrainChunk),
		pointsCache: make(map[ChunkCoord][]Vec2),
		cfg:         cfg,
		noiseParams: noiseParams,
		heightFunc:  heightFunc,
		hydro:       NewHydroManager(hydroCfg, cfg.WorldSeed),
	}
}

// HydroManager returns the hydrology manager for advanced hydrology operations.
func (cm *ChunkManager) HydroManager() *HydroManager {
	return cm.hydro
}

// NoiseParams returns the current noise parameters.
func (cm *ChunkManager) NoiseParams() NoiseParams {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return cm.noiseParams
}

// SetNoiseParams updates the noise parameters. Call ClearCaches() and regenerate chunks to apply.
func (cm *ChunkManager) SetNoiseParams(np NoiseParams) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.noiseParams = np
}

// HydroConfig returns the current hydrology configuration.
func (cm *ChunkManager) HydroConfig() HydroConfig {
	return cm.hydro.cfg
}

// SetHydroConfig updates the hydrology configuration. Call ClearCaches() and regenerate chunks to apply.
func (cm *ChunkManager) SetHydroConfig(cfg HydroConfig) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.hydro = NewHydroManager(cfg, cm.cfg.WorldSeed)
}

// ClearCaches removes all cached chunks and points, allowing regeneration with new parameters.
func (cm *ChunkManager) ClearCaches() {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.cache = make(map[ChunkCoord]*TerrainChunk)
	cm.pointsCache = make(map[ChunkCoord][]Vec2)
}

// GetOrGenerate returns a cached chunk or generates and caches a new one.
// This is the primary API for chunk access with caching optimization.
func (cm *ChunkManager) GetOrGenerate(coord ChunkCoord) (*TerrainChunk, error) {
	// Fast path: check if already cached
	cm.mu.RLock()
	if chunk, ok := cm.cache[coord]; ok {
		cm.mu.RUnlock()
		return chunk, nil
	}
	cm.mu.RUnlock()

	// Slow path: generate with cache-aware point collection
	cm.mu.Lock()
	defer cm.mu.Unlock()

	// Double-check after acquiring write lock
	if chunk, ok := cm.cache[coord]; ok {
		return chunk, nil
	}

	chunk, err := cm.generateChunkInternal(coord)
	if err != nil {
		return nil, err
	}

	cm.cache[coord] = chunk
	return chunk, nil
}

// Get returns a cached chunk without generating. Returns nil if not cached.
func (cm *ChunkManager) Get(coord ChunkCoord) *TerrainChunk {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return cm.cache[coord]
}

// Evict removes a chunk from the cache, freeing memory.
func (cm *ChunkManager) Evict(coord ChunkCoord) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	delete(cm.cache, coord)
}

// EvictOutsideRadius removes all chunks outside the given radius from center.
func (cm *ChunkManager) EvictOutsideRadius(center ChunkCoord, radius int) int {
	cm.mu.Lock()
	defer cm.mu.Unlock()

	evicted := 0
	for coord := range cm.cache {
		dx := coord.X - center.X
		dz := coord.Z - center.Z
		if dx*dx+dz*dz > radius*radius {
			delete(cm.cache, coord)
			evicted++
		}
	}
	return evicted
}

// CacheSize returns the number of cached chunks.
func (cm *ChunkManager) CacheSize() int {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return len(cm.cache)
}

// CachedCoords returns all cached chunk coordinates.
func (cm *ChunkManager) CachedCoords() []ChunkCoord {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	coords := make([]ChunkCoord, 0, len(cm.cache))
	for coord := range cm.cache {
		coords = append(coords, coord)
	}
	return coords
}

// getOrGeneratePoints returns blue noise points for a chunk, using caches when available.
// Priority: 1) Full chunk cache (Mesh.Sites), 2) Points cache, 3) Generate new.
// Caller must hold at least a read lock on cm.mu.
// If generation is needed, caller should hold a write lock.
func (cm *ChunkManager) getOrGeneratePoints(coord ChunkCoord) []Vec2 {

	// Check if points are cached.
	if pts, ok := cm.pointsCache[coord]; ok {
		fmt.Println("Grabbing point cache...")
		return pts
	}

	// Optional sanity check. Check the chunk directly for the points
	// This will never happen under normal circumstances. Caching points is an invariant to generating blue noise.
	if chunk, ok := cm.cache[coord]; ok {
		points := make([]Vec2, len(chunk.Mesh.Sites))
		for i, site := range chunk.Mesh.Sites {
			points[i] = site.Pos
		}
		fmt.Println("Uhhhh...")
		return points
	}

	// 3. Generate blue noise and cache it
	minX := float64(coord.X) * cm.cfg.ChunkSize
	minZ := float64(coord.Z) * cm.cfg.ChunkSize
	maxX := float64(coord.X+1) * cm.cfg.ChunkSize
	maxZ := float64(coord.Z+1) * cm.cfg.ChunkSize

	blueCfg := DefaultBlueNoiseConfig(cm.cfg.MinPointDist)
	seed := chunkSeed(cm.cfg.WorldSeed, coord)
	pts := GenerateBlueNoiseSeeded(seed, minX, minZ, maxX, maxZ, blueCfg)

	// Cache the generated points
	cm.pointsCache[coord] = pts
	return pts
}

// GenerateChunk creates a terrain chunk with Delaunay mesh, heights, and spatial index.
// The heightFunc provides elevation for each site position.
// This version does not use caching - for cached generation, use ChunkManager.
func GenerateChunk(coord ChunkCoord, cfg ChunkConfig, heightFunc func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64) (*TerrainChunk, error) {
	// Create a temporary ChunkManager for standalone generation (no caching benefit)
	tempManager := NewChunkManager(cfg, heightFunc)
	return tempManager.generateChunkInternal(coord)
}

// generateChunkInternal creates a terrain chunk using the ChunkManager's caches.
// Caller must hold a write lock on cm.mu.
func (cm *ChunkManager) generateChunkInternal(coord ChunkCoord) (*TerrainChunk, error) {
	chunk := &TerrainChunk{
		Coord: coord,
		Cfg:   cm.cfg,
		MinX:  float64(coord.X) * cm.cfg.ChunkSize,
		MinZ:  float64(coord.Z) * cm.cfg.ChunkSize,
		MaxX:  float64(coord.X+1) * cm.cfg.ChunkSize,
		MaxZ:  float64(coord.Z+1) * cm.cfg.ChunkSize,
	}

	// Generate points: core region + halo regions for each neighbor
	allPoints := cm.generateChunkPoints(coord)

	// Track which points are in the core region
	coreIndices := make([]int, 0, len(allPoints))
	for i, p := range allPoints {
		if p.X >= chunk.MinX && p.X < chunk.MaxX && p.Y >= chunk.MinZ && p.Y < chunk.MaxZ {
			coreIndices = append(coreIndices, i)
		}
	}
	chunk.CoreSiteIndices = coreIndices

	// Build sites with heights
	sites := make([]Site, len(allPoints))
	heights := make([]float64, len(allPoints))
	np := cm.noiseParams
	for i, p := range allPoints {
		h := 0.0
		if cm.heightFunc != nil {
			h = cm.heightFunc(p.X, p.Y, np.Octaves, np.Frequency, np.Amplitude, np.Persistence, np.Lacunarity)
		}
		sites[i] = Site{Pos: p, Height: h}
		heights[i] = h
	}
	chunk.Heights = heights

	// Build Delaunay triangulation
	tris := Triangulate(allPoints)
	mesh, err := BuildHalfEdgeMesh(sites, tris)
	if err != nil {
		return nil, err
	}
	chunk.Mesh = mesh

	// Compute face normals
	chunk.FaceNormals = mesh.AllFaceNormals(heights)

	// Build spatial index
	// Cell size roughly equal to minimum point distance for good performance
	chunk.Spatial = BuildSpatialIndex(mesh, heights, cm.cfg.MinPointDist,
		chunk.MinX-cm.cfg.HaloWidth, chunk.MinZ-cm.cfg.HaloWidth,
		chunk.MaxX+cm.cfg.HaloWidth, chunk.MaxZ+cm.cfg.HaloWidth)

	// Generate hydrology data
	chunk.Hydro = cm.generateChunkHydrology(chunk)

	return chunk, nil
}

// generateChunkHydrology computes rivers, lakes, and ocean data for a chunk.
func (cm *ChunkManager) generateChunkHydrology(chunk *TerrainChunk) *ChunkHydroData {
	hydro := &ChunkHydroData{
		Rivers:       make([]RiverSegment, 0),
		Lakes:        make([]Lake, 0),
		RiverEntries: make([]RiverExitPoint, 0),
		RiverExits:   make([]RiverExitPoint, 0),
	}

	bounds := struct{ MinX, MinZ, MaxX, MaxZ float64 }{
		chunk.MinX, chunk.MinZ, chunk.MaxX, chunk.MaxZ,
	}

	// 1. Process incoming rivers from neighbors
	entries := cm.hydro.GetRiverEntries(chunk.Coord)
	for _, entry := range entries {
		hydro.RiverEntries = append(hydro.RiverEntries, entry)

		// Find the nearest site to the entry point
		nearestSite := cm.findNearestSite(chunk.Mesh, entry.Position)
		if nearestSite < 0 {
			continue
		}

		// Continue tracing the river
		segment, exit := cm.hydro.TraceRiver(
			chunk.Mesh, chunk.Heights, nearestSite, bounds,
			entry.Width, entry.Depth, entry.Distance,
		)

		if len(segment.Vertices) > 0 {
			hydro.Rivers = append(hydro.Rivers, segment)
		}

		if exit != nil {
			exit.RiverID = entry.RiverID
			hydro.RiverExits = append(hydro.RiverExits, *exit)
			cm.hydro.RegisterRiverExit(chunk.Coord, *exit)
		}
	}
	cm.hydro.ClearRiverEntries(chunk.Coord)

	// 2. Find new river sources in this chunk
	sources := cm.hydro.FindRiverSources(chunk.Mesh, chunk.Heights, chunk.CoreSiteIndices)
	for _, source := range sources {
		segment, exit := cm.hydro.TraceRiver(
			chunk.Mesh, chunk.Heights, source, bounds,
			cm.hydro.cfg.RiverWidthBase, cm.hydro.cfg.RiverDepthBase, 0,
		)

		if len(segment.Vertices) > 0 {
			hydro.Rivers = append(hydro.Rivers, segment)
		}

		if exit != nil {
			exit.RiverID = cm.hydro.nextRiverID
			cm.hydro.nextRiverID++
			hydro.RiverExits = append(hydro.RiverExits, *exit)
			cm.hydro.RegisterRiverExit(chunk.Coord, *exit)
		}
	}

	// 3. Detect lakes at local minima
	minima := cm.hydro.FindLocalMinima(chunk.Mesh, chunk.Heights)
	for _, minimum := range minima {
		// Check if this minimum is in the core region
		isCore := false
		for _, coreIdx := range chunk.CoreSiteIndices {
			if coreIdx == minimum {
				isCore = true
				break
			}
		}
		if !isCore {
			continue
		}

		lake := cm.hydro.FloodFillLake(chunk.Mesh, chunk.Heights, minimum)
		if lake != nil {
			hydro.Lakes = append(hydro.Lakes, *lake)
		}
	}

	// 4. Compute ocean data if this chunk is in an ocean region
	hydro.Ocean = cm.hydro.ComputeOceanChunkData(chunk.Coord, chunk.Mesh, chunk.Heights)

	return hydro
}

// findNearestSite finds the mesh site closest to a given position.
func (cm *ChunkManager) findNearestSite(mesh *DelaunayMesh, pos Vec2) int {
	bestDist := float64(1e18)
	bestSite := -1

	for i, site := range mesh.Sites {
		d := site.Pos.Sub(pos).Len2()
		if d < bestDist {
			bestDist = d
			bestSite = i
		}
	}

	return bestSite
}

// generateChunkPoints generates blue noise points for a chunk including halo regions.
// Uses the point cache to avoid redundant blue noise generation for neighbors.
// Caller must hold a write lock on cm.mu.
func (cm *ChunkManager) generateChunkPoints(coord ChunkCoord) []Vec2 {

	
	minX := float64(coord.X) * cm.cfg.ChunkSize
	minZ := float64(coord.Z) * cm.cfg.ChunkSize
	maxX := float64(coord.X+1) * cm.cfg.ChunkSize
	maxZ := float64(coord.Z+1) * cm.cfg.ChunkSize

	fmt.Printf(
		"generateChunkPoints-------------\n ChunkSize\t%v\ncoord.X\t%v\ncoord.Z\t%v\nmin:\t(%v, %v)\nmax:(%v, %v)\n---------------\n", 
		cm.cfg.ChunkSize, coord.X, coord.Z, minX, minZ, maxX, maxZ)

	halo := cm.cfg.HaloWidth

	// We'll collect points from multiple regions
	allPoints := make([]Vec2, 0, 1024)
	seen := make(map[uint64]struct{}, 1024)

	// Hash a point to detect duplicates (within tolerance)
	hashPoint := func(p Vec2) uint64 {
		// Quantize to half the min distance for dedup
		scale := 2.0 / cm.cfg.MinPointDist
		qx := int64(p.X * scale)
		qz := int64(p.Y * scale)
		return uint64(qx)<<32 | uint64(qz)&0xFFFFFFFF
	}

	addPoint := func(p Vec2) bool {
		h := hashPoint(p)
		if _, exists := seen[h]; !exists {
			seen[h] = struct{}{}
			allPoints = append(allPoints, p)
			return true
		}
		return false
	}

	addPoints := func(pts []Vec2) {
		for _, p := range pts {
			addPoint(p)
		}
	}

	// 1. Get or generate core points for this chunk (uses cache if available)
	corePoints := cm.getOrGeneratePoints(coord)
	addPoints(corePoints)

	// 2. Get halo points from each neighboring chunk
	neighbors := []ChunkCoord{
		{coord.X - 1, coord.Z - 1}, {coord.X, coord.Z - 1}, {coord.X + 1, coord.Z - 1},
		{coord.X - 1, coord.Z}, {coord.X + 1, coord.Z},
		{coord.X - 1, coord.Z + 1}, {coord.X, coord.Z + 1}, {coord.X + 1, coord.Z + 1},
	}

	for _, neighbor := range neighbors {
		// The boundary region is where this chunk's halo overlaps the neighbor
		nMinX := float64(neighbor.X) * cm.cfg.ChunkSize
		nMinZ := float64(neighbor.Z) * cm.cfg.ChunkSize
		nMaxX := float64(neighbor.X+1) * cm.cfg.ChunkSize
		nMaxZ := float64(neighbor.Z+1) * cm.cfg.ChunkSize

		// Compute the overlap region between our extended bounds and neighbor's core
		overlapMinX := max(minX-halo, nMinX)
		overlapMinZ := max(minZ-halo, nMinZ)
		overlapMaxX := min(maxX+halo, nMaxX)
		overlapMaxZ := min(maxZ+halo, nMaxZ)

		if overlapMinX >= overlapMaxX || overlapMinZ >= overlapMaxZ {
			continue // No overlap
		}

		// Get neighbor's points (from chunk cache, points cache, or generate + cache)
		neighborPoints := cm.getOrGeneratePoints(neighbor)

		// Filter to points within our halo region but outside our core
		for _, p := range neighborPoints {
			inHalo := (p.X >= minX-halo && p.X < maxX+halo && p.Y >= minZ-halo && p.Y < maxZ+halo)
			inCore := (p.X >= minX && p.X < maxX && p.Y >= minZ && p.Y < maxZ)
			if inHalo && !inCore {
				addPoint(p)
			}
		}
	}

	return allPoints
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

	fmt.Println("MinX:", minX, "MaxX:", maxX)
	fmt.Println("MinZ:", minZ, "MaZZ:", maxZ)

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
	superA := Vec2{midX - 20*deltaMax, midZ - deltaMax}
	superB := Vec2{midX, midZ + 20*deltaMax}
	superC := Vec2{midX + 20*deltaMax, midZ - deltaMax}

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
			result = append(result, ensureCCW(allPts, t.a, t.b, t.c))
		}
	}

	// Log walking stats for performance validation
	if walkStats.TotalWalks > 0 {
		avgSteps := float64(walkStats.TotalSteps) / float64(walkStats.TotalWalks)
		avgBadTris := float64(walkStats.TotalBadTris) / float64(walkStats.TotalWalks)
		fmt.Printf("[Triangulate] n=%d: walks=%d, avgSteps=%.2f, maxSteps=%d, fallbacks=%d, avgBadTris=%.2f\n",
			len(points), walkStats.TotalWalks, avgSteps, walkStats.MaxSteps, walkStats.FallbackCount, avgBadTris)
	}

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
func ensureCCW(pts []Vec2, a, b, c int) Triangle {
	// Cross product of (b-a) × (c-a)
	cross := (pts[b].X-pts[a].X)*(pts[c].Y-pts[a].Y) - (pts[c].X-pts[a].X)*(pts[b].Y-pts[a].Y)
	if cross < 0 {
		return Triangle{A: a, C: b, B: c} // Swap b and c
	}
	return Triangle{A: a, B: b, C: c}
}

// SampleHeight returns the interpolated height at position (x, z).
// Returns (height, ok) where ok=false if the point is outside the mesh.
func (c *TerrainChunk) SampleHeight(x, z float64) (float64, bool) {
	if c.Spatial == nil {
		return 0, false
	}
	return c.Spatial.SampleHeight(x, z)
}

// Helper: seeded RNG for boundary regions
func boundaryRNG(worldSeed int64, coord1, coord2 ChunkCoord) *rand.Rand {
	seed := boundarySeed(worldSeed, coord1, coord2)
	return rand.New(rand.NewSource(seed))
}
