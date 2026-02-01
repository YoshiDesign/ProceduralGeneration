package duals2

import (
	"encoding/json"
	"fmt"
	"hash/fnv"
	"os"
	"procedural_generation/terrain_generation/duals2/core"
	eros "procedural_generation/terrain_generation/duals2/erosion"
	"procedural_generation/terrain_generation/duals2/hydro"
	"time"

	"sync"
)

// ChunkManager manages chunk generation and caching for optimized neighbor lookups.
type ChunkManager struct {
	mu          sync.RWMutex
	cache       map[core.ChunkCoord]*core.TerrainChunk
	pointsCache map[core.ChunkCoord][]core.Vec2 // Raw blue noise points (before full chunk is built)
	cfg         core.ChunkConfig
	noiseParams core.NoiseParams
	heightFunc  func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64
	chunkSeeds  map[core.ChunkCoord]int64

	// Hydrology system
	hydroMgr *hydro.HydroManager
	erosMgr *eros.ErosionManager
}

// DefaultNoiseParams returns sensible defaults for terrain noise generation.
func DefaultNoiseParams() core.NoiseParams {
	return core.NoiseParams{
		Octaves:     7,
		Frequency:   0.1,
		Amplitude:   3.0,
		Persistence: 0.2,
		Lacunarity:  2.0,
	}
}

// DefaultChunkConfig returns sensible defaults for a terrain chunk.
func DefaultChunkConfig() core.ChunkConfig {
	return core.ChunkConfig{
		ChunkSize:    128.0,
		MinPointDist: 8.0,
		HaloWidth:    8.0,
		WorldSeed:    42,
		ChunksX: 	  6,
		ChunksZ: 	  4,
	}
}

// NewChunkManager creates a new chunk manager with the given configuration.
func NewChunkManager(cfg core.ChunkConfig, heightFunc func(x, z float64, octaves int, frequency, amplitude, persistence, lacunarity float64) float64) *ChunkManager {
	return &ChunkManager{
		cache:       make(map[core.ChunkCoord]*core.TerrainChunk),
		pointsCache: make(map[core.ChunkCoord][]core.Vec2),
		cfg:         cfg,
		noiseParams: DefaultNoiseParams(),
		heightFunc:  heightFunc,
		hydroMgr:       hydro.NewHydroManager(hydro.DefaultHydroConfig(), cfg.WorldSeed),
	}
}

/* DEMO CODE - BEGIN */

// HydroManager returns the hydrology manager for advanced hydrology operations.
func (cm *ChunkManager) HydroManager() *hydro.HydroManager {
	return cm.hydroMgr
}

// NoiseParams returns the current noise parameters.
func (cm *ChunkManager) NoiseParams() core.NoiseParams {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return cm.noiseParams
}

// SetNoiseParams updates the noise parameters. Call ClearCaches() and regenerate chunks to apply.
func (cm *ChunkManager) SetNoiseParams(np core.NoiseParams) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.noiseParams = np
}

// hydro.HydroConfig returns the current hydrology configuration.
func (cm *ChunkManager) HydroConfig() hydro.HydroConfig {
	return cm.hydroMgr.Cfg
}

// Sethydro.HydroConfig updates the hydrology configuration. Call ClearCaches() and regenerate chunks to apply.
func (cm *ChunkManager) SetHydroConfig(cfg hydro.HydroConfig) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.hydroMgr = hydro.NewHydroManager(cfg, cm.cfg.WorldSeed)
}

// ClearCaches removes all cached chunks and points, allowing regeneration with new parameters.
func (cm *ChunkManager) ClearCaches() {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	cm.cache = make(map[core.ChunkCoord]*core.TerrainChunk)
	cm.pointsCache = make(map[core.ChunkCoord][]core.Vec2)
}

/* DEMO CODE - END */

// GetOrGenerate returns a cached chunk or generates and caches a new one.
// This is the primary API for chunk access with caching optimization.
func (cm *ChunkManager) GetOrGenerate(coord core.ChunkCoord) (*core.TerrainChunk, error) {
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

// // Get returns a cached chunk without generating. Returns nil if not cached.
// func (cm *ChunkManager) Get(coord core.ChunkCoord) *TerrainChunk {
// 	cm.mu.RLock()
// 	defer cm.mu.RUnlock()
// 	return cm.cache[coord]
// }

// // Evict removes a chunk from the cache, freeing memory.
// func (cm *ChunkManager) Evict(coord core.ChunkCoord) {
// 	cm.mu.Lock()
// 	defer cm.mu.Unlock()
// 	delete(cm.cache, coord)
// }

// // EvictOutsideRadius removes all chunks outside the given radius from center.
// func (cm *ChunkManager) EvictOutsideRadius(center core.ChunkCoord, radius int) int {
// 	cm.mu.Lock()
// 	defer cm.mu.Unlock()

// 	evicted := 0
// 	for coord := range cm.cache {
// 		dx := coord.X - center.X
// 		dz := coord.Z - center.Z
// 		if dx*dx+dz*dz > radius*radius {
// 			delete(cm.cache, coord)
// 			evicted++
// 		}
// 	}
// 	return evicted
// }

// generateChunkInternal creates a terrain chunk using the ChunkManager's caches.
// Caller must hold a write lock on cm.mu.
func (cm *ChunkManager) generateChunkInternal(coord core.ChunkCoord) (*core.TerrainChunk, error) {
	chunk := &core.TerrainChunk{
		Coord: coord,
		Cfg:   cm.cfg,
		MinX:  float64(coord.X) * cm.cfg.ChunkSize, // world-space coordinate / world units per chunk side
		MinZ:  float64(coord.Z) * cm.cfg.ChunkSize, 
		MaxX:  float64(coord.X+1) * cm.cfg.ChunkSize,
		MaxZ:  float64(coord.Z+1) * cm.cfg.ChunkSize,
	}

	// Generate points: core region + halo region of 8 neighboring chunks for this chunk
	allPoints := cm.generateChunkPoints(coord)

	// Track which points are in the core region
	coreIndices := make([]core.SiteIndex, 0, len(allPoints))
	// TODO - This might be able to occur sooner, during generateChunkPoints
	for i, p := range allPoints {
		if p.X >= chunk.MinX && p.X < chunk.MaxX && p.Y >= chunk.MinZ && p.Y < chunk.MaxZ {
			coreIndices = append(coreIndices, core.SiteIndex(i))
		}
	}
	chunk.CoreSiteIndices = coreIndices

	// Build sites with heights
	sites := make([]core.Site, len(allPoints))
	heights := make([]float64, len(allPoints))
	np := cm.noiseParams
	for i, p := range allPoints {
		// This includes the halo region
		h := 0.0
		if cm.heightFunc != nil {
			h = cm.heightFunc(p.X, p.Y, np.Octaves, np.Frequency, np.Amplitude, np.Persistence, np.Lacunarity)
		}
		sites[i] = core.Site{Pos: p, Height: h}
		heights[i] = h
	}
	chunk.Heights = heights

	// Build Delaunay triangulation
	tris := core.Triangulate(allPoints)
	mesh, err := core.BuildHalfEdgeMesh(sites, tris)
	if err != nil {
		return nil, err
	}
	chunk.Mesh = mesh

	// Compute face normals
	chunk.FaceNormals = mesh.AllFaceNormals(heights)

	// Build spatial index
	// Cell size roughly equal to minimum point distance for good performance
	chunk.Spatial = BuildSpatialGrid(mesh, heights, cm.cfg.MinPointDist,
		chunk.MinX-cm.cfg.HaloWidth, chunk.MinZ-cm.cfg.HaloWidth,
		chunk.MaxX+cm.cfg.HaloWidth, chunk.MaxZ+cm.cfg.HaloWidth)

	// chunk.Erosion = cm.GenerateChunkErosion(chunk)

	// Initialize watershed weights to 0.0 (no watershed modification)
	chunk.WatershedWeights = make([]float64, len(allPoints))

	// Generate erosion height deltas on the erosion manager for this chunk
	cm.erosMgr.HydraulicErosion(chunk, chunk.Spatial)

	// Generate hydrology data
	chunk.Hydro = cm.generateChunkHydrology(chunk)

	return chunk, nil
}

// generateChunkPoints generates blue noise points for a chunk including halo regions.
// Uses the point cache to avoid redundant blue noise generation for neighbors.
// Caller must hold a write lock on cm.mu.
func (cm *ChunkManager) generateChunkPoints(coord core.ChunkCoord) []core.Vec2 {

	// core.ChunkCoord (0, 1)
	// e.g. 1 * 128.0 = 128.0
	minX := float64(coord.X) * cm.cfg.ChunkSize
	minZ := float64(coord.Z) * cm.cfg.ChunkSize
	maxX := float64(coord.X+1) * cm.cfg.ChunkSize
	maxZ := float64(coord.Z+1) * cm.cfg.ChunkSize

	// fmt.Printf(
	// 	"[2] info-------------\n ChunkSize\t%v\ncoord.X\t%v\ncoord.Z\t%v\nmin:\t(%v, %v)\nmax:(%v, %v)\n---------------\n", 
	// 	cm.cfg.ChunkSize, coord.X, coord.Z, minX, minZ, maxX, maxZ)

	halo := cm.cfg.HaloWidth

	// We'll collect points from multiple regions
	allPoints := make([]core.Vec2, 0, 1024) // 1024 pre-allocated points
	seen := make(map[uint64]struct{}, 1024) // 1024 pre-allocated seen points

	// Hash a point to detect duplicates (within tolerance)
	hashPoint := func(p core.Vec2) uint64 {
		// Quantize to half the min distance for dedup
		scale := 2.0 / cm.cfg.MinPointDist
		qx := int64(p.X * scale)
		qz := int64(p.Y * scale)

		// Create a bitmask for this point
		return uint64(qx)<<32 | uint64(qz)&0xFFFFFFFF
	}

	addPoint := func(p core.Vec2) bool {
		h := hashPoint(p)
		if _, exists := seen[h]; !exists {
			seen[h] = struct{}{}
			allPoints = append(allPoints, p)
			return true
		}
		return false
	}

	addPoints := func(pts []core.Vec2) {
		for _, p := range pts {
			addPoint(p)
		}
	}

	// 1. Get or generate core points for this chunk (uses cache if available)
	corePoints := cm.getOrGeneratePoints(coord)
	addPoints(corePoints)

	// 2. Plot the neighboring chunk coordinates for this chunk
	neighbors := []core.ChunkCoord{
		{X: coord.X - 1, Z: coord.Z - 1}, {X: coord.X, Z: coord.Z - 1}, {X: coord.X + 1, Z: coord.Z - 1},
		{X: coord.X - 1, Z: coord.Z}, {X: coord.X + 1, Z: coord.Z},
		{X: coord.X - 1, Z: coord.Z + 1}, {X: coord.X, Z: coord.Z + 1}, {X: coord.X + 1, Z: coord.Z + 1},
	}

	// 3 Generate all neighboring chunk points and add them to the allPoints list
	for _, neighbor := range neighbors {
		// The boundary region is where this chunk's halo overlaps the neighbor (world-space coordinates)
		nMinX := float64(neighbor.X) * cm.cfg.ChunkSize
		nMinZ := float64(neighbor.Z) * cm.cfg.ChunkSize
		nMaxX := float64(neighbor.X+1) * cm.cfg.ChunkSize
		nMaxZ := float64(neighbor.Z+1) * cm.cfg.ChunkSize

		// Compute the overlap region between our extended bounds and neighbor's core
		overlapMinX := max(minX-halo, nMinX) // Hmm... I don't think nMinX will ever be less than minX-halo
		overlapMinZ := max(minZ-halo, nMinZ)
		overlapMaxX := min(maxX+halo, nMaxX)
		overlapMaxZ := min(maxZ+halo, nMaxZ)

		// Big ol' wtf right here
		if overlapMinX >= overlapMaxX || overlapMinZ >= overlapMaxZ {
			continue // No overlap
		}

		// Get neighbor's core points (from chunk cache, points cache, or generate + cache)
		coreNeighborPoints := cm.getOrGeneratePoints(neighbor)

		/*
			[Optimization]
			Definitely a spot that can be optimized by caching the padded boundary within each chunk
			that represents the halo region of neighboring chunks
		*/
		// Filter to points within our halo region but outside our core
		for _, p := range coreNeighborPoints {
			inHalo := (p.X >= minX-halo && p.X < maxX+halo && p.Y >= minZ-halo && p.Y < maxZ+halo)
			inCore := (p.X >= minX && p.X < maxX && p.Y >= minZ && p.Y < maxZ)
			if inHalo && !inCore {
				addPoint(p) // Add neighboring blue-noise to boundary points for the current chunk being processed
			}
		}
	}

	return allPoints
}

// getOrGeneratePoints returns blue noise points for a chunk, using caches when available.
// Priority: 1) Points Cache, 2) Chunk Cache (core sites), 3) Generate new.
// Caller must hold at least a read lock on cm.mu.
// If generation is needed, caller should hold a write lock.
func (cm *ChunkManager) getOrGeneratePoints(coord core.ChunkCoord) []core.Vec2 {

	// Check if points are cached.
	if pts, ok := cm.pointsCache[coord]; ok {
		return pts
	}

	// Optional sanity check. Check the chunk directly for the points
	// This will never happen under normal circumstances. Caching neighboring points affords us some invariants
	if chunk, ok := cm.cache[coord]; ok {
		points := make([]core.Vec2, len(chunk.Mesh.Sites))
		for i, site := range chunk.Mesh.Sites {
			points[i] = site.Pos
		}

		return points
	}

	// Generate blue noise and cache it
	minX := float64(coord.X) * cm.cfg.ChunkSize
	minZ := float64(coord.Z) * cm.cfg.ChunkSize
	maxX := float64(coord.X+1) * cm.cfg.ChunkSize
	maxZ := float64(coord.Z+1) * cm.cfg.ChunkSize

	blueCfg := DefaultBlueNoiseConfig(cm.cfg.MinPointDist)
	chunkseed := chunkSeed(cm.cfg.WorldSeed, coord)
	pts := GenerateBlueNoiseSeeded(chunkseed, minX, minZ, maxX, maxZ, blueCfg)

	// Cache the seed for this chunk
	cm.cache[coord].Seed = chunkseed // Cache on terrain chunk - more likely to be destroyed
	cm.chunkSeeds[coord] = chunkseed // Cache on chunk manager

	// Cache the generated points
	cm.pointsCache[coord] = pts
	return pts
}


// generateChunkHydrology computes rivers, lakes, and ocean data for a chunk.
func (cm *ChunkManager) generateChunkHydrology(chunk *core.TerrainChunk) *core.ChunkHydroData {
	hydro := &core.ChunkHydroData{
		Rivers:       make([]core.RiverSegment, 0),
		Lakes:        make([]core.Lake, 0),
		RiverEntries: make([]core.RiverInterPoint, 0),
		RiverExits:   make([]core.RiverInterPoint, 0),
		LakeEntries:  make([]core.LakeInterPoint, 0),
		LakeExits:    make([]core.LakeInterPoint, 0),
	}

	/**
		Might be wise to organize the initialization of the hydro data.
		- Hydrology doesn't follow the phased architecture we've been supporting
			Though it expresses it.
	 */

	// Print a bunch of "-"
	fmt.Println("--------------------------------")

	// Bounds in world-space
	chunkBounds :=  core.ChunkBoundaryBox{
		MinX: chunk.MinX, 
		MinZ: chunk.MinZ, 
		MaxX: chunk.MaxX, 
		MaxZ: chunk.MaxZ,
	}

	// Reset visited sites tracking for this chunk's hydrology pass
	cm.hydroMgr.ResetVisitedSites()

	// 1. Process incoming rivers from neighbors - Never occurs for the first chunk
	entries := cm.hydroMgr.GetRiverEntries(chunk.Coord)
	for _, entry := range entries {
		hydro.RiverEntries = append(hydro.RiverEntries, entry)

		// Find the nearest site to the entry point
		nearestSite := cm.findNearestSite(chunk.Mesh, entry.Position)

		//fmt.Println("[ENTRY] nearestSite", nearestSite)
		if nearestSite < 0 || cm.hydroMgr.IsSourceVisited(nearestSite) {
			//fmt.Println("[ENTRY] skipping visited entry", nearestSite)
			continue
		}

		// Continue tracing the river
		segment, exit := cm.hydroMgr.TraceRiver(
			chunk.Mesh, 
			chunk.Heights, 
			&cm.hydroMgr.VisitedSites, 
			nearestSite, 
			chunkBounds,
			entry.Distance,
		)

		if len(segment.Vertices) > 0 {
			hydro.Rivers = append(hydro.Rivers, segment)
		}

		if exit != nil {
			exit.RiverID = entry.RiverID
			hydro.RiverExits = append(hydro.RiverExits, *exit)
			cm.hydroMgr.RegisterRiverExit(chunk.Coord, *exit)
		}
	}
	cm.hydroMgr.ClearRiverEntries(chunk.Coord)

	// 2. Find new river sources in this chunk
	sources := cm.hydroMgr.FindRiverSources(chunk.Mesh, chunk.Heights, chunk.CoreSiteIndices)
	for _, source := range sources {
		// Skip sources that overlap with already-traced river paths
		if cm.hydroMgr.IsSourceVisited(source) {
			//fmt.Println("[SRC] skipping visited source", source)
			continue
		}

		segment, exit := cm.hydroMgr.TraceRiverV2(
			chunk.Mesh, 
			chunk.Heights, 
			&cm.hydroMgr.VisitedSites, 
			source, 
			chunkBounds,
			 0,
		)

		if len(segment.Vertices) > 0 {
			hydro.Rivers = append(hydro.Rivers, segment)
		}

		if exit != nil {
			exit.RiverID = cm.hydroMgr.NextRiverID
			cm.hydroMgr.NextRiverID++
			hydro.RiverExits = append(hydro.RiverExits, *exit)
			cm.hydroMgr.RegisterRiverExit(chunk.Coord, *exit)
		}
	}

	debug := make(map[core.Vec2]int)
	for _, river := range hydro.Rivers {
		for _, vertex := range river.Vertices {
			debug[vertex]++
		}
	}

	// // Print results of debug
	// for x, count := range debug {
	// 	fmt.Printf("(%f, %f): %d\n", x.X, x.Y, count)
	// }

	// 3. Process incoming lake boundaries from neighbors FIRST
	chunkLakeEntries := cm.hydroMgr.GetLakeEntries(chunk.Coord)

	// #region agent log
	debugLogPath := `c:\Users\Yoshi\go\src\procedural_generation\terrain_generation\.cursor\debug.log`
	chunkLogEntry := func(msg string, data map[string]interface{}) {
		f, err := os.OpenFile(debugLogPath, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil { return }
		defer f.Close()
		data["location"] = "chunk.go:generateChunkHydrology"
		data["message"] = msg
		data["timestamp"] = time.Now().UnixMilli()
		data["sessionId"] = "debug-session"
		data["hypothesisId"] = "B,C,E"
		jsonBytes, _ := json.Marshal(data)
		f.WriteString(string(jsonBytes) + "\n")
	}
	chunkLogEntry("HYDRO_START", map[string]interface{}{
		"chunk": fmt.Sprintf("(%d,%d)", chunk.Coord.X, chunk.Coord.Z),
		"lakeEntriesCount": len(chunkLakeEntries),
		"totalLakesInMgr": len(cm.hydroMgr.Lakes),
	})
	// #endregion

	for _, entry := range chunkLakeEntries {
		hydro.LakeEntries = append(hydro.LakeEntries, entry)

		// #region agent lo
		chunkLogEntry("LAKE_ENTRY_PROCESSING", map[string]interface{}{
			"entryLakeID": entry.LakeID, "entryWaterLevel": entry.WaterLevel,
			"entrySideIndex": entry.Side, "hypothesisJ": len(chunkLakeEntries) > 20,
		})
		// #endregion

		// Find the nearest site to the entry point
		// PROBLEM?: findNearestSite searches the entire mesh for the nearest site.
		// We built a half-edge mesh, so we should be able to find the nearest site using the twin->next->twin->next... chain.
		nearestSite := cm.findNearestSite(chunk.Mesh, entry.Position)
		if nearestSite < 0 {
			continue
		}

		// Check if this site is already part of a lake.
		existingLakeID := cm.hydroMgr.GetLakeForSite(chunk.Coord, nearestSite)

		// FIX: Skip if site is already part of THIS SAME lake (prevents circular propagation)
		if existingLakeID == entry.LakeID {
			continue
		}

		if existingLakeID != -1 && existingLakeID != entry.LakeID {
			// Two different lakes are meeting - handle merge or watershed
			existingLake := cm.hydroMgr.Lakes[existingLakeID]
			entryLake := cm.hydroMgr.Lakes[entry.LakeID]
			if existingLake != nil && entryLake != nil {
				sharedSites := cm.hydroMgr.FindSharedBoundarySites(
					chunk.Mesh,
					existingLake.SiteIndices,
					entryLake.SiteIndices,
				)
				shouldMerge, _ := cm.hydroMgr.HandleLakeOverlap(
					existingLakeID,
					entry.LakeID,
					sharedSites,
					chunk.WatershedWeights,
				)
				if !shouldMerge {
					// Watershed divide was created, skip continuing this lake
					continue
				}
			}
		}

		// Continue flood fill from this entry point
		lake, exits := cm.hydroMgr.FloodFillLakeWithBoundaries(
			chunk.Mesh,
			chunk.Heights,
			core.SiteIndex(nearestSite),
			chunkBounds, // Note: at this point this is simply the chunk's x/z plane
			chunk.Coord,
			entry.LakeID,
		)

		if lake != nil {
			hydro.Lakes = append(hydro.Lakes, *lake)
			cm.hydroMgr.RegisterLakeInChunk(lake.ID, chunk.Coord)

			// Register boundary exits for propagation to neighbors
			for _, exit := range exits {
				destChunk := cm.getNeighborChunkForEdge(chunk.Coord, exit.Side)
				cm.hydroMgr.RegisterLakeBoundary(chunk.Coord, destChunk, exit)
				hydro.LakeExits = append(hydro.LakeExits, exit)
			}
		}
	}
	cm.hydroMgr.ClearLakeEntries(chunk.Coord)

	// Detect NEW lakes at local minima (excluding sites already in lakes)
	// `minima` is a list of sites with lower height than all their neighbors
	minima := core.FindLocalMinima(chunk.Mesh, chunk.Heights)

	for _, minimum := range minima {
		// Skip if this site is already part of a lake
		// NOTE: This could happen ahead of time. Filter minima beforehand, though not terrible
		// NOTE: This implies a lake was continued from a neighbor and we either merged or watershed
		if cm.hydroMgr.GetLakeForSite(chunk.Coord, minimum) != -1 {
			continue
		}
 
		// ?? Why is it problematic if the site isn't in the core region?
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

		lake, exits := cm.hydroMgr.FloodFillLakeWithBoundaries(
			chunk.Mesh,
			chunk.Heights,
			minimum,
			chunkBounds,
			chunk.Coord,
			-1, // New lake
		)
		if lake != nil {
			hydro.Lakes = append(hydro.Lakes, *lake)
			cm.hydroMgr.RegisterLakeInChunk(lake.ID, chunk.Coord)

			// #region agent log
			chunkLogEntry("NEW_LAKE_CREATED", map[string]interface{}{
				"lakeID": lake.ID, "sitesCount": len(lake.SiteIndices), "exitsCount": len(exits),
			})
			// #endregion

			// Register boundary exits for propagation to neighbors
			for _, exit := range exits {
				destChunk := cm.getNeighborChunkForEdge(chunk.Coord, exit.Side)
				cm.hydroMgr.RegisterLakeBoundary(chunk.Coord, destChunk, exit)
				hydro.LakeExits = append(hydro.LakeExits, exit)
			}
		}
	}

	// 5. Compute ocean data if this chunk is in an ocean region
	hydro.Ocean = cm.hydroMgr.ComputeOceanChunkData(chunk.Coord, chunk.Mesh, chunk.Heights)

	return hydro
}

// getNeighborChunkForEdge returns the neighbor chunk in the direction of the given edge.
// Edge indices: 0=minX, 1=maxX, 2=minZ, 3=maxZ
func (cm *ChunkManager) getNeighborChunkForEdge(chunk core.ChunkCoord, side core.SideIndex) core.ChunkCoord {
	switch side {
	case 0: // minX edge -> neighbor to the west
		return core.ChunkCoord{X: chunk.X - 1, Z: chunk.Z}
	case 1: // maxX edge -> neighbor to the east
		return core.ChunkCoord{X: chunk.X + 1, Z: chunk.Z}
	case 2: // minZ edge -> neighbor to the south
		return core.ChunkCoord{X: chunk.X, Z: chunk.Z - 1}
	case 3: // maxZ edge -> neighbor to the north
		return core.ChunkCoord{X: chunk.X, Z: chunk.Z + 1}
	default:
		return chunk // Invalid edge, return self
	}
}

// findNearestSite finds the mesh site closest to a given position.
func (cm *ChunkManager) findNearestSite(mesh *core.DelaunayMesh, pos core.Vec2) core.SiteIndex {
	bestDist := float64(1e18)
	bestSite := core.SiteIndex(-1)

	for i, site := range mesh.Sites {
		d := site.Pos.Sub(pos).Len2()
		if d < bestDist {
			bestDist = d
			bestSite = core.SiteIndex(i)
		}
	}

	return bestSite
}

// chunkSeed computes a deterministic seed for a chunk based on world seed and coordinates.
func chunkSeed(worldSeed int64, coord core.ChunkCoord) int64 {
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

// // SampleHeight returns the interpolated height at position (x, z).
// // Returns (height, ok) where ok=false if the point is outside the mesh.
// func (c *TerrainChunk) SampleHeight(x, z float64) (float64, bool) {
// 	if c.Spatial == nil {
// 		return 0, false
// 	}
// 	return c.Spatial.SampleHeight(x, z)
// }

// // Helper: seeded RNG for boundary regions
// func boundaryRNG(worldSeed int64, coord1, coord2 core.ChunkCoord) *rand.Rand {
// 	seed := boundarySeed(worldSeed, coord1, coord2)
// 	return rand.New(rand.NewSource(seed))
// }
// boundarySeed computes a deterministic seed for the shared boundary between two chunks.
// It uses the minimum of the two chunk coords to ensure both chunks get the same seed.

// func boundarySeed(worldSeed int64, coord1, coord2 core.ChunkCoord) int64 {
// 	// Use lexicographically smaller coord first
// 	var first, second core.ChunkCoord
// 	if coord1.X < coord2.X || (coord1.X == coord2.X && coord1.Z < coord2.Z) {
// 		first, second = coord1, coord2
// 	} else {
// 		first, second = coord2, coord1
// 	}

// 	h := fnv.New64a()
// 	buf := make([]byte, 8)

// 	// World seed
// 	buf[0] = byte(worldSeed)
// 	buf[1] = byte(worldSeed >> 8)
// 	buf[2] = byte(worldSeed >> 16)
// 	buf[3] = byte(worldSeed >> 24)
// 	buf[4] = byte(worldSeed >> 32)
// 	buf[5] = byte(worldSeed >> 40)
// 	buf[6] = byte(worldSeed >> 48)
// 	buf[7] = byte(worldSeed >> 56)
// 	h.Write(buf)

// 	// First chunk
// 	buf[0] = byte(first.X)
// 	buf[1] = byte(first.X >> 8)
// 	buf[2] = byte(first.Z)
// 	buf[3] = byte(first.Z >> 8)
// 	// Second chunk
// 	buf[4] = byte(second.X)
// 	buf[5] = byte(second.X >> 8)
// 	buf[6] = byte(second.Z)
// 	buf[7] = byte(second.Z >> 8)
// 	h.Write(buf)

// 	return int64(h.Sum64())
// }

// CacheSize returns the number of cached chunks.
// func (cm *ChunkManager) CacheSize() int {
// 	cm.mu.RLock()
// 	defer cm.mu.RUnlock()
// 	return len(cm.cache)
// }

// // CachedCoords returns all cached chunk coordinates.
// func (cm *ChunkManager) CachedCoords() []core.ChunkCoord {
// 	cm.mu.RLock()
// 	defer cm.mu.RUnlock()
// 	coords := make([]core.ChunkCoord, 0, len(cm.cache))
// 	for coord := range cm.cache {
// 		coords = append(coords, coord)
// 	}
// 	return coords
// }
