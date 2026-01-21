package hydro

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// FloodFillLake fills a depression from a local minimum to find lake extent.
// Returns the lake if valid, nil otherwise.
// This is the original simple version for backward compatibility.
func (hm *HydroManager) FloodFillLake(
	mesh *core.DelaunayMesh,
	heights []float64,
	minimum int,
) *core.Lake {
	lake, _ := hm.FloodFillLakeWithBoundaries(mesh, heights, minimum, core.ChunkBounds{}, core.ChunkCoord{}, -1)
	return lake
}

// FloodFillLakeWithBoundaries fills a depression and detects boundary crossings.
// chunkBounds: the core bounds of the current chunk (for boundary detection)
// sourceChunk: the current chunk's coordinates
// existingLakeID: -1 for new lake, or lake ID when continuing from neighbor
// Returns the lake and any boundary points where the lake exits to neighbors.
func (hm *HydroManager) FloodFillLakeWithBoundaries(
	mesh *core.DelaunayMesh,
	heights []float64,
	minimum int,
	chunkBounds core.ChunkBounds,
	sourceChunk core.ChunkCoord,
	existingLakeID int,
) (*core.Lake, []core.LakeBoundaryPoint) {
	// Check probability - if 0, never create lakes (only for new lakes)
	if existingLakeID == -1 && hm.Rng.Float64() >= hm.Cfg.LakeProbability {
		return nil, nil
	}

	minHeight := heights[minimum]

	// Find all connected vertices below the spillway level
	filled := make(map[int]bool)
	frontier := []int{minimum}
	filled[minimum] = true

	var spillway *int
	spillwayHeight := math.Inf(1)

	// If continuing an existing lake, use its water level as initial spillway
	if existingLakeID != -1 {
		if existingLake, exists := hm.Lakes[existingLakeID]; exists {
			spillwayHeight = existingLake.WaterLevel
		}
	}

	// Track boundary sites for propagation
	boundarySites := make(map[int]int) // siteIndex -> edgeIndex (0=minX, 1=maxX, 2=minZ, 3=maxZ)

	// Threshold for boundary detection (slightly inside the core bounds)
	const boundaryThreshold = 1.0

	for len(frontier) > 0 {
		current := frontier[0]
		frontier = frontier[1:]

		// Check if current site is near a chunk boundary
		if chunkBounds.MaxX > chunkBounds.MinX { // bounds are set
			pos := mesh.Sites[current].Pos
			if pos.X <= chunkBounds.MinX+boundaryThreshold {
				boundarySites[current] = 0 // minX edge
			} else if pos.X >= chunkBounds.MaxX-boundaryThreshold {
				boundarySites[current] = 1 // maxX edge
			}
			if pos.Y <= chunkBounds.MinZ+boundaryThreshold {
				boundarySites[current] = 2 // minZ edge
			} else if pos.Y >= chunkBounds.MaxZ-boundaryThreshold {
				boundarySites[current] = 3 // maxZ edge
			}
		}

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
				// For NEW lakes: use depression threshold AND limit maximum size
				// For CONTINUED lakes: use existing water level as threshold
				// This allows continued lakes to fill properly across chunk boundaries
				// while preventing new lakes from consuming entire chunks
				
				canFill := false
				var threshold float64
				
				if existingLakeID != -1 {
					// Continuing existing lake: use the established water level as threshold
					// Any site below the water level that is connected should be filled
					threshold = spillwayHeight
					canFill = true // Already checked destHeight < spillwayHeight above
				} else {
					// New lake: use depression threshold from local minimum
					// Also enforce a maximum size to prevent chunk-filling lakes
					threshold = minHeight + hm.Cfg.LakeMinDepth*2
					maxNewLakeSize := len(mesh.Sites) / 3 // Max 1/3 of chunk for new lakes
					if destHeight <= threshold && len(filled) < maxNewLakeSize {
						canFill = true
					}
				}
				
				if canFill {
					filled[dest] = true
					frontier = append(frontier, dest)
				} else {
					// Potential spillway - this site blocks expansion
					if destHeight < spillwayHeight {
						spillwayHeight = destHeight
						spillway = &dest
					}
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
	if len(filled) < hm.Cfg.LakeMinArea {
		return nil, nil
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

	// Determine lake ID
	lakeID := existingLakeID
	if lakeID == -1 {
		lakeID = hm.NextLakeID
		hm.NextLakeID++
	}

	// Water level reconciliation: if continuing an existing lake and found a lower spillway,
	// update the global lake's water level
	if existingLakeID != -1 {
		if existingLake, exists := hm.Lakes[existingLakeID]; exists {
			if spillwayHeight < existingLake.WaterLevel {
				// Found a lower spillway in this chunk - update global lake
				hm.UpdateLakeWaterLevel(existingLakeID, spillwayHeight)
			}
		}
	}

	lake := &core.Lake{
		ID:          lakeID,
		CenterPos:   core.Vec2{X: centerX, Y: centerZ},
		SiteIndices: siteIndices,
		WaterLevel:  spillwayHeight,
		MinDepth:    maxDepth,
	}

	if spillway != nil {
		sp := mesh.Sites[*spillway].Pos
		lake.Spillway = &sp

		// Calculate spillway direction (away from lake center)
		lake.SpillwayDir = sp.Sub(lake.CenterPos).Normalize()

		// Also update the global lake's spillway if this is a continuation
		if existingLakeID != -1 {
			if existingLake, exists := hm.Lakes[existingLakeID]; exists {
				if spillwayHeight <= existingLake.WaterLevel {
					existingLake.Spillway = &sp
					existingLake.SpillwayDir = lake.SpillwayDir
				}
			}
		}
	}

	// Create boundary points for sites near chunk edges
	var boundaryPoints []core.LakeBoundaryPoint
	for siteIdx, edgeIdx := range boundarySites {
		if filled[siteIdx] {
			pos := mesh.Sites[siteIdx].Pos
			boundaryPoints = append(boundaryPoints, core.LakeBoundaryPoint{
				LakeID:      lakeID,
				WaterLevel:  spillwayHeight,
				Position:    pos,
				SiteIndex:   siteIdx,
				EdgeIndex:   edgeIdx,
				SourceChunk: sourceChunk,
			})
		}
	}

	// Mark all lake sites in the global tracker
	for site := range filled {
		hm.MarkLakeSite(sourceChunk, site, lakeID)
	}

	// Register lake globally if new, or update existing lake
	if existingLakeID == -1 {
		hm.Lakes[lakeID] = lake

	} else {
		// Update existing lake with new sites
		if existingLake, exists := hm.Lakes[existingLakeID]; exists {
			existingLake.SiteIndices = append(existingLake.SiteIndices, siteIndices...)
			// Update depth if this chunk has a deeper point
			if maxDepth > existingLake.MinDepth {
				existingLake.MinDepth = maxDepth
			}
		}
	}

	return lake, boundaryPoints
}

// HandleLakeOverlap determines whether to merge two overlapping lakes or create a watershed divide.
// Returns (shouldMerge, mergedIntoID).
// If shouldMerge is false, the sharedSites are marked with watershed weights.
func (hm *HydroManager) HandleLakeOverlap(
	existingLakeID int,
	newLakeID int,
	sharedSites []int,
	watershedWeights []float64,
) (shouldMerge bool, mergedIntoID int) {
	existingLake, existsExisting := hm.Lakes[existingLakeID]
	newLake, existsNew := hm.Lakes[newLakeID]

	if !existsExisting || !existsNew {
		return false, -1
	}

	levelDiff := math.Abs(existingLake.WaterLevel - newLake.WaterLevel)

	if levelDiff <= hm.Cfg.LakeWaterLevelMergeThreshold {
		// Merge: combine into the lake with the lower spillway
		if existingLake.WaterLevel <= newLake.WaterLevel {
			hm.MergeLakes(existingLakeID, newLakeID)
			return true, existingLakeID
		}
		hm.MergeLakes(newLakeID, existingLakeID)
		return true, newLakeID
	}

	// Watershed divide: mark shared boundary vertices
	for _, site := range sharedSites {
		// Weight based on how close the levels are (smoother transition)
		// Higher weight = more elevation needed
		weight := 1.0 - (levelDiff / (hm.Cfg.LakeWaterLevelMergeThreshold * 2))
		if weight < 0.5 {
			weight = 0.5 // Minimum watershed visibility
		}
		if site < len(watershedWeights) {
			watershedWeights[site] = weight
		}
	}

	return false, -1
}

// MergeLakes combines two lakes into one.
// The target lake absorbs the source lake's sites.
func (hm *HydroManager) MergeLakes(targetID, sourceID int) {
	targetLake, existsTarget := hm.Lakes[targetID]
	sourceLake, existsSource := hm.Lakes[sourceID]

	if !existsTarget || !existsSource {
		return
	}

	// Merge site indices
	targetLake.SiteIndices = append(targetLake.SiteIndices, sourceLake.SiteIndices...)

	// Note: LakeSites tracking update is skipped here because we don't have
	// chunk coordinates for historical site indices. The sites were already
	// marked when originally filled.

	// Note: Center recalculation would require mesh access.
	// This could be done by the caller if needed.

	// Use the lower water level (already done since target has lower spillway)
	if sourceLake.WaterLevel < targetLake.WaterLevel {
		targetLake.WaterLevel = sourceLake.WaterLevel
		targetLake.Spillway = sourceLake.Spillway
		targetLake.SpillwayDir = sourceLake.SpillwayDir
	}

	// Merge chunk tracking
	sourceChunks := hm.LakeChunks[sourceID]
	for _, chunk := range sourceChunks {
		found := false
		for _, existing := range hm.LakeChunks[targetID] {
			if existing == chunk {
				found = true
				break
			}
		}
		if !found {
			hm.LakeChunks[targetID] = append(hm.LakeChunks[targetID], chunk)
		}
	}

	// Remove source lake
	delete(hm.Lakes, sourceID)
	delete(hm.LakeChunks, sourceID)
}

// FindSharedBoundarySites finds sites that are at the boundary between two lakes.
// This is used to determine watershed divide locations.
// Note: Only sites valid for the given mesh are considered (lake sites may come from different chunks).
func (hm *HydroManager) FindSharedBoundarySites(
	mesh *core.DelaunayMesh,
	lake1Sites, lake2Sites []int,
) []int {
	meshSize := len(mesh.Sites)

	// Filter to only sites valid for this mesh
	lake1Set := make(map[int]bool)
	for _, s := range lake1Sites {
		if s >= 0 && s < meshSize {
			lake1Set[s] = true
		}
	}

	lake2Set := make(map[int]bool)
	for _, s := range lake2Sites {
		if s >= 0 && s < meshSize {
			lake2Set[s] = true
		}
	}

	// Find sites in lake1 that have neighbors in lake2 (and vice versa)
	sharedSites := make(map[int]bool)

	for site := range lake1Set {
		// Bounds check for SiteEdge array
		if site >= len(mesh.SiteEdge) {
			continue
		}
		startEdge := mesh.SiteEdge[site]
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

			// Bounds check for HalfEdges array
			if edge < 0 || edge >= len(mesh.HalfEdges) {
				break
			}

			dest := mesh.HalfEdges[edge].EdgeDest
			if lake2Set[dest] {
				// This site in lake1 has a neighbor in lake2
				sharedSites[site] = true
				sharedSites[dest] = true
			}

			twin := mesh.HalfEdges[edge].Twin
			if twin == -1 {
				break
			}
			if twin >= len(mesh.HalfEdges) {
				break
			}
			edge = mesh.HalfEdges[twin].Next
			if edge == startEdge {
				break
			}
		}
	}

	result := make([]int, 0, len(sharedSites))
	for site := range sharedSites {
		result = append(result, site)
	}
	return result
}