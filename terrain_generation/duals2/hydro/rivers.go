package hydro

import "procedural_generation/terrain_generation/duals2/core"

// -------------------------------------------------------------------
// River Source Detection
// -------------------------------------------------------------------

// FindRiverSources identifies potential river sources in a chunk.
// Sources are considered at each Site (a triangle vertex)
// TODO: Consider the circumcenter that we're already plotting for the Voronoi diagram
// Sources are high-elevation vertices that pass the probability check.
func (hm *HydroManager) FindRiverSources(
	mesh *core.DelaunayMesh,
	heights []float64,
	coreSiteIndices []core.SiteIndex,
) []core.SiteIndex {
	sources := make([]core.SiteIndex, 0)

	candidateCount := 0
	for _, site := range coreSiteIndices {
		h := heights[site]

		// Check elevation threshold
		if h < hm.Cfg.RiverSourceMinElev || h > hm.Cfg.RiverSourceMaxElev {
			continue
		}

		// TODO - This prevents sources in basins. See docs
		// Check if this is a local maximum or near-maximum
		// (sources shouldn't be in valleys)
		if !core.IsNearLocalMaximum(mesh, heights, site) {
			continue
		}

		candidateCount++
		randVal := hm.Rng.Float64()
		// Apply probability
		if randVal > hm.Cfg.RiverSourceProbability {
			continue
		}

		sources = append(sources, site)
	}

	return sources
}

// TraceRiver traces a river path from a starting position using gradient descent.
// Returns the river segment and any exit point if the river leaves the chunk bounds.
func (hm *HydroManager) TraceRiver(
	mesh *core.DelaunayMesh,
	heights []float64,
	visitedSites *map[core.SiteIndex]bool,
	startSite core.SiteIndex,	// source - A Site/vertex
	chunkBounds core.ChunkBoundaryBox,
	distanceTraveled float64,
) (core.RiverSegment, *core.RiverInterPoint) {

	segment := core.RiverSegment{
		Vertices: make([]core.Vec2, 0, 64),
		Widths:   make([]float64, 0, 64),
		Depths:   make([]float64, 0, 64),
	}

	current := startSite
	visited := make(map[core.SiteIndex]bool)
	distance := distanceTraveled
	width := hm.Cfg.RiverWidthBase
	depth := hm.Cfg.RiverDepthBase

	var flowDirAccum core.Vec2
	var exitPoint *core.RiverInterPoint

	for iterations := 0; iterations < 10000; iterations++ {
		if visited[current] {
			break // Cycle detected, stop
		}
		visited[current] = true
		if hm.VisitedSites != nil {
			hm.VisitedSites[current] = true
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

			exitPoint = &core.RiverInterPoint{
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
		neighbors := core.GetLowerNeighbors(mesh, hm.Cfg.RiverMinSlope, heights, current)
		if len(neighbors) == 0 {
			// Local minimum - potential lake site
			break
		}

		// Apply flow bias for tie-breaking
		bias := hm.Cfg.FlowBias(pos.Y) // pos.Y is Z in world space
		next := hm.selectNextVertexV1(mesh, heights, current, neighbors, h, bias)

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
		width = hm.Cfg.RiverWidthBase + distance*hm.Cfg.RiverWidthGrowth
		depth = hm.Cfg.RiverDepthBase + distance*hm.Cfg.RiverDepthGrowth

		current = next
	}

	// Compute average flow direction
	if len(segment.Vertices) > 1 {
		segment.FlowDir = flowDirAccum.Normalize()
	}

	return segment, exitPoint
}

func (hm *HydroManager) TraceRiverV2(
	mesh *core.DelaunayMesh,
	heights []float64,
	visitedSites *map[core.SiteIndex]bool,
	startSite core.SiteIndex,	// source - A Site/vertex
	chunkBounds core.ChunkBoundaryBox,
	distanceTraveled float64,
) (core.RiverSegment, *core.RiverInterPoint) {

	segment := core.RiverSegment{
		Vertices: make([]core.Vec2, 0, 64),
		Widths:   make([]float64, 0, 64),
		Depths:   make([]float64, 0, 64),
	}

	current := startSite
	visited := make(map[core.SiteIndex]bool)
	distance := distanceTraveled
	width := hm.Cfg.RiverWidthBase
	depth := hm.Cfg.RiverDepthBase

	var flowDirAccum core.Vec2
	var exitPoint *core.RiverInterPoint

	// Note: Iterations is an arbitrary limit
	for iterations := 0; iterations < 10000; iterations++ {
		if visited[current] {
			break // Cycle detected, stop
		}
		visited[current] = true
		if hm.VisitedSites != nil {
			hm.VisitedSites[current] = true
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

			exitPoint = &core.RiverInterPoint{
				Position: pos,
				Width:    width,
				Depth:    depth,
				FlowDir:  flowDirAccum.Normalize(),
				Distance: distance,
				ExitEdge: exitEdge,
			}
			break
		}

		bias := hm.Cfg.FlowBias(pos.Y) // pos.Y is Z in world space

		// Find neighbors with lower elevation
		neighbors := core.GetFlowBiasedNeighbors(mesh, current, bias)
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
		width = hm.Cfg.RiverWidthBase + distance*hm.Cfg.RiverWidthGrowth
		depth = hm.Cfg.RiverDepthBase + distance*hm.Cfg.RiverDepthGrowth

		current = next
	}

	// Compute average flow direction
	if len(segment.Vertices) > 1 {
		segment.FlowDir = flowDirAccum.Normalize()
	}

	return segment, exitPoint

}


// -------------------------------------------------------------------
// Cross-Chunk River Propagation
// -------------------------------------------------------------------

// RegisterRiverExit records a river exiting a chunk for propagation.
func (hm *HydroManager) RegisterRiverExit(fromChunk core.ChunkCoord, exit core.RiverInterPoint) {
	// Determine destination chunk based on exit edge
	var destChunk core.ChunkCoord
	switch exit.ExitEdge {
	case 0: // minX
		destChunk = core.ChunkCoord{X: fromChunk.X - 1, Z: fromChunk.Z}
	case 1: // maxX
		destChunk = core.ChunkCoord{X: fromChunk.X + 1, Z: fromChunk.Z}
	case 2: // minZ
		destChunk = core.ChunkCoord{X: fromChunk.X, Z: fromChunk.Z - 1}
	case 3: // maxZ
		destChunk = core.ChunkCoord{X: fromChunk.X, Z:fromChunk.Z + 1}
	}

	hm.RiverPropagation[destChunk] = append(hm.RiverPropagation[destChunk], exit)
}

// GetRiverEntries returns rivers that should enter this chunk from neighbors.
func (hm *HydroManager) GetRiverEntries(chunk core.ChunkCoord) []core.RiverInterPoint {
	return hm.RiverPropagation[chunk]
}

// ClearRiverEntries removes processed river entries for a chunk.
func (hm *HydroManager) ClearRiverEntries(chunk core.ChunkCoord) {
	delete(hm.RiverPropagation, chunk)
}
