package hydro

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// FloodFillLake fills a depression from a local minimum to find lake extent.
// Returns the lake if valid, nil otherwise.
func (hm *HydroManager) FloodFillLake(
	mesh *core.DelaunayMesh,
	heights []float64,
	minimum int,
) *core.Lake {
	// Check probability - if 0, never create lakes
	if hm.Rng.Float64() >= hm.Cfg.LakeProbability {
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
					if destHeight <= minHeight+hm.Cfg.LakeMinDepth*2 {
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
	if len(filled) < hm.Cfg.LakeMinArea {
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

	lake := &core.Lake{
		ID:          hm.NextLakeID,
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
	}

	hm.NextLakeID++
	return lake
}