package hydro

import (
	"fmt"
	"procedural_generation/terrain_generation/duals2/core"
)

// isNearLocalMaximum returns true if the site is higher than most neighbors.
func isNearLocalMaximum(mesh *core.DelaunayMesh, heights []float64, site int) bool {
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