package hydro

import (
	"fmt"
	"procedural_generation/terrain_generation/duals2/core"
)

func getFlowBiasedNeighbors(mesh *core.DelaunayMesh, site int, bias core.Vec2) []int {

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