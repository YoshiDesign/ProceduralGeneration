package hydro

import "procedural_generation/terrain_generation/duals2/core"

/**
 * getLowerNeighbors returns site indices of neighbors with lower elevation.
 * Neighbor is added to selections if `heights[dest] < currentHeight-cfg.RiverMinSlope`
 * Measurement begins with the current site's height.
 */
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

/**
 *	getFlowBiasedNeighbors returns site indices of neighbors aligned with the flow bias.
 *  Neighbor is added to selections if `biasAlignment > 0.1`
 */
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