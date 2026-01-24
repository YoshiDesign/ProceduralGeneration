package core

// ChunkSiteKey uniquely identifies a site across all chunks.
// This is a general purpose key for tracking sites across chunks.
type ChunkSiteKey struct {
	Coordinates     ChunkCoord
	SiteIndex SiteIndex
}

// FindLocalMinima identifies vertices that are lower than all their neighbors.
// These are potential lake sites.
func FindLocalMinima(mesh *DelaunayMesh, heights []float64) []SiteIndex {
	minima := make([]SiteIndex, 0)

	for site := range mesh.Sites {
		if isLocalMinimum(mesh, heights, SiteIndex(site)) {
			minima = append(minima, SiteIndex(site))
		}
	}

	return minima
}

// isLocalMinimum returns true if the site is lower than all its neighbors.
func isLocalMinimum(mesh *DelaunayMesh, heights []float64, site SiteIndex) bool {
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

