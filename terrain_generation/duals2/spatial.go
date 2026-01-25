package duals2

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// BuildSpatialGrid constructs a spatial grid for fast triangle lookup.
// cellSize should be roughly the average triangle edge length for best performance.
func BuildSpatialGrid(mesh *core.DelaunayMesh, heights []float64, cellSize, minX, minZ, maxX, maxZ float64) *core.SpatialGrid {
	if mesh == nil || len(mesh.Tris) == 0 {
		return nil
	}

	width := maxX - minX
	height := maxZ - minZ
	if width <= 0 || height <= 0 || cellSize <= 0 {
		return nil
	}

	gridW := int(math.Ceil(width / cellSize))
	gridH := int(math.Ceil(height / cellSize))
	if gridW <= 0 {
		gridW = 1
	}
	if gridH <= 0 {
		gridH = 1
	}

	sg := &core.SpatialGrid{
		Mesh:     mesh,
		Heights:  heights,
		CellSize: cellSize, // make sure SpatialGrid.CellSize is tuned appropriately
		// A good heuristic is 2-4x your average triangle edge length, but it depends on your query patterns.
		MinX:     minX,
		MinZ:     minZ,
		MaxX:     maxX,
		MaxZ:     maxZ,
		GridW:    gridW,
		GridH:    gridH,
		Cells:    make([][]int, gridW*gridH),
	}

	// For each triangle, find which cells it overlaps and add it to those cells
	for ti, t := range mesh.Tris {
		// Get triangle vertices
		a := mesh.Sites[t.A].Pos
		b := mesh.Sites[t.B].Pos
		c := mesh.Sites[t.C].Pos

		// Find bounding box of triangle
		triMinX := min(a.X, min(b.X, c.X))
		triMaxX := max(a.X, max(b.X, c.X))
		triMinZ := min(a.Y, min(b.Y, c.Y))
		triMaxZ := max(a.Y, max(b.Y, c.Y))

		// Convert to grid cells
		cellMinX := int((triMinX - minX) / cellSize)
		cellMaxX := int((triMaxX - minX) / cellSize)
		cellMinZ := int((triMinZ - minZ) / cellSize)
		cellMaxZ := int((triMaxZ - minZ) / cellSize)

		// Clamp to grid bounds
		if cellMinX < 0 {
			cellMinX = 0
		}
		if cellMaxX >= gridW {
			cellMaxX = gridW - 1
		}
		if cellMinZ < 0 {
			cellMinZ = 0
		}
		if cellMaxZ >= gridH {
			cellMaxZ = gridH - 1
		}

		// Add triangle to all overlapping cells
		for cz := cellMinZ; cz <= cellMaxZ; cz++ {
			for cx := cellMinX; cx <= cellMaxX; cx++ {
				cellIdx := cz*gridW + cx
				sg.Cells[cellIdx] = append(sg.Cells[cellIdx], ti)
			}
		}
	}

	return sg
}
