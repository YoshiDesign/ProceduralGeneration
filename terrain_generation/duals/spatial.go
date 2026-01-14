package duals

import "math"

// SpatialGrid provides O(1) point location for height queries on a Delaunay mesh.
type SpatialGrid struct {
	mesh     *DelaunayMesh
	heights  []float64
	cellSize float64
	minX, minZ float64
	maxX, maxZ float64
	gridW, gridH int

	// Each cell stores a list of triangle indices that overlap it
	cells [][]int
}

// BuildSpatialIndex constructs a spatial grid for fast triangle lookup.
// cellSize should be roughly the average triangle edge length for best performance.
func BuildSpatialIndex(mesh *DelaunayMesh, heights []float64, cellSize, minX, minZ, maxX, maxZ float64) *SpatialGrid {
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

	sg := &SpatialGrid{
		mesh:     mesh,
		heights:  heights,
		cellSize: cellSize,
		minX:     minX,
		minZ:     minZ,
		maxX:     maxX,
		maxZ:     maxZ,
		gridW:    gridW,
		gridH:    gridH,
		cells:    make([][]int, gridW*gridH),
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
				sg.cells[cellIdx] = append(sg.cells[cellIdx], ti)
			}
		}
	}

	return sg
}

// LocateTriangle finds the triangle containing point (x, z).
// Returns (triangleIndex, ok) where ok=false if the point is outside all triangles.
func (sg *SpatialGrid) LocateTriangle(x, z float64) (int, bool) {
	if sg == nil {
		return -1, false
	}

	// Find the cell
	cx := int((x - sg.minX) / sg.cellSize)
	cz := int((z - sg.minZ) / sg.cellSize)

	// Out of bounds check
	if cx < 0 || cx >= sg.gridW || cz < 0 || cz >= sg.gridH {
		return -1, false
	}

	cellIdx := cz*sg.gridW + cx
	candidates := sg.cells[cellIdx]

	p := Vec2{x, z}

	// Test each candidate triangle
	for _, ti := range candidates {
		if sg.pointInTriangle(ti, p) {
			return ti, true
		}
	}

	return -1, false
}

// pointInTriangle tests if point p is inside triangle ti using barycentric coordinates.
func (sg *SpatialGrid) pointInTriangle(ti int, p Vec2) bool {
	t := sg.mesh.Tris[ti]
	a := sg.mesh.Sites[t.A].Pos
	b := sg.mesh.Sites[t.B].Pos
	c := sg.mesh.Sites[t.C].Pos

	// Barycentric test
	v0 := b.Sub(a)
	v1 := c.Sub(a)
	v2 := p.Sub(a)

	denom := cross2(v0, v1)
	if math.Abs(denom) < 1e-12 {
		return false // Degenerate triangle
	}

	wb := cross2(v2, v1) / denom
	wc := cross2(v0, v2) / denom
	wa := 1.0 - wb - wc

	// Point is inside if all weights are non-negative (with small epsilon for edge cases)
	const eps = -1e-9
	return wa >= eps && wb >= eps && wc >= eps
}

// SampleHeight returns the interpolated height at position (x, z).
// Returns (height, ok) where ok=false if the point is outside all triangles.
func (sg *SpatialGrid) SampleHeight(x, z float64) (float64, bool) {
	ti, ok := sg.LocateTriangle(x, z)
	if !ok {
		return 0, false
	}

	p := Vec2{x, z}
	t := sg.mesh.Tris[ti]
	a := sg.mesh.Sites[t.A].Pos
	b := sg.mesh.Sites[t.B].Pos
	c := sg.mesh.Sites[t.C].Pos

	// Compute barycentric weights
	v0 := b.Sub(a)
	v1 := c.Sub(a)
	v2 := p.Sub(a)

	denom := cross2(v0, v1)
	if math.Abs(denom) < 1e-12 {
		return 0, false
	}

	wb := cross2(v2, v1) / denom
	wc := cross2(v0, v2) / denom
	wa := 1.0 - wb - wc

	// Interpolate height
	ha := sg.heights[t.A]
	hb := sg.heights[t.B]
	hc := sg.heights[t.C]

	return wa*ha + wb*hb + wc*hc, true
}

// GetTriangleNormal returns the face normal for the triangle at position (x, z).
// Returns (normal, ok) where ok=false if the point is outside all triangles.
func (sg *SpatialGrid) GetTriangleNormal(x, z float64, normals []Vec3) (Vec3, bool) {
	ti, ok := sg.LocateTriangle(x, z)
	if !ok || ti >= len(normals) {
		return Vec3{}, false
	}
	return normals[ti], true
}

// TrianglesInBounds returns all triangle indices whose bounding boxes overlap the given region.
func (sg *SpatialGrid) TrianglesInBounds(minX, minZ, maxX, maxZ float64) []int {
	if sg == nil {
		return nil
	}

	// Convert to grid cells
	cellMinX := int((minX - sg.minX) / sg.cellSize)
	cellMaxX := int((maxX - sg.minX) / sg.cellSize)
	cellMinZ := int((minZ - sg.minZ) / sg.cellSize)
	cellMaxZ := int((maxZ - sg.minZ) / sg.cellSize)

	// Clamp
	if cellMinX < 0 {
		cellMinX = 0
	}
	if cellMaxX >= sg.gridW {
		cellMaxX = sg.gridW - 1
	}
	if cellMinZ < 0 {
		cellMinZ = 0
	}
	if cellMaxZ >= sg.gridH {
		cellMaxZ = sg.gridH - 1
	}

	// Collect unique triangles
	seen := make(map[int]struct{})
	result := make([]int, 0)

	for cz := cellMinZ; cz <= cellMaxZ; cz++ {
		for cx := cellMinX; cx <= cellMaxX; cx++ {
			cellIdx := cz*sg.gridW + cx
			for _, ti := range sg.cells[cellIdx] {
				if _, exists := seen[ti]; !exists {
					seen[ti] = struct{}{}
					result = append(result, ti)
				}
			}
		}
	}

	return result
}

// Raycast performs a vertical ray test at position (x, z) and returns the height.
// This is an alias for SampleHeight, provided for semantic clarity in collision code.
func (sg *SpatialGrid) Raycast(x, z float64) (float64, bool) {
	return sg.SampleHeight(x, z)
}
