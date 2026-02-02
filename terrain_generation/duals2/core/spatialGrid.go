package core

// SpatialGrid provides O(1) point location for height queries on a Delaunay mesh.
type SpatialGrid struct {
	Mesh     *DelaunayMesh
	Heights  []float64
	CellSize float64
	MinX, MinZ float64
	MaxX, MaxZ float64
	GridW, GridH int

	// Each cell stores a list of triangle indices that overlap it
	Cells [][]int
}

// LocateTriangle finds the triangle containing point (x, z).
// Returns (triangleIndex, ok) where ok=false if the point is outside all triangles.
func (sg *SpatialGrid) LocateTriangle(x, z float64) (int, bool) {
	if sg == nil {
		return -1, false
	}

	// Find the cell
	cx := int((x - sg.MinX) / sg.CellSize) // Note: int() will floor
	cz := int((z - sg.MinZ) / sg.CellSize)

	// Out of bounds check
	if cx < 0 || cx >= sg.GridW || cz < 0 || cz >= sg.GridH {
		return -1, false
	}

	// row * width + offset (everyone's favorite formula)
	cellIdx := cz*sg.GridW + cx
	candidates := sg.Cells[cellIdx]

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

	wa, wb, wc, ok := sg.Mesh.Barycentric(ti, p)
	if !ok {
		return false
	}

	// Point is inside if all weights are non-negative (with small epsilon for edge cases)
	const eps = -1e-9
	return wa >= eps && wb >= eps && wc >= eps
}

// SEE ALSO: Mesh.SampleScalar() for when you already know the triangle index
// SampleHeight returns the interpolated height at position (x, z).
// Returns (height, ok) where ok=false if the point is outside all triangles.
func (sg *SpatialGrid) SampleHeight(x, z float64) (float64, bool) {
	ti, ok := sg.LocateTriangle(x, z)
	if !ok {
		return 0, false
	}

	p := Vec2{x, z}
	t := sg.Mesh.Tris[ti]
	wa, wb, wc, ok := sg.Mesh.Barycentric(ti, p)
	if !ok {
		return 0, false
	}

	// Discreet samples
	ha := sg.Heights[t.A]
	hb := sg.Heights[t.B]
	hc := sg.Heights[t.C]

	// Interpolated heights
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
// Providing min and max values allow us to cover more than 1x1 cells, we can select entire rectangular regions.
func (sg *SpatialGrid) TrianglesInBounds(minX, minZ, maxX, maxZ float64) []int {
	if sg == nil {
		return nil
	}

	// Convert to grid cells
	cellMinX := int((minX - sg.MinX) / sg.CellSize)
	cellMaxX := int((maxX - sg.MinX) / sg.CellSize)
	cellMinZ := int((minZ - sg.MinZ) / sg.CellSize)
	cellMaxZ := int((maxZ - sg.MinZ) / sg.CellSize)

	// Clamp
	if cellMinX < 0 {
		cellMinX = 0
	}
	if cellMaxX >= sg.GridW {
		cellMaxX = sg.GridW - 1
	}
	if cellMinZ < 0 {
		cellMinZ = 0
	}
	if cellMaxZ >= sg.GridH {
		cellMaxZ = sg.GridH - 1
	}

	// Collect unique triangles
	seen := make(map[int]struct{})
	result := make([]int, 0)

	for cz := cellMinZ; cz <= cellMaxZ; cz++ {
		for cx := cellMinX; cx <= cellMaxX; cx++ {
			cellIdx := cz*sg.GridW + cx
			for _, ti := range sg.Cells[cellIdx] {
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
