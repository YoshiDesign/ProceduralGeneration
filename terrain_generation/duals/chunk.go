package duals

import (
	"hash/fnv"
	"math/rand"
	"sync"
)

// ChunkCoord identifies a chunk by integer grid coordinates.
type ChunkCoord struct {
	X, Z int
}

// ChunkManager manages chunk generation and caching for optimized neighbor lookups.
type ChunkManager struct {
	mu         sync.RWMutex
	cache      map[ChunkCoord]*TerrainChunk
	cfg        ChunkConfig
	heightFunc func(x, z float64) float64
}

// NewChunkManager creates a new chunk manager with the given configuration.
func NewChunkManager(cfg ChunkConfig, heightFunc func(x, z float64) float64) *ChunkManager {
	return &ChunkManager{
		cache:      make(map[ChunkCoord]*TerrainChunk),
		cfg:        cfg,
		heightFunc: heightFunc,
	}
}

// GetOrGenerate returns a cached chunk or generates and caches a new one.
// This is the primary API for chunk access with caching optimization.
func (cm *ChunkManager) GetOrGenerate(coord ChunkCoord) (*TerrainChunk, error) {
	// Fast path: check if already cached
	cm.mu.RLock()
	if chunk, ok := cm.cache[coord]; ok {
		cm.mu.RUnlock()
		return chunk, nil
	}
	cm.mu.RUnlock()

	// Slow path: generate with cache-aware point collection
	cm.mu.Lock()
	defer cm.mu.Unlock()

	// Double-check after acquiring write lock
	if chunk, ok := cm.cache[coord]; ok {
		return chunk, nil
	}

	chunk, err := generateChunkWithCache(coord, cm.cfg, cm.heightFunc, cm.cache)
	if err != nil {
		return nil, err
	}

	cm.cache[coord] = chunk
	return chunk, nil
}

// Get returns a cached chunk without generating. Returns nil if not cached.
func (cm *ChunkManager) Get(coord ChunkCoord) *TerrainChunk {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return cm.cache[coord]
}

// Evict removes a chunk from the cache, freeing memory.
func (cm *ChunkManager) Evict(coord ChunkCoord) {
	cm.mu.Lock()
	defer cm.mu.Unlock()
	delete(cm.cache, coord)
}

// EvictOutsideRadius removes all chunks outside the given radius from center.
func (cm *ChunkManager) EvictOutsideRadius(center ChunkCoord, radius int) int {
	cm.mu.Lock()
	defer cm.mu.Unlock()

	evicted := 0
	for coord := range cm.cache {
		dx := coord.X - center.X
		dz := coord.Z - center.Z
		if dx*dx+dz*dz > radius*radius {
			delete(cm.cache, coord)
			evicted++
		}
	}
	return evicted
}

// CacheSize returns the number of cached chunks.
func (cm *ChunkManager) CacheSize() int {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	return len(cm.cache)
}

// CachedCoords returns all cached chunk coordinates.
func (cm *ChunkManager) CachedCoords() []ChunkCoord {
	cm.mu.RLock()
	defer cm.mu.RUnlock()
	coords := make([]ChunkCoord, 0, len(cm.cache))
	for coord := range cm.cache {
		coords = append(coords, coord)
	}
	return coords
}

// ChunkConfig holds parameters for terrain chunk generation.
type ChunkConfig struct {
	ChunkSize   float64 // World units per chunk side (e.g., 256.0)
	MinPointDist float64 // Minimum distance between blue noise points
	HaloWidth   float64 // Boundary overlap region width (typically = MinPointDist)
	WorldSeed   int64   // Global world seed
}

// DefaultChunkConfig returns sensible defaults for a terrain chunk.
func DefaultChunkConfig() ChunkConfig {
	return ChunkConfig{
		ChunkSize:    256.0,
		MinPointDist: 8.0,
		HaloWidth:    8.0,
		WorldSeed:    42,
	}
}

// TerrainChunk represents a generated terrain chunk with mesh data.
type TerrainChunk struct {
	Coord ChunkCoord
	Cfg   ChunkConfig

	// Core bounds (what this chunk "owns")
	MinX, MinZ float64
	MaxX, MaxZ float64

	// The Delaunay mesh (includes halo points for boundary continuity)
	Mesh *DelaunayMesh

	// Height values per site (parallel to Mesh.Sites)
	Heights []float64

	// Face normals per triangle (parallel to Mesh.Tris)
	FaceNormals []Vec3

	// Spatial index for fast point location
	Spatial *SpatialGrid

	// Which sites are in the core region (not halo)
	CoreSiteIndices []int
}

// chunkSeed computes a deterministic seed for a chunk based on world seed and coordinates.
func chunkSeed(worldSeed int64, coord ChunkCoord) int64 {
	h := fnv.New64a()
	// Write world seed
	buf := make([]byte, 8)
	buf[0] = byte(worldSeed)
	buf[1] = byte(worldSeed >> 8)
	buf[2] = byte(worldSeed >> 16)
	buf[3] = byte(worldSeed >> 24)
	buf[4] = byte(worldSeed >> 32)
	buf[5] = byte(worldSeed >> 40)
	buf[6] = byte(worldSeed >> 48)
	buf[7] = byte(worldSeed >> 56)
	h.Write(buf)

	// Write chunk X
	buf[0] = byte(coord.X)
	buf[1] = byte(coord.X >> 8)
	buf[2] = byte(coord.X >> 16)
	buf[3] = byte(coord.X >> 24)
	buf[4] = 0
	buf[5] = 0
	buf[6] = 0
	buf[7] = 0
	h.Write(buf)

	// Write chunk Z
	buf[0] = byte(coord.Z)
	buf[1] = byte(coord.Z >> 8)
	buf[2] = byte(coord.Z >> 16)
	buf[3] = byte(coord.Z >> 24)
	h.Write(buf)

	return int64(h.Sum64())
}

// boundarySeed computes a deterministic seed for the shared boundary between two chunks.
// It uses the minimum of the two chunk coords to ensure both chunks get the same seed.
func boundarySeed(worldSeed int64, coord1, coord2 ChunkCoord) int64 {
	// Use lexicographically smaller coord first
	var first, second ChunkCoord
	if coord1.X < coord2.X || (coord1.X == coord2.X && coord1.Z < coord2.Z) {
		first, second = coord1, coord2
	} else {
		first, second = coord2, coord1
	}

	h := fnv.New64a()
	buf := make([]byte, 8)

	// World seed
	buf[0] = byte(worldSeed)
	buf[1] = byte(worldSeed >> 8)
	buf[2] = byte(worldSeed >> 16)
	buf[3] = byte(worldSeed >> 24)
	buf[4] = byte(worldSeed >> 32)
	buf[5] = byte(worldSeed >> 40)
	buf[6] = byte(worldSeed >> 48)
	buf[7] = byte(worldSeed >> 56)
	h.Write(buf)

	// First chunk
	buf[0] = byte(first.X)
	buf[1] = byte(first.X >> 8)
	buf[2] = byte(first.Z)
	buf[3] = byte(first.Z >> 8)
	// Second chunk
	buf[4] = byte(second.X)
	buf[5] = byte(second.X >> 8)
	buf[6] = byte(second.Z)
	buf[7] = byte(second.Z >> 8)
	h.Write(buf)

	return int64(h.Sum64())
}

// GenerateChunk creates a terrain chunk with Delaunay mesh, heights, and spatial index.
// The heightFunc provides elevation for each site position.
// This version does not use caching - for cached generation, use ChunkManager.
func GenerateChunk(coord ChunkCoord, cfg ChunkConfig, heightFunc func(x, z float64) float64) (*TerrainChunk, error) {
	return generateChunkWithCache(coord, cfg, heightFunc, nil)
}

// generateChunkWithCache creates a terrain chunk, optionally using a cache for neighbor lookups.
// When cache is provided, it extracts halo points from cached neighbors instead of regenerating them.
func generateChunkWithCache(coord ChunkCoord, cfg ChunkConfig, heightFunc func(x, z float64) float64, cache map[ChunkCoord]*TerrainChunk) (*TerrainChunk, error) {
	chunk := &TerrainChunk{
		Coord: coord,
		Cfg:   cfg,
		MinX:  float64(coord.X) * cfg.ChunkSize,
		MinZ:  float64(coord.Z) * cfg.ChunkSize,
		MaxX:  float64(coord.X+1) * cfg.ChunkSize,
		MaxZ:  float64(coord.Z+1) * cfg.ChunkSize,
	}

	// Generate points: core region + halo regions for each neighbor
	allPoints := generateChunkPointsWithCache(coord, cfg, cache)

	// Track which points are in the core region
	coreIndices := make([]int, 0, len(allPoints))
	for i, p := range allPoints {
		if p.X >= chunk.MinX && p.X < chunk.MaxX && p.Y >= chunk.MinZ && p.Y < chunk.MaxZ {
			coreIndices = append(coreIndices, i)
		}
	}
	chunk.CoreSiteIndices = coreIndices

	// Build sites with heights
	sites := make([]Site, len(allPoints))
	heights := make([]float64, len(allPoints))
	for i, p := range allPoints {
		h := 0.0
		if heightFunc != nil {
			h = heightFunc(p.X, p.Y) // p.Y is Z in 2D
		}
		sites[i] = Site{Pos: p, Height: h}
		heights[i] = h
	}
	chunk.Heights = heights

	// Build Delaunay triangulation
	tris := Triangulate(allPoints)
	mesh, err := BuildHalfEdgeMesh(sites, tris)
	if err != nil {
		return nil, err
	}
	chunk.Mesh = mesh

	// Compute face normals
	chunk.FaceNormals = computeAllFaceNormals(mesh, heights)

	// Build spatial index
	// Cell size roughly equal to minimum point distance for good performance
	chunk.Spatial = BuildSpatialIndex(mesh, heights, cfg.MinPointDist,
		chunk.MinX-cfg.HaloWidth, chunk.MinZ-cfg.HaloWidth,
		chunk.MaxX+cfg.HaloWidth, chunk.MaxZ+cfg.HaloWidth)

	return chunk, nil
}

// generateChunkPoints generates blue noise points for a chunk including halo regions.
// Points in boundary regions are generated with shared seeds to ensure continuity.
// Deprecated: Use generateChunkPointsWithCache for better performance.
func generateChunkPoints(coord ChunkCoord, cfg ChunkConfig) []Vec2 {
	return generateChunkPointsWithCache(coord, cfg, nil)
}

// generateChunkPointsWithCache generates blue noise points for a chunk including halo regions.
// When cache is provided, it extracts halo points from cached neighbors (O(n) filter)
// instead of regenerating them (O(n) blue noise generation per neighbor).
func generateChunkPointsWithCache(coord ChunkCoord, cfg ChunkConfig, cache map[ChunkCoord]*TerrainChunk) []Vec2 {
	minX := float64(coord.X) * cfg.ChunkSize
	minZ := float64(coord.Z) * cfg.ChunkSize
	maxX := float64(coord.X+1) * cfg.ChunkSize
	maxZ := float64(coord.Z+1) * cfg.ChunkSize

	halo := cfg.HaloWidth
	blueCfg := DefaultBlueNoiseConfig(cfg.MinPointDist)

	// We'll collect points from multiple regions
	allPoints := make([]Vec2, 0, 1024)
	seen := make(map[uint64]struct{}, 1024)

	// Hash a point to detect duplicates (within tolerance)
	hashPoint := func(p Vec2) uint64 {
		// Quantize to half the min distance for dedup
		scale := 2.0 / cfg.MinPointDist
		qx := int64(p.X * scale)
		qz := int64(p.Y * scale)
		return uint64(qx)<<32 | uint64(qz)&0xFFFFFFFF
	}

	addPoint := func(p Vec2) bool {
		h := hashPoint(p)
		if _, exists := seen[h]; !exists {
			seen[h] = struct{}{}
			allPoints = append(allPoints, p)
			return true
		}
		return false
	}

	addPoints := func(pts []Vec2) {
		for _, p := range pts {
			addPoint(p)
		}
	}

	// 1. Generate core points for this chunk
	coreSeed := chunkSeed(cfg.WorldSeed, coord)
	corePoints := GenerateBlueNoiseSeeded(coreSeed, minX, minZ, maxX, maxZ, blueCfg)
	addPoints(corePoints)

	// 2. Get halo points from each neighboring chunk
	neighbors := []ChunkCoord{
		{coord.X - 1, coord.Z - 1}, {coord.X, coord.Z - 1}, {coord.X + 1, coord.Z - 1},
		{coord.X - 1, coord.Z}, {coord.X + 1, coord.Z},
		{coord.X - 1, coord.Z + 1}, {coord.X, coord.Z + 1}, {coord.X + 1, coord.Z + 1},
	}

	for _, neighbor := range neighbors {
		// The boundary region is where this chunk's halo overlaps the neighbor
		nMinX := float64(neighbor.X) * cfg.ChunkSize
		nMinZ := float64(neighbor.Z) * cfg.ChunkSize
		nMaxX := float64(neighbor.X+1) * cfg.ChunkSize
		nMaxZ := float64(neighbor.Z+1) * cfg.ChunkSize

		// Compute the overlap region between our extended bounds and neighbor's core
		overlapMinX := max(minX-halo, nMinX)
		overlapMinZ := max(minZ-halo, nMinZ)
		overlapMaxX := min(maxX+halo, nMaxX)
		overlapMaxZ := min(maxZ+halo, nMaxZ)

		if overlapMinX >= overlapMaxX || overlapMinZ >= overlapMaxZ {
			continue // No overlap
		}

		// Check if neighbor is cached - if so, extract points directly (fast path)
		if cache != nil {
			if cachedNeighbor, ok := cache[neighbor]; ok {
				// Extract halo points from cached neighbor's sites
				for _, site := range cachedNeighbor.Mesh.Sites {
					p := site.Pos
					// Check if point is in our halo region but not our core
					inHalo := (p.X >= minX-halo && p.X < maxX+halo && p.Y >= minZ-halo && p.Y < maxZ+halo)
					inCore := (p.X >= minX && p.X < maxX && p.Y >= minZ && p.Y < maxZ)
					if inHalo && !inCore {
						addPoint(p)
					}
				}
				continue // Skip regeneration for this neighbor
			}
		}

		// Fallback: regenerate neighbor's points (slow path)
		neighborSeed := chunkSeed(cfg.WorldSeed, neighbor)
		neighborPoints := GenerateBlueNoiseSeeded(neighborSeed, nMinX, nMinZ, nMaxX, nMaxZ, blueCfg)

		// Filter to points within our extended bounds but outside our core
		for _, p := range neighborPoints {
			inHalo := (p.X >= minX-halo && p.X < maxX+halo && p.Y >= minZ-halo && p.Y < maxZ+halo)
			inCore := (p.X >= minX && p.X < maxX && p.Y >= minZ && p.Y < maxZ)
			if inHalo && !inCore {
				addPoint(p)
			}
		}
	}

	return allPoints
}

// Triangulate performs Delaunay triangulation on a set of 2D points.
// This is a simple implementation using the Bowyer-Watson algorithm.
func Triangulate(points []Vec2) []Triangle {
	if len(points) < 3 {
		return nil
	}

	// Find bounding box
	minX, maxX := points[0].X, points[0].X
	minZ, maxZ := points[0].Y, points[0].Y
	for _, p := range points {
		if p.X < minX {
			minX = p.X
		}
		if p.X > maxX {
			maxX = p.X
		}
		if p.Y < minZ {
			minZ = p.Y
		}
		if p.Y > maxZ {
			maxZ = p.Y
		}
	}

	// Create super-triangle that encompasses all points
	dx := maxX - minX
	dz := maxZ - minZ
	deltaMax := dx
	if dz > dx {
		deltaMax = dz
	}
	midX := (minX + maxX) / 2
	midZ := (minZ + maxZ) / 2

	// Super-triangle vertices (large enough to contain all points)
	superA := Vec2{midX - 20*deltaMax, midZ - deltaMax}
	superB := Vec2{midX, midZ + 20*deltaMax}
	superC := Vec2{midX + 20*deltaMax, midZ - deltaMax}

	// Working list of points including super-triangle
	allPts := make([]Vec2, len(points)+3)
	copy(allPts, points)
	allPts[len(points)] = superA
	allPts[len(points)+1] = superB
	allPts[len(points)+2] = superC

	superI := [3]int{len(points), len(points) + 1, len(points) + 2}

	// Triangle type for the algorithm
	type tri struct {
		a, b, c int
	}

	triangles := []tri{{superI[0], superI[1], superI[2]}}

	// Insert points one at a time
	for pi := 0; pi < len(points); pi++ {
		p := allPts[pi]

		// Find triangles whose circumcircle contains p
		badTris := make([]int, 0)
		for ti, t := range triangles {
			if inCircumcircle(allPts[t.a], allPts[t.b], allPts[t.c], p) {
				badTris = append(badTris, ti)
			}
		}

		// Find the boundary polygon of the bad triangles
		type edge struct{ a, b int }
		edgeCount := make(map[edge]int)
		for _, ti := range badTris {
			t := triangles[ti]
			edges := []edge{{t.a, t.b}, {t.b, t.c}, {t.c, t.a}}
			for _, e := range edges {
				// Normalize edge direction for counting
				if e.a > e.b {
					e.a, e.b = e.b, e.a
				}
				edgeCount[e]++
			}
		}

		// Polygon edges are those that appear exactly once
		polygon := make([]edge, 0)
		for _, ti := range badTris {
			t := triangles[ti]
			edges := []edge{{t.a, t.b}, {t.b, t.c}, {t.c, t.a}}
			for _, e := range edges {
				ne := e
				if ne.a > ne.b {
					ne.a, ne.b = ne.b, ne.a
				}
				if edgeCount[ne] == 1 {
					polygon = append(polygon, e)
				}
			}
		}

		// Remove bad triangles (in reverse order to preserve indices)
		for i := len(badTris) - 1; i >= 0; i-- {
			ti := badTris[i]
			triangles[ti] = triangles[len(triangles)-1]
			triangles = triangles[:len(triangles)-1]
		}

		// Create new triangles from polygon edges to the new point
		for _, e := range polygon {
			triangles = append(triangles, tri{e.a, e.b, pi})
		}
	}

	// Remove triangles that share vertices with the super-triangle
	result := make([]Triangle, 0, len(triangles))
	for _, t := range triangles {
		usesSuperVertex := false
		for _, si := range superI {
			if t.a == si || t.b == si || t.c == si {
				usesSuperVertex = true
				break
			}
		}
		if !usesSuperVertex {
			// Ensure CCW winding
			result = append(result, ensureCCW(allPts, t.a, t.b, t.c))
		}
	}

	return result
}

// inCircumcircle returns true if point p is inside the circumcircle of triangle (a, b, c).
func inCircumcircle(a, b, c, p Vec2) bool {
	// Using the determinant method
	ax, ay := a.X-p.X, a.Y-p.Y
	bx, by := b.X-p.X, b.Y-p.Y
	cx, cy := c.X-p.X, c.Y-p.Y

	det := (ax*ax+ay*ay)*(bx*cy-cx*by) -
		(bx*bx+by*by)*(ax*cy-cx*ay) +
		(cx*cx+cy*cy)*(ax*by-bx*ay)

	// For CCW triangles, det > 0 means p is inside
	// Handle both windings by checking the sign of the triangle area
	area := (b.X-a.X)*(c.Y-a.Y) - (c.X-a.X)*(b.Y-a.Y)
	if area < 0 {
		return det < 0
	}
	return det > 0
}

// ensureCCW returns a Triangle with vertices in counter-clockwise order.
func ensureCCW(pts []Vec2, a, b, c int) Triangle {
	// Cross product of (b-a) × (c-a)
	cross := (pts[b].X-pts[a].X)*(pts[c].Y-pts[a].Y) - (pts[c].X-pts[a].X)*(pts[b].Y-pts[a].Y)
	if cross < 0 {
		return Triangle{A: a, C: b, B: c} // Swap b and c
	}
	return Triangle{A: a, B: b, C: c}
}

// computeAllFaceNormals computes face normals for all triangles in the mesh.
func computeAllFaceNormals(mesh *DelaunayMesh, heights []float64) []Vec3 {
	normals := make([]Vec3, len(mesh.Tris))
	for i, t := range mesh.Tris {
		// Get 3D positions
		a := Vec3{mesh.Sites[t.A].Pos.X, heights[t.A], mesh.Sites[t.A].Pos.Y}
		b := Vec3{mesh.Sites[t.B].Pos.X, heights[t.B], mesh.Sites[t.B].Pos.Y}
		c := Vec3{mesh.Sites[t.C].Pos.X, heights[t.C], mesh.Sites[t.C].Pos.Y}

		// Two edges
		ab := b.Sub(a)
		ac := c.Sub(a)

		// Cross product gives normal (CCW winding means this points "up")
		n := ab.Cross(ac).Normalize()
		normals[i] = n
	}
	return normals
}

// SampleHeight returns the interpolated height at position (x, z).
// Returns (height, ok) where ok=false if the point is outside the mesh.
func (c *TerrainChunk) SampleHeight(x, z float64) (float64, bool) {
	if c.Spatial == nil {
		return 0, false
	}
	return c.Spatial.SampleHeight(x, z)
}

// Helper: seeded RNG for boundary regions
func boundaryRNG(worldSeed int64, coord1, coord2 ChunkCoord) *rand.Rand {
	seed := boundarySeed(worldSeed, coord1, coord2)
	return rand.New(rand.NewSource(seed))
}
