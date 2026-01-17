# Terrain Generation Study Guide

This document provides a comprehensive overview of the procedural terrain generation system implemented in the `duals` package. The system generates infinite, seamless terrain using Delaunay triangulation over blue noise distributed points.

## Table of Contents

1. [System Overview](#system-overview)
2. [Blue Noise Point Distribution](#blue-noise-point-distribution)
3. [Delaunay Triangulation](#delaunay-triangulation)
4. [The Walking Optimization](#the-walking-optimization)
5. [Chunk System and Streaming](#chunk-system-and-streaming)
6. [Half-Edge Mesh Topology](#half-edge-mesh-topology)
7. [Height Sampling and Interpolation](#height-sampling-and-interpolation)
8. [Spatial Indexing](#spatial-indexing)
9. [Voronoi Duality](#voronoi-duality)
10. [Performance Characteristics](#performance-characteristics)

---

## System Overview

The terrain generation pipeline follows these stages:

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   Blue Noise    │────▶│    Delaunay     │────▶│   Half-Edge     │
│   Generation    │     │  Triangulation  │     │   Mesh Build    │
└─────────────────┘     └─────────────────┘     └─────────────────┘
                                                        │
┌─────────────────┐     ┌─────────────────┐             ▼
│  Height Query   │◀────│  Spatial Index  │◀────────────┘
│  (Runtime)      │     │    (Grid)       │
└─────────────────┘     └─────────────────┘
```

**Key Design Goals:**

- **Deterministic**: Same seed produces identical terrain
- **Streamable**: Chunks can be generated independently and on-demand
- **Seamless**: No visible seams between adjacent chunks
- **Efficient**: O(n√n) triangulation, O(1) height queries

---

## Blue Noise Point Distribution

### What is Blue Noise?

Blue noise is a sampling pattern where points are randomly distributed but maintain a **minimum distance** from each other. Unlike pure random (white noise) which creates clumps and gaps, blue noise produces aesthetically pleasing, evenly-spaced, non-repeating distributions.

```
White Noise (clumpy):        Blue Noise:
  •    ••                      •    •   •
    •••    •                     •    •
  •      ••                    •    •   •
      •                          •    •
```

### Bridson's Algorithm

We use Bridson's algorithm for Poisson disk sampling, which generates blue noise in O(n) time:

```go
// From bluenoise.go
func GenerateBlueNoise(rng *rand.Rand, minX, minZ, maxX, maxZ float64, cfg BlueNoiseConfig) []Vec2
```

**Algorithm Steps:**

1. Create a spatial grid with cell size `r/√2` (guarantees at most one point per cell)
2. Start with a random seed point
3. For each "active" point:
   - Try up to `MaxTries` random candidates in the annulus `[r, 2r]`
   - If a valid candidate (far enough from all neighbors) is found, add it
   - If no valid candidate found, deactivate the point
4. Repeat until no active points remain

**Why Blue Noise for Terrain?**

- Avoids axis-aligned artifacts that regular grids produce
- Naturally adapts to any region shape
- The minimum distance constraint (`MinPointDist`) controls mesh density
- Produces well-shaped triangles (no slivers)

---

## Delaunay Triangulation

### What is Delaunay Triangulation?

Given a set of points, the Delaunay triangulation connects them into triangles such that:

- No point lies inside the circumcircle of any triangle
- Maximizes the minimum angle (produces "fat" triangles, avoiding slivers)

```
Circumcircle Property:
       ○───────○
      ╱ ╲     ╱
     ╱   ╲   ╱
    ○─────○─────○
         │
    No point inside
    any circumcircle
```

### Bowyer-Watson Algorithm

We implement the incremental Bowyer-Watson algorithm:

```go
// From chunk.go
func Triangulate(points []Vec2) []Triangle
```

**Basic Steps (Naive O(n²)):**

1. Create a super-triangle encompassing all points
2. For each point p:
   - Find all triangles whose circumcircle contains p ("bad triangles")
   - Remove bad triangles, leaving a star-shaped polygon hole
   - Create new triangles connecting p to each edge of the polygon
3. Remove triangles connected to the super-triangle

**The Naive Bottleneck:**

```go
// This was O(n) per point insertion → O(n²) total
for ti, t := range triangles {
    if inCircumcircle(allPts[t.a], allPts[t.b], allPts[t.c], p) {
        badTris = append(badTris, ti)
    }
}
```

---

## The Walking Optimization

### Problem: O(n²) Triangle Search

The naive algorithm checks every triangle for each point insertion. For 20,000 points, this means ~400 million circumcircle tests.

### Solution: Jump-and-Walk with Adjacency

We maintain triangle adjacency information and use it to "walk" toward the query point:

```go
// From chunk.go
type adjTri struct {
    a, b, c    int  // vertex indices (CCW)
    n0, n1, n2 int  // neighbor opposite to vertex a, b, c
    alive      bool // false = deleted (lazy deletion)
}
```

**Key Components:**

#### 1. Walking Point Location

```go
func walkToPoint(triangles []adjTri, pts []Vec2, startTri int, p Vec2) int
```

Starting from a known triangle, use orientation tests to determine which edge to cross:

```
Current Triangle:          Point p is "outside" edge BC
      A                    (negative orientation)
     ╱ ╲
    ╱   ╲    • p           Cross to neighbor n0
   ╱     ╲                 (opposite vertex A)
  B───────C
     n0
```

**Expected Steps:** O(√n) for spatially random points

#### 2. Flood-Fill Bad Triangles

Once we find one bad triangle, we flood through adjacency to find all connected bad triangles:

```go
func floodFillBadTris(triangles []adjTri, pts []Vec2, startTri int, p Vec2) []int
```

This exploits the fact that bad triangles form a connected region (star-shaped hole).

#### 3. Adjacency Maintenance

When inserting new triangles, we update adjacency in O(1) per triangle using an edge→triangle map:

```go
edgeToTri := make(map[edgeKey]triEdgeRef)
```

### Complexity Improvement

| Operation          | Naive     | Optimized                |
| ------------------ | --------- | ------------------------ |
| Point location     | O(n)      | O(√n) expected           |
| Find bad triangles | O(n)      | O(k) where k = bad count |
| Adjacency updates  | N/A       | O(1) per triangle        |
| **Total**          | **O(n²)** | **O(n√n)**               |

**Real-World Results (n=242 points):**

```
avgSteps=7.66, maxSteps=38, fallbacks=0, avgBadTris=3.98
```

For 20,000 points: ~140x speedup over naive approach.

---

## Chunk System and Streaming

### The Chunking Problem

For infinite terrain, we can't generate all points at once. We divide the world into chunks that can be generated independently.

```go
type ChunkCoord struct {
    X, Z int
}

type ChunkConfig struct {
    ChunkSize    float64 // World units per chunk (e.g., 256.0)
    MinPointDist float64 // Blue noise minimum distance
    HaloWidth    float64 // Overlap region for seamless boundaries
    WorldSeed    int64   // Global deterministic seed
}
```

### The Boundary Problem

If chunks are generated independently, triangles at chunk boundaries won't match:

```
Chunk A          │         Chunk B
                 │
    •──────•     │     •──────•
     ╲    ╱      │      ╲    ╱
      ╲  ╱       │       ╲  ╱
       ╲╱        │        ╲╱
        • ???    │    ??? •
                 │
          GAP/MISMATCH
```

### Solution: Halo Regions

Each chunk includes a "halo" of points from neighboring chunks:

```
┌─────────────────────────────────┐
│           HALO REGION           │
│   ┌───────────────────────┐     │
│   │                       │     │
│   │     CORE REGION       │     │
│   │   (chunk owns these)  │     │
│   │                       │     │
│   └───────────────────────┘     │
│                                 │
└─────────────────────────────────┘
```

**Key Insight:** Triangles are only rendered if they touch at least one core vertex, but halo vertices ensure consistent connectivity at boundaries.

### Deterministic Seeding

Each chunk's points are generated with a deterministic seed based on world seed and chunk coordinates:

```go
func chunkSeed(worldSeed int64, coord ChunkCoord) int64 {
    // FNV-1a hash of (worldSeed, coord.X, coord.Z)
    // Same inputs always produce same output
}
```

This ensures:

- Chunk (0,0) always has the same points regardless of generation order
- Neighbor chunks can regenerate each other's points identically

### Point Cache Optimization

Without caching, generating chunk A requires regenerating points for all 8 neighbors (to get halo data). With caching:

```go
type ChunkManager struct {
    cache       map[ChunkCoord]*TerrainChunk
    pointsCache map[ChunkCoord][]Vec2  // Lightweight point-only cache
}
```

**Cache Priority:**

1. Check if full chunk exists → use `Mesh.Sites`
2. Check if points cached → use cached points
3. Generate blue noise → cache immediately

**Result:** Each chunk's points generated at most once, regardless of access order.

---

## Half-Edge Mesh Topology

### What is a Half-Edge?

A half-edge is a directed edge that stores connectivity information:

```go
type HalfEdge struct {
    Origin   int // Site index where this edge starts
    Tri      int // Triangle this edge belongs to
    Next     int // Next edge in the triangle (CCW)
    Twin     int // Opposite direction edge (-1 if boundary)
    Prev     int // Previous edge in the triangle
    EdgeDest int // Destination site (cached for convenience)
}
```

**Visual Representation:**

```
        A
       ╱│╲
      ╱ │ ╲
     ╱  │  ╲
    ╱   │   ╲
   B────┼────C
        │
   Edge A→B has Twin B→A
   Edge A→B's Next is B→C
```

### Why Half-Edges?

Half-edge structure enables O(1) operations:

- **Find neighbors:** Follow `Twin` pointers
- **Walk around a vertex:** `Twin.Next` cycle
- **Find adjacent triangles:** `Twin.Tri`

### Building the Mesh

```go
func BuildHalfEdgeMesh(sites []Site, tris []Triangle) (*DelaunayMesh, error)
```

1. Create 3 half-edges per triangle (A→B, B→C, C→A)
2. Build edge map: (origin, dest) → half-edge index
3. Link twins by looking up reverse edge
4. Compute circumcenters for Voronoi dual

---

## Height Sampling and Interpolation

### Barycentric Coordinates

Any point inside a triangle can be expressed as a weighted sum of the vertices:

```
p = wa·A + wb·B + wc·C
where wa + wb + wc = 1
```

**Geometric Interpretation:** The weights represent the ratio of sub-triangle areas.

```go
func (m *DelaunayMesh) Barycentric(triID int, p Vec2) (wa, wb, wc float64, ok bool)
```

### Height Interpolation

Given heights at each vertex, interpolate smoothly across the triangle:

```go
height = wa·heightA + wb·heightB + wc·heightC
```

This produces a piecewise-linear surface (each triangle is a flat plane).

### Gradient and Normals

The gradient is constant across each triangle (it's a plane!):

```go
func (m *DelaunayMesh) TriangleGradient(triID int, values []float64) (dhdx, dhdz float64, ok bool)
```

The face normal is derived from the gradient:

```go
normal = normalize(-dhdx, 1, -dhdz)
```

---

## Spatial Indexing

### The Problem

Given a query point (x, z), which triangle contains it? Naive search is O(n).

### Solution: Uniform Grid

```go
type SpatialGrid struct {
    cells    [][]int  // Each cell stores overlapping triangle indices
    cellSize float64  // Should match ~MinPointDist
}
```

**Construction:**

1. For each triangle, compute bounding box
2. Add triangle index to all cells it overlaps

**Query:**

1. Map query point to cell: O(1)
2. Test candidates in that cell: O(small constant)

```go
func (sg *SpatialGrid) LocateTriangle(x, z float64) (int, bool)
func (sg *SpatialGrid) SampleHeight(x, z float64) (float64, bool)
```

**Result:** O(1) amortized height queries.

---

## Voronoi Duality

### Delaunay ↔ Voronoi

The Delaunay triangulation and Voronoi diagram are **duals**:

- Delaunay vertices → Voronoi cells
- Delaunay triangles → Voronoi vertices (circumcenters)
- Delaunay edges → Voronoi edges

```
Delaunay:                 Voronoi:
    A─────B                   │
   ╱ ╲   ╱ ╲              ────┼────
  ╱   ╲ ╱   ╲                 │
 C─────D─────E            ────┼────
                              │
```

### Computing Voronoi Cells

```go
func (m *DelaunayMesh) VoronoiCellForSite(site int) VoronoiCell
```

Walk around the site using half-edge topology, collecting circumcenters:

```go
type VoronoiCell struct {
    Site     int     // Which site this cell belongs to
    Vertices []Vec2  // Circumcenters of adjacent triangles
    Closed   bool    // False if cell is unbounded (boundary site)
}
```

**Applications:**

- Natural region boundaries for biomes
- Influence zones for game mechanics
- Procedural region generation

---

## Performance Characteristics

### Generation Complexity

| Component              | Complexity | Notes                |
| ---------------------- | ---------- | -------------------- |
| Blue noise generation  | O(n)       | Bridson's algorithm  |
| Triangulation          | O(n√n)     | Walking optimization |
| Half-edge construction | O(n)       | Linear in triangles  |
| Spatial index build    | O(n)       | Linear in triangles  |

### Runtime Query Complexity

| Operation         | Complexity                   |
| ----------------- | ---------------------------- |
| Height at point   | O(1) amortized               |
| Triangle location | O(1) amortized               |
| Normal at point   | O(1) amortized               |
| Voronoi cell      | O(k) where k = cell vertices |

### Memory Usage

For n points:

- ~2n triangles (Euler's formula)
- ~6n half-edges
- Point cache: O(chunks × points/chunk)
- Spatial grid: O(n) triangle references

### Measured Results

For n=242 points per chunk:

```
Walking stats:
  avgSteps = 7.66 (expected: √242 ≈ 15.6)
  maxSteps = 38
  fallbacks = 0 (walking always succeeds)
  avgBadTris = 3.98 (typical for Delaunay)
```

---

## Summary

This terrain generation system combines several elegant algorithms:

1. **Blue Noise** provides aesthetically pleasing, evenly-distributed vertices
2. **Delaunay Triangulation** connects them optimally (maximizing minimum angles)
3. **Walking Optimization** achieves O(n√n) instead of O(n²)
4. **Chunk System** enables infinite, streamable terrain
5. **Point Caching** avoids redundant generation
6. **Half-Edge Topology** enables efficient mesh traversal
7. **Spatial Indexing** provides O(1) runtime queries

The result is a production-ready terrain system suitable for games and simulations, with clear pathways to C++ optimization for maximum performance.
