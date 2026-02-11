# C++ Vulkan Porting Guide: Terrain Generation System

This document outlines the strategy for porting the Go terrain generation prototype to a C++ Vulkan engine, with emphasis on collision detection, compute shaders for vertex normals, and cache-coherent SoA data layouts.

---

## Table of Contents

1. [Data Layout Philosophy (SoA)](#data-layout-philosophy-soa)
2. [SoA Structures Reference](#soa-structures-reference)
3. [Vertex Buffers](#vertex-buffers)
4. [Index Buffers](#index-buffers)
5. [Compute Shader: Vertex Normals](#compute-shader-vertex-normals)
6. [Collision Detection System](#collision-detection-system)
7. [Projectile Collision at Distance](#projectile-collision-at-distance)
8. [Batch Collision Queries](#batch-collision-queries)
9. [Cache Coherence Guidelines](#cache-coherence-guidelines)
10. [Runtime Terrain Modification](#runtime-terrain-modification)

---

## Data Layout Philosophy (SoA)

Structure of Arrays (SoA) is preferred over Array of Structures (AoS) for:

- **Cache efficiency**: Sequential memory access when iterating one attribute across many elements
- **SIMD friendliness**: Contiguous data enables vectorized operations
- **GPU compatibility**: Maps naturally to SSBOs and vertex attributes

**Exception**: Hybrid SoAoS (Structure of Arrays of Structures) is used when a small, fixed-size grouping improves access patterns (e.g., VertexAdjacency).

---

## SoA Structures Reference

### Sites (Vertices)

The fundamental terrain vertices. Parallel arrays indexed by `SiteIndex`.

```cpp
// Per-site data (parallel arrays, count = numSites)
struct SitesData {
    std::vector<float> posX;      // X position (world space)
    std::vector<float> posZ;      // Z position (world space)
    std::vector<float> height;    // Y position (mutable - can change at runtime)
    
    // Computed by compute shader (output)
    std::vector<float> normalX;   // Vertex normal X
    std::vector<float> normalY;   // Vertex normal Y
    std::vector<float> normalZ;   // Vertex normal Z
};
```

**Role**: Source geometry for vertex buffers, collision height queries, and normal computation.

---

### Triangles (Tris)

Triangle connectivity. Parallel arrays indexed by `TriIndex`.

```cpp
// Per-triangle data (parallel arrays, count = numTriangles)
struct TrianglesData {
    std::vector<uint32_t> siteA;  // Index into Sites for vertex A
    std::vector<uint32_t> siteB;  // Index into Sites for vertex B
    std::vector<uint32_t> siteC;  // Index into Sites for vertex C
    
    // Precomputed geometry (immutable after triangulation)
    std::vector<float> baryDenom; // 1.0 / cross(AB, AC) - barycentric denominator
    std::vector<float> abX, abY;  // Edge vector A->B (XZ plane)
    std::vector<float> acX, acY;  // Edge vector A->C (XZ plane)
    
    // Face normals (recomputed when heights change)
    std::vector<float> faceNormalX;
    std::vector<float> faceNormalY;
    std::vector<float> faceNormalZ;
};
```

**Role**: Index buffer source, collision triangle testing, face normal storage.

---

### HalfEdges

Topological connectivity for mesh traversal.

```cpp
// Per-half-edge data (parallel arrays, count = numHalfEdges = numTriangles * 3)
struct HalfEdgesData {
    std::vector<uint32_t> origin;   // Site index of edge origin
    std::vector<uint32_t> tri;      // Triangle index this edge belongs to
    std::vector<int32_t>  twin;     // Opposite half-edge index (-1 if boundary)
    std::vector<uint32_t> next;     // Next half-edge in triangle (CCW)
};
```

**Role**: Mesh traversal for adjacency queries, Voronoi cell computation, boundary detection.

---

### VertexAdjacency (SoAoS for Compute Shader)

Maps each vertex to its adjacent triangles. Used by compute shader to average face normals into vertex normals.

```cpp
// Hybrid SoAoS: small fixed-size structure stored in array
// Delaunay meshes typically have 5-7 triangles per vertex (average ~6)
constexpr uint32_t MAX_ADJACENT_TRIS = 12;  // Conservative upper bound

struct alignas(16) VertexAdjacency {
    uint32_t triangleIndices[MAX_ADJACENT_TRIS];  // Which triangles share this vertex
    uint32_t count;                                // Actual count (1-12 typically)
    uint32_t _padding[3];                          // Align to 64 bytes for GPU
};

// Packed into SSBO
std::vector<VertexAdjacency> adjacency;  // One per site, parallel to SitesData
```

**Why SoAoS here**: Each compute shader invocation needs ALL adjacent triangle indices for ONE vertex. Grouping them together avoids scattered reads.

**Triangle count per vertex**: For Delaunay triangulations from blue noise:
- Interior vertices: typically 5-7 adjacent triangles (hexagonal neighborhood)
- Boundary vertices: 3-5 adjacent triangles
- Maximum observed: rarely exceeds 10

---

### SpatialGrid (Collision Acceleration)

Grid-based spatial index for O(1) triangle lookup.

```cpp
struct SpatialGrid {
    // Grid metadata
    float cellSize;
    float minX, minZ;
    float maxX, maxZ;
    int32_t gridW, gridH;
    
    // Cell contents (flattened SoA)
    std::vector<uint32_t> cellOffsets;   // Where each cell's data starts (size = gridW * gridH + 1)
    std::vector<uint32_t> cellTriangles; // Triangle indices, all cells concatenated
    
    // Query: cell index = (z / cellSize) * gridW + (x / cellSize)
    // Triangles for cell i are: cellTriangles[cellOffsets[i] .. cellOffsets[i+1])
};
```

**Role**: Fast triangle lookup for collision queries. One per chunk.

---

### ChunkManager

World-level organization for infinite terrain.

```cpp
struct ChunkManager {
    std::unordered_map<ChunkCoord, TerrainChunk*, ChunkCoordHash> chunks;
    float chunkSize;  // World units per chunk side
    
    TerrainChunk* getChunk(float worldX, float worldZ) {
        ChunkCoord coord = {
            static_cast<int>(std::floor(worldX / chunkSize)),
            static_cast<int>(std::floor(worldZ / chunkSize))
        };
        auto it = chunks.find(coord);
        return (it != chunks.end()) ? it->second : nullptr;
    }
};
```

**Role**: O(1) chunk lookup for any world position.

---

## Vertex Buffers

Terrain vertex buffers are constructed from the SoA site data.

### Vertex Attributes

| Attribute | Type | Source | Notes |
|-----------|------|--------|-------|
| Position | `vec3` | `posX`, `height`, `posZ` | Interleaved at buffer creation |
| Normal | `vec3` | Compute shader output | Written by compute, read by vertex shader |
| TexCoord | `vec2` | Derived from `posX`, `posZ` | Optional, for texture mapping |

### Buffer Layout Options

**Option A: Separate Buffers (SoA on GPU)**
```cpp
// Binding 0: Positions
VkBuffer positionBuffer;  // vec3[] - can be immutable for XZ, mutable for Y

// Binding 1: Normals (compute shader output)
VkBuffer normalBuffer;    // vec3[] - written by compute shader

// Binding 2: TexCoords
VkBuffer texCoordBuffer;  // vec2[]
```

**Option B: Interleaved Buffer (AoS on GPU)**
```cpp
struct TerrainVertex {
    float posX, posY, posZ;
    float normalX, normalY, normalZ;
    float texU, texV;
};
VkBuffer vertexBuffer;  // TerrainVertex[]
```

**Recommendation**: Option A (separate buffers) if heights change frequently, since you only update the position buffer's Y components. Option B if terrain is mostly static after generation.

---

## Index Buffers

Triangle indices for indexed drawing.

```cpp
// Index buffer construction from TrianglesData
std::vector<uint32_t> indices;
indices.reserve(numTriangles * 3);

for (uint32_t i = 0; i < numTriangles; i++) {
    indices.push_back(tris.siteA[i]);
    indices.push_back(tris.siteB[i]);
    indices.push_back(tris.siteC[i]);
}

// Upload to VkBuffer with VK_BUFFER_USAGE_INDEX_BUFFER_BIT
```

**Note**: Use `uint32_t` indices if site count exceeds 65535, otherwise `uint16_t` saves bandwidth.

---

## Compute Shader: Vertex Normals

Compute vertex normals by averaging adjacent face normals.

### SSBO Layout

```glsl
// Binding 0: Face normals (input, one per triangle)
layout(std430, binding = 0) readonly buffer FaceNormals {
    vec4 faceNormals[];  // vec4 for alignment, .w unused
};

// Binding 1: Vertex adjacency (input, one per vertex)
layout(std430, binding = 1) readonly buffer Adjacency {
    VertexAdjacency adjacency[];
};

struct VertexAdjacency {
    uint triangleIndices[12];
    uint count;
    uint _padding[3];
};

// Binding 2: Vertex normals (output, one per vertex)
layout(std430, binding = 2) writeonly buffer VertexNormals {
    vec4 vertexNormals[];
};
```

### Compute Shader Code

```glsl
#version 450

layout(local_size_x = 256, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding = 0) readonly buffer FaceNormals {
    vec4 faceNormals[];
};

layout(std430, binding = 1) readonly buffer Adjacency {
    uint adjacencyData[];  // Packed: [tri0, tri1, ..., tri11, count, pad, pad, pad] per vertex
};

layout(std430, binding = 2) writeonly buffer VertexNormals {
    vec4 vertexNormals[];
};

layout(push_constant) uniform PushConstants {
    uint numVertices;
};

void main() {
    uint vertexId = gl_GlobalInvocationID.x;
    if (vertexId >= numVertices) return;
    
    // Read adjacency for this vertex
    uint baseIdx = vertexId * 16;  // 12 tris + count + 3 padding = 16 uints
    uint count = adjacencyData[baseIdx + 12];
    
    vec3 normalSum = vec3(0.0);
    
    for (uint i = 0; i < count; i++) {
        uint triId = adjacencyData[baseIdx + i];
        normalSum += faceNormals[triId].xyz;
    }
    
    // Normalize and store
    vec3 result = (length(normalSum) > 0.0001) 
        ? normalize(normalSum) 
        : vec3(0.0, 1.0, 0.0);
    
    vertexNormals[vertexId] = vec4(result, 0.0);
}
```

### Dispatch

```cpp
// After face normals are updated (initial generation or height modification)
void computeVertexNormals(VkCommandBuffer cmd, uint32_t numVertices) {
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, vertexNormalPipeline);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ...);
    
    // Push constant
    vkCmdPushConstants(cmd, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &numVertices);
    
    // Dispatch: ceil(numVertices / 256) workgroups
    uint32_t groupCount = (numVertices + 255) / 256;
    vkCmdDispatch(cmd, groupCount, 1, 1);
    
    // Memory barrier before vertex shader reads normals
    VkMemoryBarrier barrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER };
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_VERTEX_ATTRIBUTE_READ_BIT;
    vkCmdPipelineBarrier(cmd, 
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
        VK_PIPELINE_STAGE_VERTEX_INPUT_BIT,
        0, 1, &barrier, 0, nullptr, 0, nullptr);
}
```

---

## Collision Detection System

### Overview

The collision system must efficiently detect when an AABB (the player's fighter jet) intersects terrain triangles. The hierarchy is:

```
World Position
    │
    ▼
ChunkManager::getChunk()  ─────────────────────  O(1) hash lookup
    │
    ▼
SpatialGrid::getCandidateTriangles()  ─────────  O(1) grid cell lookup
    │
    ▼
AABB-Triangle Intersection Tests  ─────────────  O(k) where k = triangles in cells
    │
    ▼
Collision Response (optional physics)
```

### AABB-Terrain Collision Algorithm

For a fighter jet with AABB bounds, we need to detect if any part of the AABB penetrates the terrain surface.

#### Step 1: Determine Relevant Chunks

```cpp
struct AABB {
    float minX, minY, minZ;
    float maxX, maxY, maxZ;
};

std::vector<TerrainChunk*> getChunksOverlappingAABB(const AABB& aabb) {
    std::vector<TerrainChunk*> result;
    
    int minChunkX = static_cast<int>(std::floor(aabb.minX / chunkSize));
    int maxChunkX = static_cast<int>(std::floor(aabb.maxX / chunkSize));
    int minChunkZ = static_cast<int>(std::floor(aabb.minZ / chunkSize));
    int maxChunkZ = static_cast<int>(std::floor(aabb.maxZ / chunkSize));
    
    for (int cz = minChunkZ; cz <= maxChunkZ; cz++) {
        for (int cx = minChunkX; cx <= maxChunkX; cx++) {
            if (TerrainChunk* chunk = chunkManager.get({cx, cz})) {
                result.push_back(chunk);
            }
        }
    }
    return result;
}
```

#### Step 2: Get Candidate Triangles from SpatialGrid

```cpp
std::vector<uint32_t> getCandidateTriangles(const TerrainChunk* chunk, const AABB& aabb) {
    const SpatialGrid& grid = chunk->spatial;
    
    // Convert AABB to grid cells
    int cellMinX = static_cast<int>((aabb.minX - grid.minX) / grid.cellSize);
    int cellMaxX = static_cast<int>((aabb.maxX - grid.minX) / grid.cellSize);
    int cellMinZ = static_cast<int>((aabb.minZ - grid.minZ) / grid.cellSize);
    int cellMaxZ = static_cast<int>((aabb.maxZ - grid.minZ) / grid.cellSize);
    
    // Clamp to grid bounds
    cellMinX = std::max(0, cellMinX);
    cellMaxX = std::min(grid.gridW - 1, cellMaxX);
    cellMinZ = std::max(0, cellMinZ);
    cellMaxZ = std::min(grid.gridH - 1, cellMaxZ);
    
    // Collect unique triangles
    std::vector<uint32_t> candidates;
    std::unordered_set<uint32_t> seen;
    
    for (int cz = cellMinZ; cz <= cellMaxZ; cz++) {
        for (int cx = cellMinX; cx <= cellMaxX; cx++) {
            int cellIdx = cz * grid.gridW + cx;
            uint32_t start = grid.cellOffsets[cellIdx];
            uint32_t end = grid.cellOffsets[cellIdx + 1];
            
            for (uint32_t i = start; i < end; i++) {
                uint32_t triId = grid.cellTriangles[i];
                if (seen.insert(triId).second) {
                    candidates.push_back(triId);
                }
            }
        }
    }
    return candidates;
}
```

#### Step 3: AABB-Triangle Intersection Test

There are two approaches, depending on precision needs:

**Approach A: Height Sample Method (Fast, Sufficient for Most Cases)**

Sample terrain height at multiple points within the AABB footprint. If any sample point has terrain height above the AABB's minY, there's a collision.

```cpp
struct CollisionResult {
    bool collided;
    float penetrationDepth;
    vec3 contactNormal;      // Face normal at deepest penetration
    vec3 contactPoint;       // World position of deepest penetration
};

CollisionResult checkAABBTerrainCollision(const AABB& aabb, const TerrainChunk* chunk) {
    CollisionResult result = { false, 0.0f, {0,1,0}, {0,0,0} };
    
    // Sample grid within AABB footprint
    const float sampleSpacing = 0.5f;  // Adjust based on triangle size
    
    for (float z = aabb.minZ; z <= aabb.maxZ; z += sampleSpacing) {
        for (float x = aabb.minX; x <= aabb.maxX; x += sampleSpacing) {
            
            // Find triangle at this XZ position
            int triId;
            if (!chunk->spatial.locateTriangle(x, z, triId)) continue;
            
            // Interpolate terrain height
            float terrainHeight = interpolateHeight(chunk, triId, x, z);
            
            // Check if AABB bottom is below terrain
            if (aabb.minY < terrainHeight) {
                float depth = terrainHeight - aabb.minY;
                
                if (depth > result.penetrationDepth) {
                    result.collided = true;
                    result.penetrationDepth = depth;
                    result.contactNormal = getFaceNormal(chunk, triId);
                    result.contactPoint = { x, terrainHeight, z };
                }
            }
        }
    }
    
    // Also sample AABB corners specifically
    vec2 corners[4] = {
        {aabb.minX, aabb.minZ}, {aabb.maxX, aabb.minZ},
        {aabb.minX, aabb.maxZ}, {aabb.maxX, aabb.maxZ}
    };
    for (const auto& corner : corners) {
        // Same height check as above...
    }
    
    return result;
}
```

**Approach B: Separating Axis Theorem (SAT) - Precise**

For exact AABB-triangle intersection, use SAT. This is more expensive but handles edge cases (e.g., AABB corner grazing triangle edge).

```cpp
bool aabbIntersectsTriangle(const AABB& aabb, vec3 v0, vec3 v1, vec3 v2) {
    // Translate triangle to AABB center
    vec3 center = {
        (aabb.minX + aabb.maxX) * 0.5f,
        (aabb.minY + aabb.maxY) * 0.5f,
        (aabb.minZ + aabb.maxZ) * 0.5f
    };
    vec3 halfSize = {
        (aabb.maxX - aabb.minX) * 0.5f,
        (aabb.maxY - aabb.minY) * 0.5f,
        (aabb.maxZ - aabb.minZ) * 0.5f
    };
    
    v0 = v0 - center;
    v1 = v1 - center;
    v2 = v2 - center;
    
    // Edge vectors
    vec3 e0 = v1 - v0;
    vec3 e1 = v2 - v1;
    vec3 e2 = v0 - v2;
    
    // Test 13 separating axes:
    // - 3 AABB face normals (X, Y, Z axes)
    // - 1 triangle face normal
    // - 9 cross products of AABB edges x triangle edges
    
    // Axis: AABB X
    if (!testAxisAABB(v0.x, v1.x, v2.x, halfSize.x)) return false;
    // Axis: AABB Y
    if (!testAxisAABB(v0.y, v1.y, v2.y, halfSize.y)) return false;
    // Axis: AABB Z
    if (!testAxisAABB(v0.z, v1.z, v2.z, halfSize.z)) return false;
    
    // Axis: Triangle normal
    vec3 triNormal = cross(e0, e1);
    if (!testAxisTriangleNormal(v0, triNormal, halfSize)) return false;
    
    // 9 cross-product axes (AABB edges x triangle edges)
    // ... (omitted for brevity, but required for complete SAT)
    
    return true;  // No separating axis found = intersection
}
```

#### Step 4: Collision Response Data

For a fighter jet crash, you need:

| Data | Source | Purpose |
|------|--------|---------|
| `penetrationDepth` | Height difference | How deep into terrain |
| `contactNormal` | `FaceNormals[triId]` | Impact angle, explosion direction |
| `contactPoint` | Interpolated position | Spawn explosion VFX here |

### Complete Collision Query Function

```cpp
CollisionResult checkFighterJetCollision(const AABB& jetAABB) {
    CollisionResult result = { false, 0.0f, {0,1,0}, {0,0,0} };
    
    // Get overlapping chunks
    auto chunks = getChunksOverlappingAABB(jetAABB);
    
    for (const TerrainChunk* chunk : chunks) {
        // Get candidate triangles
        auto candidates = getCandidateTriangles(chunk, jetAABB);
        
        for (uint32_t triId : candidates) {
            // Get triangle vertices in 3D
            vec3 v0 = getTriangleVertex(chunk, triId, 0);
            vec3 v1 = getTriangleVertex(chunk, triId, 1);
            vec3 v2 = getTriangleVertex(chunk, triId, 2);
            
            // Option A: Fast height-based check
            // (Check if AABB footprint overlaps triangle XZ projection,
            //  then compare heights)
            
            // Option B: Precise SAT
            if (aabbIntersectsTriangle(jetAABB, v0, v1, v2)) {
                // Compute penetration depth and contact info
                // Update result if this is deeper penetration
            }
        }
    }
    
    return result;
}
```

---

## Projectile Collision at Distance

For long-range projectile hits (player shooting enemies 1km+ away), the approach differs from AABB collision.

### Projectile vs Entity Collision

Projectiles hitting distant enemies don't need terrain collision at the enemy's location - they need **ray vs entity AABB** tests.

```
Projectile Ray: origin + t * direction
    │
    ▼
Spatial partitioning of entities (separate from terrain)
    │
    ▼
Ray-AABB intersection tests for candidate entities
```

### Projectile vs Terrain Collision

For projectiles that might hit terrain (missed shots, ground-targeted weapons):

**Option A: Raycast (Continuous)**

Step along the ray, checking terrain height at intervals.

```cpp
bool raycastTerrain(vec3 origin, vec3 direction, float maxDist, RayHit& outHit) {
    const float stepSize = 1.0f;  // Tune based on terrain roughness
    
    for (float t = 0; t < maxDist; t += stepSize) {
        vec3 pos = origin + direction * t;
        
        TerrainChunk* chunk = chunkManager.getChunk(pos.x, pos.z);
        if (!chunk) continue;
        
        float terrainHeight;
        if (chunk->spatial.sampleHeight(pos.x, pos.z, terrainHeight)) {
            if (pos.y <= terrainHeight) {
                // Hit terrain - refine with binary search
                outHit = refineRayHit(origin, direction, t - stepSize, t, chunk);
                return true;
            }
        }
    }
    return false;
}
```

**Option B: DDA Grid Traversal (Faster for Long Rays)**

Use Digital Differential Analyzer to step through spatial grid cells along the ray.

```cpp
bool raycastTerrainDDA(vec3 origin, vec3 direction, float maxDist, RayHit& outHit) {
    // Initialize DDA for chunk grid
    // Step through chunks along ray
    // Within each chunk, DDA through spatial grid cells
    // Test triangles in each cell
    
    // (Implementation is more complex but O(cells traversed) instead of O(distance))
}
```

### Tiered Collision for Distant Entities

For entities (not projectiles), tiered collision still makes sense:

| Distance | Collision Frequency | Method |
|----------|--------------------| -------|
| Near (< 100m) | Every frame | Full AABB-terrain |
| Medium (100-500m) | Every 4-8 frames | Full AABB-terrain, staggered |
| Far (> 500m) | Every 16+ frames | Simple height sample at center |

**Key insight**: Projectiles use raycasting (continuous along path), entities use AABB collision (discrete positions). Different problems.

---

## Batch Collision Queries

For thousands of entities, batch processing improves cache utilization.

### Batch Query Structure

```cpp
struct BatchCollisionQuery {
    // Input: entity AABBs (SoA)
    std::vector<float> aabbMinX, aabbMinY, aabbMinZ;
    std::vector<float> aabbMaxX, aabbMaxY, aabbMaxZ;
    
    // Output: collision results (SoA)
    std::vector<bool>  collided;
    std::vector<float> penetrationDepth;
    std::vector<float> normalX, normalY, normalZ;
};
```

### Batch Processing Strategy

**Step 1: Sort by Chunk**

Group entities by which chunk(s) they overlap. This ensures we load each chunk's data once.

```cpp
void batchCollisionCheck(BatchCollisionQuery& query) {
    // Build chunk -> entity mapping
    std::unordered_map<ChunkCoord, std::vector<uint32_t>> chunkEntities;
    
    for (uint32_t i = 0; i < query.size(); i++) {
        AABB aabb = getAABB(query, i);
        auto chunks = getChunkCoordsOverlappingAABB(aabb);
        for (const auto& coord : chunks) {
            chunkEntities[coord].push_back(i);
        }
    }
    
    // Process each chunk's entities together
    for (const auto& [coord, entityIndices] : chunkEntities) {
        TerrainChunk* chunk = chunkManager.get(coord);
        if (!chunk) continue;
        
        // Prefetch chunk data into cache
        prefetchChunkData(chunk);
        
        for (uint32_t entityIdx : entityIndices) {
            // Process collision for this entity
            // Chunk data is now hot in cache
        }
    }
}
```

**Step 2: SIMD Where Possible**

Height interpolation and point-in-triangle tests can be vectorized.

```cpp
// Process 8 entities at once using AVX2
void checkHeightsAVX(
    const float* posX, const float* posZ,
    const float* aabbMinY,
    float* outCollided,  // 1.0 = collided, 0.0 = no
    uint32_t count
) {
    for (uint32_t i = 0; i < count; i += 8) {
        __m256 px = _mm256_load_ps(&posX[i]);
        __m256 pz = _mm256_load_ps(&posZ[i]);
        __m256 minY = _mm256_load_ps(&aabbMinY[i]);
        
        // Sample terrain heights (this part is harder to vectorize due to lookups)
        __m256 terrainH = sampleHeightsBatch(px, pz);  // Custom implementation
        
        // Compare
        __m256 collided = _mm256_cmp_ps(minY, terrainH, _CMP_LT_OQ);
        _mm256_store_ps(&outCollided[i], collided);
    }
}
```

---

## Cache Coherence Guidelines

### Memory Access Patterns

| Pattern | Cache Behavior | When to Use |
|---------|---------------|-------------|
| Sequential array access | Excellent (prefetcher works) | Iterating all triangles, all sites |
| Strided access | Good if stride is small | Accessing every Nth element |
| Random access | Poor (cache misses) | Avoid if possible |
| Pointer chasing | Very poor | Avoid; use indices instead |

### Data Locality Strategies

**1. Keep hot data together**

```cpp
// Good: Collision query only touches what it needs
struct CollisionData {
    // Hot: accessed every query
    std::vector<float> heights;
    SpatialGrid grid;
    
    // Warm: accessed for matching triangles
    std::vector<float> faceNormalX, faceNormalY, faceNormalZ;
    std::vector<uint32_t> triSiteA, triSiteB, triSiteC;
    
    // Cold: rarely accessed during collision
    std::vector<int32_t> halfEdgeTwin;  // Topology queries only
};
```

**2. Avoid interleaving mutable and immutable data**

```cpp
// Bad: mutable height mixed with immutable position
struct Vertex {
    float x, y, z;     // x, z immutable; y mutable
    float nx, ny, nz;  // mutable
};

// Good: separate arrays
std::vector<float> posX, posZ;  // Immutable after generation
std::vector<float> height;       // Mutable
std::vector<float> normalX, normalY, normalZ;  // Mutable
```

**3. Process in chunk order**

```cpp
// Bad: random chunk access
for (Entity& e : entities) {
    TerrainChunk* chunk = getChunk(e.pos);  // Random access
    checkCollision(e, chunk);
}

// Good: sorted by chunk
std::sort(entities.begin(), entities.end(), [](auto& a, auto& b) {
    return getChunkCoord(a.pos) < getChunkCoord(b.pos);
});
for (Entity& e : entities) {
    // Sequential entities likely share same chunk (hot in cache)
}
```

### Prefetching

```cpp
// Explicit prefetch for known access pattern
void processTriangles(const TrianglesData& tris, uint32_t start, uint32_t end) {
    for (uint32_t i = start; i < end; i++) {
        // Prefetch next iteration's data
        if (i + 8 < end) {
            _mm_prefetch(&tris.siteA[i + 8], _MM_HINT_T0);
            _mm_prefetch(&tris.siteB[i + 8], _MM_HINT_T0);
            _mm_prefetch(&tris.siteC[i + 8], _MM_HINT_T0);
        }
        
        // Process current triangle
        uint32_t a = tris.siteA[i];
        uint32_t b = tris.siteB[i];
        uint32_t c = tris.siteC[i];
        // ...
    }
}
```

---

## Runtime Terrain Modification

Since terrain heights can change at runtime (e.g., meteor crater), gradients are NOT precomputed.

### What Must Update When Heights Change

| Data | Update Required | Method |
|------|-----------------|--------|
| `heights[]` | Direct write | Modify values in crater radius |
| `FaceNormals[]` | Recompute | For affected triangles only |
| `VertexNormals[]` | Recompute | Dispatch compute shader |
| Spatial grid | No change | Grid structure is based on XZ, not height |
| Barycentric denom | No change | Based on XZ positions only |

### Efficient Partial Update

```cpp
void applyMeteorCrater(vec3 impactPos, float radius, float depth) {
    // 1. Modify heights in radius
    std::vector<uint32_t> affectedSites;
    std::vector<uint32_t> affectedTriangles;
    
    for (uint32_t i = 0; i < numSites; i++) {
        float dx = sites.posX[i] - impactPos.x;
        float dz = sites.posZ[i] - impactPos.z;
        float dist = std::sqrt(dx*dx + dz*dz);
        
        if (dist < radius) {
            // Crater profile (e.g., parabolic depression)
            float t = dist / radius;
            float depression = depth * (1.0f - t*t);
            sites.height[i] -= depression;
            
            affectedSites.push_back(i);
        }
    }
    
    // 2. Find affected triangles
    for (uint32_t i = 0; i < numTriangles; i++) {
        if (siteAffected(tris.siteA[i]) || 
            siteAffected(tris.siteB[i]) || 
            siteAffected(tris.siteC[i])) {
            affectedTriangles.push_back(i);
        }
    }
    
    // 3. Recompute face normals for affected triangles only
    for (uint32_t triId : affectedTriangles) {
        recomputeFaceNormal(triId);
    }
    
    // 4. Upload face normals to GPU, dispatch compute shader for vertex normals
    updateFaceNormalsBuffer(affectedTriangles);
    dispatchVertexNormalCompute();
    
    // 5. Update vertex buffer heights
    updateHeightBuffer(affectedSites);
}
```

### Why No Gradient Precomputation

The gradient inverse optimization precomputes:
```
invM = inverse([[ab.x, ab.y], [ac.x, ac.y]]) / det
```

This is based on vertex **XZ positions**, which are indeed immutable. However:

1. The optimization saves ~10 ops per gradient query
2. Each triangle would need 4 extra floats stored
3. Memory bandwidth for loading those floats may negate the compute savings
4. For your use case (collision detection, not per-frame gradient queries for every droplet), the base `TriangleGradient` computation is fast enough

**Decision**: Skip gradient precomputation. Use `TriangleGradient()` as-is. If profiling later shows it's a bottleneck, the optimization is straightforward to add.

---

## Summary

| Component | Strategy |
|-----------|----------|
| Data layout | SoA for cache coherence, SoAoS for VertexAdjacency |
| Vertex normals | Compute shader averaging face normals |
| Index buffers | Direct from Tris.siteA/B/C |
| Collision hierarchy | ChunkManager → SpatialGrid → Triangle tests |
| AABB collision | Height sampling (fast) or SAT (precise) |
| Projectile collision | Raycast with DDA grid traversal |
| Distant entities | Tiered frequency, not tiered precision |
| Batch queries | Sort by chunk, SIMD where possible |
| Runtime modification | Partial updates to heights, face normals, vertex normals |
