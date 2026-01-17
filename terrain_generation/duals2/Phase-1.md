Core Data Structures
Water Body Types:
RiverSegment — Path vertices with width/depth arrays, ordered upstream to downstream
River — Complete river with segments, source type, tributaries, and destination reference
Lake — Filled depression with site indices, water level, spillway, and flow connections
OceanRegion — Multi-chunk ocean defined by bounding chunk coordinates
Configuration:
HydroConfig — All tunable parameters with DefaultHydroConfig() providing sensible defaults
Per-Chunk Data:
ChunkHydroData — Rivers, lakes, ocean data, plus entry/exit points for cross-chunk continuity
River System
Path Tracing Algorithm:
TraceRiver(mesh, heights, startSite, bounds, width, depth, distance)
Uses gradient descent — follows lower neighbors
Applies equator-based flow bias — prevents contradicting flow directions
Adds meandering noise — ±5% random score adjustment
Tracks width/depth growth — rivers widen and deepen as they flow
Detects chunk boundary exits — caches exit points for neighbor propagation
Source Detection:
Sources must be above RiverSourceMinElev (default: 50 units)
Sources must be near local maxima (≤1/3 neighbors higher)
Sources pass probability check (RiverSourceProbability: 30% default)
Cross-Chunk Continuity:
RiverExitPoint stores position, width, depth, flow direction at boundary
HydroManager.riverPropagation map queues entries for destination chunks
When chunk B generates, it picks up rivers from GetRiverEntries(B)
Lake System
Depression Filling Algorithm:
FloodFillLake(mesh, heights, minimum)
Starts at a local minimum (vertex lower than all neighbors)
Flood-fills upward through connected low vertices
Identifies spillway (lowest exit point)
Returns Lake with boundary vertices, water level, depth
Constraints:
LakeMinDepth: 2 units minimum depression
LakeMinArea: 3+ vertices required
LakeProbability: 50% of eligible depressions become lakes
Ocean System
How Ocean Area Is Implemented:
type OceanRegion struct {    ID        int    MinChunkX int     // Inclusive bounding box    MinChunkZ int    MaxChunkX int     // Inclusive    MaxChunkZ int    SeaLevel  float64}
Generation (GenerateOceanRegion):
Uses deterministic seed from world seed + chunk coordinates
Checks OceanProbability (default: 10%)
Randomly determines extent within OceanMaxChunksX × OceanMaxChunksZ (default: 8×8)
Centers on seed chunk with random offset
Registers all chunks in the region to oceanChunks map
Per-Chunk Ocean Data (ComputeOceanChunkData):
SubmergedTris: Triangles where all vertices are below sea level
CoastlinePts: Edge crossing points where sea level intersects triangles (for rendering)
Key Design Point: Oceans only appear in chunks registered as ocean regions. Terrain can go below sea level elsewhere without becoming ocean — this is intentional per your requirement that oceans require the ocean biome designation.
Integration with Chunk System
TerrainChunk now includes:
Hydro *ChunkHydroData
ChunkManager now includes:
hydro *HydroManager
Generation Flow (in generateChunkInternal):
Generate terrain mesh (unchanged)
New: Call generateChunkHydrology(chunk) which:
Processes incoming rivers from neighbors
Finds new river sources
Traces all rivers
Detects lakes at local minima
Computes ocean data if in ocean region
Biome Integration (Future Refinement)
The current implementation uses static probabilities for hydrology formation:
Parameter	Default	Purpose
RiverSourceProbability	0.3	Chance a valid high-elevation point spawns a river
LakeProbability	0.5	Chance an eligible depression becomes a lake
OceanProbability	0.1	Chance a seed point generates an ocean region
To integrate biomes, you would modify these checks to query biome data:
// Example future integration:func (hm *HydroManager) shouldSpawnRiver(site int, biome BiomeType) bool {    switch biome {    case BiomeMountainWet, BiomeMountainCold:        return hm.rng.Float64() < 0.6  // Higher probability    case BiomeDesert:        return hm.rng.Float64() < 0.05 // Very low    default:        return hm.rng.Float64() < hm.cfg.RiverSourceProbability    }}
The HydroConfig parameters provide the baseline that biome modifiers would adjust.


River Rendering Detail
Rivers are rendered as connected line segments with width:

Use vector.StrokeLine with width from `segment.Widths[i]`
Color gradient from light blue (narrow) to deep blue (wide)
Draw flow direction indicators optionally
Lake Rendering Detail
Lakes are rendered by:

Finding triangles whose vertices are all in lake.SiteIndices
Filling those triangles with water color
Drawing a subtle outline for the shore
Ocean Rendering Detail
Oceans use pre-computed OceanChunkData:

Fill SubmergedTris with ocean color (darker for deeper)
Draw CoastlinePts as a connected line or individual markers