package hydro

import (
	"math/rand"
	"procedural_generation/terrain_generation/duals2/core"
)

// -------------------------------------------------------------------
// Ocean Region System
// -------------------------------------------------------------------

// OceanSeedPoint represents a potential ocean origin.
type OceanSeedPoint struct {
	ChunkX int
	ChunkZ int
	Seed   int64
}

// GenerateOceanRegion creates an ocean region starting from a seed chunk.
// The ocean expands from the seed up to the configured maximum size.
func (hm *HydroManager) GenerateOceanRegion(seedChunk core.ChunkCoord, worldSeed int64) *core.OceanRegion {
	// Create deterministic RNG for this ocean
	oceanSeed := combineSeeds(worldSeed, int64(seedChunk.X), int64(seedChunk.Z))
	rng := rand.New(rand.NewSource(oceanSeed))

	// Check probability
	if rng.Float64() > hm.Cfg.OceanProbability {
		return nil
	}

	// Determine ocean extent (random within max bounds)
	extentX := rng.Intn(hm.Cfg.OceanMaxChunksX) + 1
	extentZ := rng.Intn(hm.Cfg.OceanMaxChunksZ) + 1

	// Center the ocean on the seed chunk (with some randomness)
	offsetX := rng.Intn(extentX)
	offsetZ := rng.Intn(extentZ)

	ocean := &core.OceanRegion{
		ID:        hm.NextOceanID,
		MinChunkX: seedChunk.X - offsetX,
		MinChunkZ: seedChunk.Z - offsetZ,
		MaxChunkX: seedChunk.X - offsetX + extentX - 1,
		MaxChunkZ: seedChunk.Z - offsetZ + extentZ - 1,
		SeaLevel:  hm.Cfg.SeaLevel,
	}

	// Register chunks as ocean
	for x := ocean.MinChunkX; x <= ocean.MaxChunkX; x++ {
		for z := ocean.MinChunkZ; z <= ocean.MaxChunkZ; z++ {
			hm.OceanChunks[core.ChunkCoord{X: x, Z: z}] = ocean.ID
		}
	}

	hm.Oceans = append(hm.Oceans, *ocean)
	hm.NextOceanID++

	return ocean
}

// IsOceanChunk returns true if the chunk is part of an ocean region.
func (hm *HydroManager) IsOceanChunk(coord core.ChunkCoord) bool {
	_, ok := hm.OceanChunks[coord]
	return ok
}

// GetOceanForChunk returns the ocean region for a chunk, or nil if not ocean.
func (hm *HydroManager) GetOceanForChunk(coord core.ChunkCoord) *core.OceanRegion {
	oceanID, ok := hm.OceanChunks[coord]
	if !ok {
		return nil
	}
	for i := range hm.Oceans {
		if hm.Oceans[i].ID == oceanID {
			return &hm.Oceans[i]
		}
	}
	return nil
}

// ComputeOceanChunkData calculates ocean data for a specific chunk.
func (hm *HydroManager) ComputeOceanChunkData(
	coord core.ChunkCoord,
	mesh *core.DelaunayMesh,
	heights []float64,
) core.OceanChunkData {
	ocean := hm.GetOceanForChunk(coord)
	if ocean == nil {
		return core.OceanChunkData{IsOcean: false, OceanID: -1}
	}

	data := core.OceanChunkData{
		IsOcean:       true,
		OceanID:       ocean.ID,
		SubmergedTris: make([]int, 0),
		CoastlinePts:  make([]core.Vec2, 0),
	}

	// Find submerged triangles and coastline
	for triID, tri := range mesh.Tris {
		hA := heights[tri.A]
		hB := heights[tri.B]
		hC := heights[tri.C]

		belowA := hA <= ocean.SeaLevel
		belowB := hB <= ocean.SeaLevel
		belowC := hC <= ocean.SeaLevel

		if belowA && belowB && belowC {
			// Fully submerged
			data.SubmergedTris = append(data.SubmergedTris, triID)
		} else if belowA || belowB || belowC {
			// Partial - coastline passes through this triangle
			// Record the crossing points for coastline rendering
			data.SubmergedTris = append(data.SubmergedTris, triID)

			// Find edge crossings
			pA := mesh.Sites[tri.A].Pos
			pB := mesh.Sites[tri.B].Pos
			pC := mesh.Sites[tri.C].Pos

			if belowA != belowB {
				t := (ocean.SeaLevel - hA) / (hB - hA)
				crossing := core.Vec2{X: pA.X + t*(pB.X-pA.X), Y: pA.Y + t*(pB.Y-pA.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
			if belowB != belowC {
				t := (ocean.SeaLevel - hB) / (hC - hB)
				crossing := core.Vec2{X: pB.X + t*(pC.X-pB.X), Y: pB.Y + t*(pC.Y-pB.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
			if belowC != belowA {
				t := (ocean.SeaLevel - hC) / (hA - hC)
				crossing := core.Vec2{X: pC.X + t*(pA.X-pC.X), Y: pC.Y + t*(pA.Y-pC.Y)}
				data.CoastlinePts = append(data.CoastlinePts, crossing)
			}
		}
	}

	return data
}
