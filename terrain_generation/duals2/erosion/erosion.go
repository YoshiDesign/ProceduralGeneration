package eros

/** Some Notes:

### Direction Calculation
- In [Vol] a random direction is chosen if the new directions magnitude is below a
threshold. This happens most often when a drop is initialized inside a flat cell. My
algorithm only chooses a random direction if the direction otherwise would be 0. The
direction of motion is then normalized.

This direction is a blended value between g and dirold defined by a parameter
pinertia with a value between 0 and 1. 1 means g is not taken into account and the
direction never changes and 0 means the previous direction dirold is ignored and the
new direction is along the negative gradient.

dirNew = dirOld · pinertia − g · (1 − pinertia)
posNew += dirNew


*/

import (
	"fmt"
	"math/rand"
	"procedural_generation/terrain_generation/duals2/core"
)

type ErosionManager struct {
	ErosionHeightDeltas map[core.ChunkCoord][/*core.SiteIndex*/]float64
	Droplets []Droplet
	cfg ErosionConfig
	ChunkConfig core.ChunkConfig
}

type ErosionConfig struct {
	Pinertia float64 // Must always be clamped between 0 and 1
	PCapacity float64
	PDeposition float64
	PErosion float64
	PEvaporation float64
	PMinSlope float64
	Gravity float64
	NumDroplets int
	NumSteps int
}

type Droplet struct {
	Pos core.Vec2
	Vel core.Vec2
	Dir core.Vec2 // According to gradient
	Capacity float64
	Water float64
	Sediment float64
}

func (em *ErosionManager) NewErosionManager(chunkConfig core.ChunkConfig) *ErosionManager {
	return &ErosionManager{
		ErosionHeightDeltas: make(map[core.ChunkCoord][/*core.SiteIndex*/]float64), // Erosion output for every chunk generated
		cfg: ErosionConfig{
			Pinertia: 0.1,
			PCapacity: 8.0,
			PDeposition: 0.1,
			PErosion: 0.1,
			PEvaporation: 0.01,
			PMinSlope: -0.001,
			Gravity: 0.1,
			NumDroplets: 2000,
			NumSteps: 20,
		},
		ChunkConfig: chunkConfig,
	}
}

func (em *ErosionManager) HydraulicErosion(chunk *core.TerrainChunk, sg *core.SpatialGrid) error {

	// The same seed that generated the blue noise for this chunk
	rng := rand.New(rand.NewSource(chunk.Seed))

	// To be our erosion output
	em.ErosionHeightDeltas[chunk.Coord] = make([/*core.SiteIndex*/]float64, len(em.ErosionHeightDeltas))

	// Initialize
	for i := 0; i < em.cfg.NumDroplets; i++ {
		// Init new droplet at a random position in world space within the chunk's core + halo
		droplet := em.NewDroplet(
			rng.Intn(int(chunk.Cfg.ChunkSize) + int(chunk.Cfg.HaloWidth)), 	// x
			rng.Intn(int(chunk.Cfg.ChunkSize) + int(chunk.Cfg.HaloWidth))) 	// z

		// Get the triangle this droplet "fell" on
		em.Droplets = append(em.Droplets, droplet)

	}

	// Simulate
	for i := 0; i < em.cfg.NumSteps; i++ {

		/*
			Note to self - You're using the Barycentric weights here across multiple
			places in the code. In the production version, create function signatures
			that accept the weights in the event that you already have them available
			for any given triangle. Think of how this situation might occur... 
		*/

		for di := range em.Droplets {

			droplet := &em.Droplets[di]

			t1Idx, ok := sg.LocateTriangle(droplet.Pos.X, droplet.Pos.Y)
			if !ok {
				fmt.Println("Droplet out of bounds: ", droplet.Pos.X, ", ", droplet.Pos.Y)
				continue
			}
			// Triangle under droplet's initial position
			t1 := chunk.Mesh.Tris[t1Idx]
			siteA := chunk.Mesh.Sites[t1.A]
			siteB := chunk.Mesh.Sites[t1.B]
			siteC := chunk.Mesh.Sites[t1.C]

			// Height at initial position
			h1, ok := sg.SampleHeight(droplet.Pos.X, droplet.Pos.Y)
			if !ok {
				fmt.Println("Height not found for triangle 1: ", t1Idx)
				continue
			}
			// x (u) and z (v) - constant gradients of the triangle
			u, v, ok := chunk.Mesh.TriangleGradient(t1Idx, chunk.Heights)
			if !ok {
				fmt.Println("Gradient not found for triangle 1: ", t1Idx)
				continue
			}
			
			// Pack it
			grad := core.Vec2{X: -u, Y: -v}

			// Move it
			em.MoveDroplet(droplet, grad)
						
			// Bop it
			droplet.Pos = droplet.Pos.Add(droplet.Dir)
			
			// Get the triangle at the new position (might be the same as the initial triangle)
			t2Idx, ok := sg.LocateTriangle(droplet.Pos.X, droplet.Pos.Y)
			if !ok {
				fmt.Println("Droplet out of bounds: ", droplet.Pos.X, ", ", droplet.Pos.Y)
				continue
			}
			// Height at the new position
			h2, ok := sg.SampleHeight(droplet.Pos.X, droplet.Pos.Y)
			if !ok {
				fmt.Println("Height not found for triangle 2: ", t2Idx)
				continue
			}

			// Downhill movement (erosion)
			if h1 > h2 {
				em.Erode(chunk.Coord, siteA, droplet, t1.A)
				em.Erode(chunk.Coord, siteB, droplet, t1.B)
				em.Erode(chunk.Coord, siteC, droplet, t1.C)
			} else { // Uphill (deposition)
				em.Deposit(chunk.Coord, siteA, droplet, t1.A)
				em.Deposit(chunk.Coord, siteB, droplet, t1.B)
				em.Deposit(chunk.Coord, siteC, droplet, t1.C)
			}
	
		}

	}

	return nil
}

func (em *ErosionManager) MoveDroplet(droplet *Droplet, _grad core.Vec2) error {

	newDirection := func (dirOld, grad core.Vec2) core.Vec2 {
		dirNew := dirOld.Mul(em.cfg.Pinertia).Sub(grad.Mul(1 - em.cfg.Pinertia))
		if dirNew.Len2() < 1e-20 {
			// Either keep old dir, or randomize a tiny direction
			return dirOld // or some fallback
		}
		return dirNew.Normalize()
	}

	droplet.Dir = newDirection(droplet.Dir, _grad)

	return nil
	
}

func (em *ErosionManager) Erode(chunk_coord core.ChunkCoord, site core.Site, droplet *Droplet, idx core.SiteIndex) {
	em.ErosionHeightDeltas[chunk_coord][idx] += 123 // TODO some result
}

func (em *ErosionManager) Deposit(chunk_coord core.ChunkCoord,site core.Site, droplet *Droplet, idx core.SiteIndex) {
	
}

func (em *ErosionManager) NewDroplet(x, z int) Droplet {

	pos := core.Vec2{
		X: float64(x),
		Y: float64(z),
	}

	return Droplet{
		Pos: pos,
		Capacity: em.cfg.PCapacity,
		Vel: core.Vec2{},
		Dir: core.Vec2{},
		Water: 0.0,
		Sediment: 0.0,
	}
}
