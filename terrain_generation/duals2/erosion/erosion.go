package eros

/** Some Notes:

- The dimensionality of droplet.Sediment/Capacity implies height units, not an abstract unit of sediment.

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

TODO: We're not yet taking into account an erosion radius.
	  We just erode individual points without considering neighboring points.
	  Note that this only applies to erosion, not deposition.

*/

import (
	"fmt"
	"math"
	"math/rand"
	"procedural_generation/terrain_generation/duals2/core"
)

type ErosionManager struct {
	ErosionHeightDeltas map[core.ChunkCoord]map[core.SiteIndex]float64
	Droplets []Droplet
	cfg ErosionConfig
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
	Vel float64
	Dir core.Vec2 // According to gradient
	Capacity float64
	Water float64
	Sediment float64
	IsAlive bool
}

func NewErosionManager() *ErosionManager {
	return &ErosionManager{
		ErosionHeightDeltas: make(map[core.ChunkCoord]map[core.SiteIndex]float64), // Erosion output for every chunk generated
		cfg: HeavyErosion(),
	}
}

func (em *ErosionManager) HydraulicErosion(chunk *core.TerrainChunk, sg *core.SpatialGrid, chunkSeed int64) error {

	oob := 0

	em.ClearHeightDeltas()
	em.ClearDroplets()

	// The same seed that generated the blue noise for this chunk
	rng := rand.New(rand.NewSource(chunk.Seed))

	// To be our erosion output
	em.ErosionHeightDeltas[chunk.Coord] = make(map[core.SiteIndex]float64)

	// Initialize
	for i := 0; i < em.cfg.NumDroplets; i++ {
		// Init new droplet at a random position in world space within the chunk's core + halo
		droplet := em.NewDroplet(
			rng.Intn(int(chunk.Cfg.ChunkSize)), 	// x
			rng.Intn(int(chunk.Cfg.ChunkSize)), 	// z
			rng) // seeded rng for deterministic direction

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

			if !em.Droplets[di].IsAlive {
				continue
			}

		droplet := &em.Droplets[di]

		// Bounds check on local position before triangle lookup
		chunkSize := chunk.Cfg.ChunkSize
		if droplet.Pos.X < 0 || droplet.Pos.X >= chunkSize ||
			droplet.Pos.Y < 0 || droplet.Pos.Y >= chunkSize {
			oob++
			em.KillDroplet(droplet)
			continue
		}

		worldPos1 := core.Vec2{X: droplet.Pos.X + chunk.MinX, Y: droplet.Pos.Y + chunk.MinZ}

		// Explicitly locate the triangle so we can have site indices available here for the erosion/deposition operations
		t1Idx, ok := sg.LocateTriangle(worldPos1.X, worldPos1.Y)
		if !ok {
			oob++
			em.KillDroplet(droplet)
			continue
		}

		// Triangle & its Sites (vertices) at droplet's initial position
		t1 := chunk.Mesh.Tris[t1Idx]

		// Weights for the triangle at the initial position - these will be used for the erosion/deposition operations
		// Use world coordinates since the mesh is in world space
		wa, wb, wc, ok := sg.Mesh.Barycentric(t1Idx, worldPos1)
		if !ok {
			return fmt.Errorf("barycentric weights not found!")
		}


		// Height at initial position (use world coordinates)
		h1, ok := sg.Mesh.SampleScalar(t1Idx, worldPos1, chunk.Heights)
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
			
			// Pack it - Note: This value will be normalized so every droplet moves the same amount per step, regardless of its velocity
			grad := core.Vec2{X: -u, Y: -v}

			// Move it
			em.MoveDroplet(droplet, grad)
	
		// Bop it (record the update)
		droplet.Pos = droplet.Pos.Add(droplet.Dir)

		// Bounds check after movement
		if droplet.Pos.X < 0 || droplet.Pos.X >= chunkSize ||
			droplet.Pos.Y < 0 || droplet.Pos.Y >= chunkSize {
			oob++
			em.KillDroplet(droplet)
			continue
		}

		worldPos2 := core.Vec2{X: droplet.Pos.X + chunk.MinX, Y: droplet.Pos.Y + chunk.MinZ}

		// Get the new position
		t2Idx, ok := sg.LocateTriangle(worldPos2.X, worldPos2.Y)
		if !ok {
			oob++
			em.KillDroplet(droplet)
			continue
		}

			t2 := chunk.Mesh.Tris[t2Idx]

		// NOTE: SampleHeight is the same as SampleScalar (w/ heights) if we don't care to keep track of the TriID.
		// New Height after movement (use world coordinates)
		h2, ok := sg.Mesh.SampleScalar(t2Idx, worldPos2, chunk.Heights)
		if !ok {
			fmt.Println("Height not found for triangle 2")
			continue
		}

		// Barycentric weights for the new position (use world coordinates)
		wa2, wb2, wc2, ok := sg.Mesh.Barycentric(t2Idx, worldPos2)
			if !ok {
				fmt.Println("Barycentric weights not found for triangle 2: ", t2Idx)
				continue
			}

			hdiff := h2 - h1

			// Downhill movement (erosion routine)
			if hdiff < 0.0 {
				em.Erode(
					chunk.Coord, 
					[]core.SiteIndex{t1.A, t1.B, t1.C}, // Invariant ordering
					[]float64{wa, wb, wc},
					droplet, 
					hdiff,
				)
			} else if hdiff > 0.0 { // Uphill (deposition routine)
				em.UphillDeposit(
					chunk.Coord, 
					[]core.SiteIndex{t1.A, t1.B, t1.C}, 
					[]core.SiteIndex{t2.A, t2.B, t2.C},
					[]float64{wa, wb, wc},
					[]float64{wa2, wb2, wc2},
					[]float64{chunk.Heights[t2.A], chunk.Heights[t2.B], chunk.Heights[t2.C]},
					droplet, 
					hdiff,
				)
			} else { // Flat - no change
				continue
			}

			em.UpdateVelocity(droplet, hdiff)
			em.Evaporate(droplet)
	
		}

	}

	fmt.Println("Done: Out of bounds: ", oob)
	return nil
}

func (em *ErosionManager) MoveDroplet(droplet *Droplet, _grad core.Vec2) {

	newDirection := func (dirOld, grad core.Vec2) core.Vec2 {
		dirNew := dirOld.Mul(em.cfg.Pinertia).Sub(grad.Mul(1 - em.cfg.Pinertia))
		if dirNew.Len2() < 1e-20 {
			// Either keep old dir, or randomize a tiny direction
			return dirOld // or some fallback
		}
		return dirNew.Normalize() // Always has magnitude 1
	}

	droplet.Dir = newDirection(droplet.Dir, _grad)
	
}

// Erosion, but also covers conditional deposition
func (em *ErosionManager) Erode(
	chunk_coord core.ChunkCoord, // So we can index into chunks -> height deltas. Alternatively we could pass in the current ChunkCoord being processed from the beginning
	siteIndices []core.SiteIndex, 
	w []float64,
	droplet *Droplet, 
	hdiff float64,
	) {

	// Update the carry-capacity of the droplet: c = max(−hdiff, pminSlope) · vel · water · pcapacity
	/*
	 	If the height difference converges to 0 the capacity also would converge to 0, which
		leads to less erosion and more deposition in flatter areas. Sometimes however it is more
		aesthetic to erode flatter terrain too.
	*/
	// Note: Capacity will always end up being reduced
	slopeTerm := math.Max(-hdiff, em.cfg.PMinSlope)
	droplet.Capacity = slopeTerm * droplet.Vel * droplet.Water * em.cfg.PCapacity

	excess := false
	toDeposit := 0.0
	toErode := 0.0

	// If the droplet carries more sediment than its capacity, drop a percentage of the sediment based on PDeposition
	if droplet.Sediment > droplet.Capacity {
		toDeposit = (droplet.Sediment - droplet.Capacity) * em.cfg.PDeposition
		droplet.Sediment -= toDeposit
		excess = true
	}

	// Add the eroded sediment to the droplet
	if droplet.Capacity > droplet.Sediment {
		// increase sediment by: min((c −sediment) · perosion,−hdiff)
		// Note: We're not performing this step on a per-site basis here.
		toErode = (droplet.Capacity - droplet.Sediment) * em.cfg.PErosion
		droplet.Sediment += toErode
		excess = false
	}

	// Deposit if excess sediment, otherwise erode
	if excess {

		for i := 0; i < len(siteIndices); i++ {
			em.ErosionHeightDeltas[chunk_coord][siteIndices[i]] += w[i] * toDeposit
		}

	} else {
		
		for i := 0; i < len(siteIndices); i++ {
			em.ErosionHeightDeltas[chunk_coord][siteIndices[i]] -= w[i] * toErode
		}

	}
}

func (em *ErosionManager) UphillDeposit(
	chunk_coord core.ChunkCoord, // So we can index into chunks -> height deltas. Alternatively we could pass in the current ChunkCoord being processed from the beginning
	srcIndices []core.SiteIndex, 
	dstIndices []core.SiteIndex,
	w []float64,
	dstWeights []float64,
	dstHeights []float64,
	droplet *Droplet, 
	hdiff float64,
	) {

	// TODO: This is dimensional weirdness...
	if hdiff > droplet.Sediment {
		// Drop all sediment at the source site
		for i := 0; i < len(srcIndices); i++ {
			// Note: We could also use the dst weights here, experiment w/ that as well
			em.ErosionHeightDeltas[chunk_coord][srcIndices[i]] += droplet.Sediment * w[i]
		}

		em.KillDroplet(droplet)
		return

	} else {

		// Fill as much as we can
		deposit := math.Min(droplet.Sediment, hdiff)    

		for i := 0; i < len(srcIndices); i++ {
			em.ErosionHeightDeltas[chunk_coord][srcIndices[i]] += w[i] * deposit
		}

		droplet.Sediment -= deposit
		if droplet.Sediment < 1e-6 {
			// This produces little alluvial fans / puddle deposits, in theory
			droplet.Sediment = 0.0
			droplet.Vel = 0.0
			droplet.Dir = core.Vec2{X: 0.0, Y: 0.0}
		}
	}

}

func (em *ErosionManager) UpdateVelocity(droplet *Droplet, hdiff float64) {
	droplet.Vel = math.Sqrt(droplet.Vel * droplet.Vel + hdiff * em.cfg.Gravity)
}

func (em *ErosionManager) Evaporate(droplet *Droplet) {
	droplet.Water *= (1 - em.cfg.PEvaporation)
	if droplet.Water < 1e-6 {
		em.KillDroplet(droplet)
	}
}

func (em *ErosionManager) NewDroplet(x, z int, rng *rand.Rand) Droplet {

	pos := core.Vec2{
		X: float64(x),
		Y: float64(z),
	}

	// Initialize a random direction unit vector
	dir := core.Vec2{
		X: rng.Float64()*2 - 1,
		Y: rng.Float64()*2 - 1,
	}
	dir = dir.Normalize()

	return Droplet{
		Pos: pos,
		Capacity: em.cfg.PCapacity,
		Vel: 1.0,
		Dir: dir,
		Water: 0.8,
		Sediment: 0.0,
		IsAlive: true,
	}
}

func (em *ErosionManager) KillDroplet(droplet *Droplet) {
	droplet.IsAlive = false
}

func (em *ErosionManager) ClearDroplets() {
	em.Droplets = em.Droplets[:0]
}

func (em *ErosionManager) ClearHeightDeltas() {
	em.ErosionHeightDeltas = make(map[core.ChunkCoord]map[core.SiteIndex]float64)
}