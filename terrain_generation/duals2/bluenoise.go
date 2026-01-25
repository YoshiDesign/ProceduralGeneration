package duals2

import (
	"math"
	"math/rand"
	"procedural_generation/terrain_generation/duals2/core"
)

// BlueNoiseConfig configures Poisson disk sampling.
type BlueNoiseConfig struct {
	MinDist float64 // Minimum distance between points
	MaxTries int    // Attempts per active point before giving up (typically 30)
}

// DefaultBlueNoiseConfig returns sensible defaults for terrain generation.
func DefaultBlueNoiseConfig(minDist float64) BlueNoiseConfig {
	return BlueNoiseConfig{
		MinDist:  minDist,
		MaxTries: 30,
	}
}

// GenerateBlueNoise generates Poisson disk–distributed points within a rectangular region
// using Bridson's algorithm. The region spans [minX, maxX) × [minZ, maxZ).
//
// rng is a seeded random source for determinism.
// Returns a slice of core.Vec2 points.
func GenerateBlueNoise(rng *rand.Rand, minX, minZ, maxX, maxZ float64, cfg BlueNoiseConfig) []core.Vec2 {
	if cfg.MinDist <= 0 {
		return nil
	}
	if cfg.MaxTries <= 0 {
		cfg.MaxTries = 30
	}

	width := maxX - minX
	height := maxZ - minZ
	if width <= 0 || height <= 0 {
		return nil
	}

	// fmt.Printf("GenerateBlueNoise---\nmax: (%v, %v)\nmin: (%v, %v)\nwidth: %v\theight: %v ----", maxX, maxZ, minX, minZ, width, height)

	// Cell size for the background grid: r / sqrt(2) guarantees at most one point per cell
	cellSize := cfg.MinDist / math.Sqrt(2)
	gridW := int(math.Ceil(width / cellSize))
	gridH := int(math.Ceil(height / cellSize))

	// Grid stores index into points slice, or -1 if empty
	grid := make([]int, gridW*gridH)
	for i := range grid {
		grid[i] = -1
	}

	points := make([]core.Vec2, 0, gridW*gridH/4)
	active := make([]int, 0, 128)

	// Helper: convert world coords to a local grid cell
	toGrid := func(p core.Vec2) (int, int) {
		gx := int((p.X - minX) / cellSize) 
		gz := int((p.Y - minZ) / cellSize)
		// Clamp to valid range
		if gx < 0 {
			gx = 0
		}
		if gx >= gridW {
			gx = gridW - 1
		}
		if gz < 0 {
			gz = 0
		}
		if gz >= gridH {
			gz = gridH - 1
		}
		return gx, gz
	}

	// Helper: check if point is valid (in bounds and far enough from neighbors)
	isValid := func(p core.Vec2) bool {
		if p.X < minX || p.X >= maxX || p.Y < minZ || p.Y >= maxZ {
			return false
		}
		gx, gz := toGrid(p)

		// Check surrounding cells (5x5 neighborhood is sufficient for r/sqrt(2) cell size)
		r2 := cfg.MinDist * cfg.MinDist
		for dz := -2; dz <= 2; dz++ {
			for dx := -2; dx <= 2; dx++ {
				nx, nz := gx+dx, gz+dz
				if nx < 0 || nx >= gridW || nz < 0 || nz >= gridH {
					continue
				}
				idx := grid[nz*gridW+nx]
				if idx != -1 {
					diff := points[idx].Sub(p)
					if diff.Len2() < r2 {
						return false
					}
				}
			}
		}
		return true
	}

	// Insert a point
	insert := func(p core.Vec2) {
		idx := len(points)
		points = append(points, p)
		active = append(active, idx)
		gx, gz := toGrid(p) // The grid cell this point belongs to
		grid[gz*gridW+gx] = idx // Store the index of this point in the 1D array (grid)
	}

	// Start with a random initial point
	startX := minX + rng.Float64()*width // somewhere between min and a fraction of max
	startZ := minZ + rng.Float64()*height
	insert(core.Vec2{X: startX, Y: startZ})

	// Main loop
	for len(active) > 0 {
		// Pick a random active point
		ai := rng.Intn(len(active)) // Always 0 on first iteration
		pi := active[ai]
		p := points[pi]

		found := false
		for k := 0; k < cfg.MaxTries; k++ {
			// Generate a random point in the annulus [r, 2r] around p
			angle := rng.Float64() * 2 * math.Pi
			dist := cfg.MinDist + rng.Float64()*cfg.MinDist // min dist plus a fraction of the min dist
			candidate := core.Vec2{
				X: p.X + dist*math.Cos(angle),
				Y: p.Y + dist*math.Sin(angle),
			}

			if isValid(candidate) {
				insert(candidate)
				found = true
				break
			}
		}

		if !found {
			// Remove from active list (swap with last, then pop)
			active[ai] = active[len(active)-1]
			active = active[:len(active)-1]
		}
	}

	return points
}

// GenerateBlueNoiseSeeded is a convenience wrapper that creates a seeded RNG.
func GenerateBlueNoiseSeeded(seed int64, minX, minZ, maxX, maxZ float64, cfg BlueNoiseConfig) []core.Vec2 {
	rng := rand.New(rand.NewSource(seed))
	return GenerateBlueNoise(rng, minX, minZ, maxX, maxZ, cfg)
}
