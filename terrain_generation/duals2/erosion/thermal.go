package eros

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// ThermalConfig holds parameters for thermal erosion simulation.
type ThermalConfig struct {
	Enabled       bool    // Whether thermal erosion is active
	TalusThreshold float64 // Slope threshold as tan(angle). 0.57 ≈ 30°, 1.0 ≈ 45°
	TransferRate  float64 // Fraction of excess material to transfer per iteration (0.1-0.5)
	Iterations    int     // Number of thermal erosion passes
	Timing        string  // "before", "after", or "interleaved"
}

// DefaultThermalConfig returns a balanced thermal erosion configuration.
func DefaultThermalConfig() ThermalConfig {
	return ThermalConfig{
		Enabled:       true,
		TalusThreshold: 0.57, // ~30 degrees angle of repose
		TransferRate:  0.25,  // Transfer 25% of excess per iteration
		Iterations:    20,
		Timing:        "after",
	}
}

// SubtleThermal returns a light thermal erosion configuration.
func SubtleThermal() ThermalConfig {
	return ThermalConfig{
		Enabled:       true,
		TalusThreshold: 0.7,  // ~35 degrees - only very steep slopes
		TransferRate:  0.1,   // Gentle transfer
		Iterations:    10,
		Timing:        "after",
	}
}

// HeavyThermal returns an aggressive thermal erosion configuration.
func HeavyThermal() ThermalConfig {
	return ThermalConfig{
		Enabled:       true,
		TalusThreshold: 0.4,  // ~22 degrees - triggers on gentler slopes
		TransferRate:  0.5,   // Aggressive transfer
		Iterations:    50,
		Timing:        "after",
	}
}

// ThermalErosion simulates talus/scree accumulation by moving material
// from steep slopes to neighboring lower areas.
//
// The algorithm:
// 1. For each site, compute slope to each neighbor
// 2. If slope exceeds talus threshold, transfer material proportionally
// 3. Material flows from higher site to lower neighbor
// 4. Repeat for configured number of iterations
func (em *ErosionManager) ThermalErosion(chunk *core.TerrainChunk, cfg ThermalConfig) {
	if !cfg.Enabled || cfg.Iterations <= 0 {
		return
	}

	numSites := len(chunk.Heights)
	
	// We'll accumulate deltas and apply them after each full pass
	// to avoid order-dependent artifacts
	deltas := make([]float64, numSites)

	for iter := 0; iter < cfg.Iterations; iter++ {
		// Clear deltas for this iteration
		for i := range deltas {
			deltas[i] = 0
		}

		// Process each site
		for siteIdx := 0; siteIdx < numSites; siteIdx++ {
			site := core.SiteIndex(siteIdx)
			sitePos := chunk.Mesh.Sites[site].Pos
			siteHeight := chunk.Heights[site]

			neighbors := core.GetAllNeighbors(chunk.Mesh, site)
			if len(neighbors) == 0 {
				continue
			}

			for _, neighbor := range neighbors {
				neighborPos := chunk.Mesh.Sites[neighbor].Pos
				neighborHeight := chunk.Heights[neighbor]

				// Only transfer downhill
				if neighborHeight >= siteHeight {
					continue
				}

				// Calculate horizontal distance
				dx := neighborPos.X - sitePos.X
				dz := neighborPos.Y - sitePos.Y // Vec2.Y is Z in world space
				distance := math.Sqrt(dx*dx + dz*dz)
				if distance < 1e-9 {
					continue
				}

				// Calculate slope (rise/run = tan of angle)
				heightDiff := siteHeight - neighborHeight
				slope := heightDiff / distance

				// Only transfer if slope exceeds threshold
				if slope <= cfg.TalusThreshold {
					continue
				}

				// Calculate transfer amount
				// Transfer proportional to how much the slope exceeds the threshold
				// excessSlope := slope - cfg.TalusThreshold
				
				// The maximum we could transfer is half the height difference
				// (to reach equilibrium at the threshold slope)
				maxTransfer := (heightDiff - cfg.TalusThreshold*distance) / 2.0
				
				// Apply transfer rate
				transfer := maxTransfer * cfg.TransferRate
				if transfer < 0 {
					transfer = 0
				}

				// Accumulate deltas (will apply after full pass)
				deltas[siteIdx] -= transfer
				deltas[neighbor] += transfer
			}
		}

		// Apply deltas after full iteration to avoid order-dependent results
		for i := range chunk.Heights {
			chunk.Heights[i] += deltas[i]
			chunk.Mesh.Sites[i].Height += deltas[i]
		}
	}
}

// ThermalErosionInterleaved runs a single iteration of thermal erosion.
// Useful for interleaving with hydraulic erosion steps.
func (em *ErosionManager) ThermalErosionInterleaved(chunk *core.TerrainChunk, cfg ThermalConfig) {
	if !cfg.Enabled {
		return
	}

	// Run single iteration by temporarily setting iterations to 1
	singleCfg := cfg
	singleCfg.Iterations = 1
	em.ThermalErosion(chunk, singleCfg)
}
