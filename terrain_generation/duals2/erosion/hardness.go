package eros

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// HardnessConfig controls how rock hardness is computed across the terrain.
// Hardness values range from 0.0 (soft, easily eroded) to 1.0 (hard, resists erosion).
type HardnessConfig struct {
	Enabled         bool    // Whether hardness affects erosion
	ElevationWeight float64 // How much elevation affects hardness (0-1)
	NoiseWeight     float64 // How much noise variation to add (0-1)
	NoiseFrequency  float64 // Frequency of hardness noise (lower = larger features)
	BaseHardness    float64 // Minimum hardness for all sites (0-1)
	ElevationPower  float64 // Exponent for elevation curve (1=linear, 2=quadratic)
}

// DefaultHardnessConfig returns a balanced hardness configuration.
func DefaultHardnessConfig() HardnessConfig {
	return HardnessConfig{
		Enabled:         true,
		ElevationWeight: 0.7,  // Elevation is primary factor
		NoiseWeight:     0.3,  // Some random variation
		NoiseFrequency:  0.01, // Large-scale variation
		BaseHardness:    0.1,  // Even valleys have some resistance
		ElevationPower:  2.0,  // Quadratic - peaks are much harder
	}
}

// SubtleHardness returns a configuration with minimal hardness variation.
func SubtleHardness() HardnessConfig {
	return HardnessConfig{
		Enabled:         true,
		ElevationWeight: 0.5,
		NoiseWeight:     0.2,
		NoiseFrequency:  0.02,
		BaseHardness:    0.2,
		ElevationPower:  1.5,
	}
}

// StrongHardness returns a configuration where peaks strongly resist erosion.
func StrongHardness() HardnessConfig {
	return HardnessConfig{
		Enabled:         true,
		ElevationWeight: 0.8,
		NoiseWeight:     0.4,
		NoiseFrequency:  0.005,
		BaseHardness:    0.05,
		ElevationPower:  3.0, // Cubic - very strong peak protection
	}
}

// SimplexNoise2D is a reference to the noise function (injected to avoid import cycle)
var SimplexNoise2D func(x, z float64) float64

// ComputeHardnessMap calculates hardness values for all sites in a chunk.
// Returns a slice where hardness[siteIndex] is the hardness value (0-1).
func ComputeHardnessMap(chunk *core.TerrainChunk, cfg HardnessConfig) []float64 {
	numSites := len(chunk.Heights)
	hardness := make([]float64, numSites)

	if !cfg.Enabled {
		// Return zero hardness (no resistance) when disabled
		return hardness
	}

	// Find min/max heights for normalization
	minHeight, maxHeight := chunk.Heights[0], chunk.Heights[0]
	for _, h := range chunk.Heights {
		if h < minHeight {
			minHeight = h
		}
		if h > maxHeight {
			maxHeight = h
		}
	}

	heightRange := maxHeight - minHeight
	if heightRange < 1e-9 {
		// Flat terrain - use base hardness everywhere
		for i := range hardness {
			hardness[i] = cfg.BaseHardness
		}
		return hardness
	}

	// Compute hardness for each site
	for i := 0; i < numSites; i++ {
		site := chunk.Mesh.Sites[i]
		height := chunk.Heights[i]

		// Normalize elevation to 0-1
		normalizedHeight := (height - minHeight) / heightRange

		// Apply power curve (makes peaks much harder than mid-elevations)
		elevationFactor := math.Pow(normalizedHeight, cfg.ElevationPower)

		// Sample noise for variation (if noise function is available)
		noiseFactor := 0.0
		if SimplexNoise2D != nil {
			// Noise returns [-1, 1], normalize to [0, 1]
			noiseValue := SimplexNoise2D(site.Pos.X*cfg.NoiseFrequency, site.Pos.Y*cfg.NoiseFrequency)
			noiseFactor = (noiseValue + 1.0) / 2.0
		}

		// Combine factors with weights
		weightedHardness := elevationFactor*cfg.ElevationWeight + noiseFactor*cfg.NoiseWeight

		// Normalize by total weight and add base hardness
		totalWeight := cfg.ElevationWeight + cfg.NoiseWeight
		if totalWeight > 0 {
			weightedHardness /= totalWeight
		}

		// Final hardness = base + weighted contribution
		hardness[i] = cfg.BaseHardness + weightedHardness*(1.0-cfg.BaseHardness)

		// Clamp to [0, 1]
		if hardness[i] < 0 {
			hardness[i] = 0
		}
		if hardness[i] > 1 {
			hardness[i] = 1
		}
	}

	return hardness
}

// GetHardnessAt returns the hardness value at a specific site.
// If hardnessMap is nil, returns 0 (no resistance).
func GetHardnessAt(hardnessMap []float64, site core.SiteIndex) float64 {
	if hardnessMap == nil || int(site) >= len(hardnessMap) {
		return 0
	}
	return hardnessMap[site]
}

// AverageHardness returns the average hardness of multiple sites.
// Useful when erosion affects multiple vertices with barycentric weights.
func AverageHardness(hardnessMap []float64, sites []core.SiteIndex, weights []float64) float64 {
	if hardnessMap == nil || len(sites) == 0 {
		return 0
	}

	if len(weights) != len(sites) {
		// Uniform weights if not provided correctly
		total := 0.0
		for _, s := range sites {
			if int(s) < len(hardnessMap) {
				total += hardnessMap[s]
			}
		}
		return total / float64(len(sites))
	}

	// Weighted average
	total := 0.0
	weightSum := 0.0
	for i, s := range sites {
		if int(s) < len(hardnessMap) {
			total += hardnessMap[s] * weights[i]
			weightSum += weights[i]
		}
	}

	if weightSum < 1e-9 {
		return 0
	}
	return total / weightSum
}
