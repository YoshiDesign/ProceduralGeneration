package eros

import (
	"math"
	"procedural_generation/terrain_generation/duals2/core"
)

// RidgeConfig controls ridge detection and enhancement.
type RidgeConfig struct {
	Enabled       bool    // Whether ridge enhancement is active
	Threshold     float64 // Minimum ridgeness score to enhance (0-1)
	BoostAmount   float64 // Maximum height boost in world units
	NoiseAmount   float64 // Jaggedness noise amplitude
	NoiseFreq     float64 // Jaggedness noise frequency
	Iterations    int     // Number of enhancement passes
	MinHeight     float64 // Minimum elevation for ridge enhancement (world units)
	MinHeightMode string  // "absolute" = world height, "normalized" = 0-1 within chunk
}

// DefaultRidgeConfig returns a balanced ridge enhancement configuration.
func DefaultRidgeConfig() RidgeConfig {
	return RidgeConfig{
		Enabled:       true,
		Threshold:     0.3,  // Moderate threshold
		BoostAmount:   1.5,  // Noticeable but not extreme
		NoiseAmount:   0.3,  // Some jaggedness
		NoiseFreq:     0.05, // Medium-scale noise
		Iterations:    1,
		MinHeight:     0.5,        // Only enhance upper half of terrain
		MinHeightMode: "normalized", // Use normalized height (0-1 within chunk)
	}
}

// SubtleRidges returns a configuration for gentle ridge enhancement.
func SubtleRidges() RidgeConfig {
	return RidgeConfig{
		Enabled:       true,
		Threshold:     0.4,
		BoostAmount:   0.8,
		NoiseAmount:   0.2,
		NoiseFreq:     0.03,
		Iterations:    1,
		MinHeight:     0.3,
		MinHeightMode: "normalized",
	}
}

// DramaticRidges returns a configuration for pronounced, jagged ridges.
func DramaticRidges() RidgeConfig {
	return RidgeConfig{
		Enabled:       true,
		Threshold:     0.2,
		BoostAmount:   3.0,
		NoiseAmount:   0.8,
		NoiseFreq:     0.1,
		Iterations:    2,
		MinHeight:     0.4,
		MinHeightMode: "normalized",
	}
}

// ComputeRidgeness calculates how "ridge-like" a site is.
// Returns a value from 0 (not a ridge) to 1 (strong ridge).
//
// A ridge is characterized by:
// - Being higher than most neighbors (but not all - that's a peak)
// - Having neighbors that are lower in opposing directions
func ComputeRidgeness(mesh *core.DelaunayMesh, heights []float64, site core.SiteIndex) float64 {
	neighbors := core.GetAllNeighbors(mesh, site)
	if len(neighbors) < 3 {
		return 0
	}

	siteHeight := heights[site]
	lowerCount := 0
	higherCount := 0

	for _, neighbor := range neighbors {
		if heights[neighbor] < siteHeight {
			lowerCount++
		} else {
			higherCount++
		}
	}

	totalNeighbors := len(neighbors)

	// Pure peaks (all neighbors lower) aren't ridges
	if higherCount == 0 {
		return 0
	}

	// Pure valleys (all neighbors higher) aren't ridges
	if lowerCount == 0 {
		return 0
	}

	// Ridgeness is maximized when most neighbors are lower but some are higher
	// This creates a curve that peaks when lowerCount/total is around 0.7-0.8
	lowerRatio := float64(lowerCount) / float64(totalNeighbors)

	// Asymmetric scoring: prefer more lower neighbors (ridge-like)
	// but require at least some higher neighbors (continuation along ridge)
	if lowerRatio < 0.5 {
		return 0 // More neighbors are higher - this is a slope, not a ridge
	}

	// Score based on how ridge-like the configuration is
	// Maximum ridgeness when ~75% of neighbors are lower
	optimalRatio := 0.75
	deviation := math.Abs(lowerRatio - optimalRatio)
	ridgeness := 1.0 - (deviation / 0.25) // Linear falloff from optimal

	if ridgeness < 0 {
		ridgeness = 0
	}

	// Also consider height variance - ridges should have significant drops
	heightVariance := computeHeightVariance(heights, site, neighbors)
	varianceFactor := math.Min(1.0, heightVariance/2.0) // Normalize by expected variance

	return ridgeness * varianceFactor
}

// computeHeightVariance calculates the variance of neighbor heights relative to the site.
func computeHeightVariance(heights []float64, site core.SiteIndex, neighbors []core.SiteIndex) float64 {
	if len(neighbors) == 0 {
		return 0
	}

	siteHeight := heights[site]
	sumSqDiff := 0.0

	for _, neighbor := range neighbors {
		diff := heights[neighbor] - siteHeight
		sumSqDiff += diff * diff
	}

	return math.Sqrt(sumSqDiff / float64(len(neighbors)))
}

// EnhanceRidges detects ridge lines and enhances them for more dramatic peaks.
func EnhanceRidges(chunk *core.TerrainChunk, cfg RidgeConfig) {
	if !cfg.Enabled || cfg.Iterations <= 0 {
		return
	}

	numSites := len(chunk.Heights)

	// Compute min/max heights for normalized mode - ONLY REQUIRED FOR NORMALIZED MODE
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

	// Helper to check if a site meets minimum height requirement
	meetsMinHeight := func(siteHeight float64) bool {
		if cfg.MinHeightMode == "absolute" {
			return siteHeight >= cfg.MinHeight
		}
		// Normalized mode (default)
		if heightRange < 1e-9 {
			return true // Flat terrain, apply everywhere
		}
		normalizedHeight := (siteHeight - minHeight) / heightRange
		return normalizedHeight >= cfg.MinHeight
	}

	for iter := 0; iter < cfg.Iterations; iter++ {
		// Compute ridgeness for all sites first (avoid order-dependent issues)
		ridgeness := make([]float64, numSites)
		for i := 0; i < numSites; i++ {
			ridgeness[i] = ComputeRidgeness(chunk.Mesh, chunk.Heights, core.SiteIndex(i))
		}

		// Apply enhancement
		for i := 0; i < numSites; i++ {
			if ridgeness[i] < cfg.Threshold {
				continue
			}

			// Check minimum height requirement
			if !meetsMinHeight(chunk.Heights[i]) {
				continue
			}

			site := chunk.Mesh.Sites[i]

			// Calculate boost proportional to ridgeness
			boost := ridgeness[i] * cfg.BoostAmount

			// Add noise for jaggedness (if noise function is available)
			noise := 0.0
			if SimplexNoise2D != nil {
				// Use a different frequency/offset than hardness for variety
				noiseValue := SimplexNoise2D(
					site.Pos.X*cfg.NoiseFreq+1000.0,
					site.Pos.Y*cfg.NoiseFreq+1000.0,
				)
				noise = noiseValue * cfg.NoiseAmount
			}

			// Apply enhancement
			enhancement := boost + noise
			chunk.Heights[i] += enhancement
			chunk.Mesh.Sites[i].Height += enhancement
		}
	}
}

// FindRidgeLines identifies connected ridge sites for visualization or further processing.
// Returns groups of connected ridge sites.
func FindRidgeLines(mesh *core.DelaunayMesh, heights []float64, threshold float64) [][]core.SiteIndex {
	numSites := len(heights)
	ridgeness := make([]float64, numSites)

	// Compute ridgeness for all sites
	for i := 0; i < numSites; i++ {
		ridgeness[i] = ComputeRidgeness(mesh, heights, core.SiteIndex(i))
	}

	// Find ridge sites
	isRidge := make([]bool, numSites)
	for i := 0; i < numSites; i++ {
		isRidge[i] = ridgeness[i] >= threshold
	}

	// Group connected ridge sites using flood fill
	visited := make([]bool, numSites)
	var ridgeLines [][]core.SiteIndex

	for i := 0; i < numSites; i++ {
		if !isRidge[i] || visited[i] {
			continue
		}

		// Start a new ridge line
		var line []core.SiteIndex
		queue := []core.SiteIndex{core.SiteIndex(i)}

		for len(queue) > 0 {
			current := queue[0]
			queue = queue[1:]

			if visited[current] {
				continue
			}
			visited[current] = true

			if !isRidge[current] {
				continue
			}

			line = append(line, current)

			// Add unvisited ridge neighbors to queue
			neighbors := core.GetAllNeighbors(mesh, current)
			for _, neighbor := range neighbors {
				if !visited[neighbor] && isRidge[neighbor] {
					queue = append(queue, neighbor)
				}
			}
		}

		if len(line) > 0 {
			ridgeLines = append(ridgeLines, line)
		}
	}

	return ridgeLines
}
