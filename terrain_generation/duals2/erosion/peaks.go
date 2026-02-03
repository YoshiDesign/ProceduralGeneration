package eros

// PeakConfig combines hardness and ridge enhancement settings.
// This provides a unified configuration for peak preservation and enhancement.
type PeakConfig struct {
	Hardness HardnessConfig
	Ridge    RidgeConfig
}

// DefaultPeakConfig returns a balanced peak enhancement configuration.
// Provides natural-looking peaks with moderate protection and enhancement.
func DefaultPeakConfig() PeakConfig {
	return PeakConfig{
		Hardness: DefaultHardnessConfig(),
		Ridge:    DefaultRidgeConfig(),
	}
}

// SubtlePeaks returns a configuration for gentle peak preservation.
// Good for rolling hills or terrain where peaks shouldn't be too prominent.
func SubtlePeaks() PeakConfig {
	return PeakConfig{
		Hardness: SubtleHardness(),
		Ridge:    SubtleRidges(),
	}
}

// DramaticPeaks returns a configuration for prominent, jagged peaks.
// Creates mountain ranges with sharp ridges and significant elevation contrast.
func DramaticPeaks() PeakConfig {
	return PeakConfig{
		Hardness: StrongHardness(),
		Ridge:    DramaticRidges(),
	}
}

// HardnessOnly returns a configuration with hardness but no ridge enhancement.
// Peaks are protected during erosion but not artificially enhanced afterward.
func HardnessOnly() PeakConfig {
	return PeakConfig{
		Hardness: DefaultHardnessConfig(),
		Ridge: RidgeConfig{
			Enabled: false,
		},
	}
}

// RidgesOnly returns a configuration with ridge enhancement but no hardness.
// Erosion proceeds normally, then ridges are enhanced afterward.
func RidgesOnly() PeakConfig {
	return PeakConfig{
		Hardness: HardnessConfig{
			Enabled: false,
		},
		Ridge: DefaultRidgeConfig(),
	}
}

// NoPeaks returns a configuration with both systems disabled.
// Use this for flat terrain or when peak enhancement is not desired.
func NoPeaks() PeakConfig {
	return PeakConfig{
		Hardness: HardnessConfig{Enabled: false},
		Ridge:    RidgeConfig{Enabled: false},
	}
}

// AlpinePeaks returns a configuration mimicking alpine mountain ranges.
// Strong hardness at elevation with dramatic, jagged ridges.
func AlpinePeaks() PeakConfig {
	return PeakConfig{
		Hardness: HardnessConfig{
			Enabled:         true,
			ElevationWeight: 0.85,
			NoiseWeight:     0.25,
			NoiseFrequency:  0.008,
			BaseHardness:    0.05,
			ElevationPower:  2.5,
		},
		Ridge: RidgeConfig{
			Enabled:       true,
			Threshold:     0.25,
			BoostAmount:   2.5,
			NoiseAmount:   0.6,
			NoiseFreq:     0.08,
			Iterations:    2,
			MinHeight:     0.6,          // Only enhance upper 40% of terrain
			MinHeightMode: "normalized",
		},
	}
}

// RollingHills returns a configuration for gentle, rounded terrain.
// Light hardness prevents complete flattening, minimal ridge enhancement.
func RollingHills() PeakConfig {
	return PeakConfig{
		Hardness: HardnessConfig{
			Enabled:         true,
			ElevationWeight: 0.4,
			NoiseWeight:     0.3,
			NoiseFrequency:  0.02,
			BaseHardness:    0.3,
			ElevationPower:  1.2,
		},
		Ridge: RidgeConfig{
			Enabled:       true,
			Threshold:     0.5,
			BoostAmount:   0.5,
			NoiseAmount:   0.1,
			NoiseFreq:     0.02,
			Iterations:    1,
			MinHeight:     0.2, // Enhance most of the terrain gently
			MinHeightMode: "normalized",
		},
	}
}

// VolcanicPeaks returns a configuration for volcanic-style terrain.
// Very hard peaks with moderate ridge enhancement.
func VolcanicPeaks() PeakConfig {
	return PeakConfig{
		Hardness: HardnessConfig{
			Enabled:         true,
			ElevationWeight: 0.9,
			NoiseWeight:     0.2,
			NoiseFrequency:  0.015,
			BaseHardness:    0.0,
			ElevationPower:  3.0,
		},
		Ridge: RidgeConfig{
			Enabled:       true,
			Threshold:     0.35,
			BoostAmount:   2.0,
			NoiseAmount:   0.4,
			NoiseFreq:     0.06,
			Iterations:    1,
			MinHeight:     0.7, // Only enhance the very top (volcanic peaks)
			MinHeightMode: "normalized",
		},
	}
}
