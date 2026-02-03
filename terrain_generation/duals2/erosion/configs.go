package eros

// SubtleErosion creates gentle, minimal terrain modification.
// Good for adding slight weathering to otherwise pristine terrain.
func SubtleErosion() ErosionConfig {
	return ErosionConfig{
		Pinertia:     0.05,   // Low inertia - droplets follow terrain closely
		PCapacity:    4.0,    // Lower capacity - less sediment transport
		PDeposition:  0.0005, // Very slow deposition
		PErosion:     0.0005, // Very slow erosion
		PEvaporation: 0.02,   // Faster evaporation - shorter droplet lifetime
		PMinSlope:    0.001,
		Gravity:      0.05,   // Gentle gravity
		NumDroplets:  2000,   // Fewer droplets
		NumSteps:     30,     // Shorter simulation
	}
}

// AverageErosion provides balanced, natural-looking erosion.
// Suitable for general terrain generation.
func AverageErosion() ErosionConfig {
	return ErosionConfig{
		Pinertia:     0.1,
		PCapacity:    10.0,
		PDeposition:  0.001,
		PErosion:     0.001,
		PEvaporation: 0.01,
		PMinSlope:    0.07,
		Gravity:      0.1,
		NumDroplets:  4000,
		NumSteps:     50,
	}
}

// HeavyErosion creates deep valleys and significant terrain carving.
// Simulates long-term geological erosion or high rainfall environments.
func HeavyErosion() ErosionConfig {
	return ErosionConfig{
		Pinertia:     0.3,   // Slightly more inertia - creates smoother channels
		PCapacity:    8.0,   // High capacity - can transport lots of sediment
		PDeposition:  0.2,  // Faster deposition in valleys
		PErosion:     0.7,  // Aggressive erosion
		PEvaporation: 0.02,  // Slow evaporation - droplets travel far
		PMinSlope:    0.01, // Erodes even flatter areas
		Gravity:      10.0,    // Strong gravity - faster downhill movement
		NumDroplets:  100_000,   // Many droplets
		NumSteps:     200,     // Long simulation
	}
}

// ArtisticErosion creates stylized, exaggerated terrain features.
// High inertia causes swirling patterns, extreme capacity creates
// dramatic sediment deposits, and inverted gravity/erosion ratios
// produce otherworldly landscapes.
func ArtisticErosion() ErosionConfig {
	return ErosionConfig{
		Pinertia:     0.7,    // Very high inertia - droplets resist turning, create swirls
		PCapacity:    50.0,   // Extreme capacity - massive sediment transport
		PDeposition:  0.01,   // Fast deposition - creates visible sediment fans
		PErosion:     0.0002, // Very slow erosion relative to deposition
		PEvaporation: 0.002,  // Very slow evaporation - droplets travel extremely far
		PMinSlope:    0.01,   // Only erodes steeper slopes
		Gravity:      0.02,   // Weak gravity - floaty, dreamlike movement
		NumDroplets:  6000,
		NumSteps:     100,    // Long simulation to see full effect
	}
}

// ============================================================================
// Combined Erosion Presets (Hydraulic + Thermal)
// ============================================================================

// CombinedErosionConfig holds both hydraulic and thermal erosion settings.
type CombinedErosionConfig struct {
	Hydraulic ErosionConfig
	Thermal   ThermalConfig
}

// NaturalErosion returns a balanced combination of hydraulic and thermal erosion.
// Produces realistic terrain with carved channels and talus deposits.
func NaturalErosion() CombinedErosionConfig {
	return CombinedErosionConfig{
		Hydraulic: AverageErosion(),
		Thermal:   DefaultThermalConfig(),
	}
}

// AggressiveErosion returns heavy hydraulic erosion with aggressive thermal weathering.
// Creates deeply carved terrain with significant talus accumulation.
func AggressiveErosion() CombinedErosionConfig {
	return CombinedErosionConfig{
		Hydraulic: HeavyErosion(),
		Thermal:   HeavyThermal(),
	}
}

// GentleWeathering returns subtle erosion for light terrain modification.
// Good for adding natural character without dramatic changes.
func GentleWeathering() CombinedErosionConfig {
	return CombinedErosionConfig{
		Hydraulic: SubtleErosion(),
		Thermal:   SubtleThermal(),
	}
}

// HydraulicOnly returns hydraulic erosion with thermal disabled.
// Use when you only want water-based erosion.
func HydraulicOnly() CombinedErosionConfig {
	return CombinedErosionConfig{
		Hydraulic: HeavyErosion(),
		Thermal: ThermalConfig{
			Enabled: false,
		},
	}
}

// ThermalOnly returns thermal erosion with minimal hydraulic erosion.
// Use when you only want slope-based material movement.
func ThermalOnly() CombinedErosionConfig {
	return CombinedErosionConfig{
		Hydraulic: SubtleErosion(),
		Thermal:   HeavyThermal(),
	}
}
