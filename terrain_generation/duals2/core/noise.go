package core

// NoiseParams holds configurable parameters for fractal noise terrain generation.
type NoiseParams struct {
	Octaves     int
	Frequency   float64
	Amplitude   float64
	Persistence float64
	Lacunarity  float64
}
