package hydro

// combineSeeds creates a deterministic seed from multiple values.
func combineSeeds(a, b, c int64) int64 {
	// Simple hash combination
	h := a
	h = h*31 + b
	h = h*31 + c
	return h
}