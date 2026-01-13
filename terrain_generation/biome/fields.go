package biome

/* Meant to influence terrain characteristics */
type GlobalBiomeInput struct {
	hills       float32
	mountains   float32
	rivers      float32
	coastline   float32
	volcano     float32
}

type MacroBiomeInput struct {
	temperature float32
	elevation	float32
}