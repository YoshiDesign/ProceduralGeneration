package hexa

type Tile struct {
	biome         int
	slope         float32
	slopeGradient Vec3
}

type TileSet struct {
	tiles []Tile
	count int
}