package squares

import "procedural_generation/terrain_generation/hexa"

type Tile struct {
	biome         int
	slope         float32
	slopeGradient hexa.Vec3
}

type TileSet struct {
	tiles []Tile
	count int
}

func MakeChunkTiles() TileSet {

	tiles_ := []Tile{}

	return TileSet{
		tiles: tiles_,
		count: len(tiles_),
	}

}