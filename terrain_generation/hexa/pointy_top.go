package hexa

func MakeChunkTiles() TileSet {

	tiles_ := []Tile{


	}

	return TileSet{
		tiles: tiles_,
		count: len(tiles_),
	}

}