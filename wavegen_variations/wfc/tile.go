package wfc

import "image/color"

type Dir uint8

const (
	N Dir = iota
	E
	S
	W
)

func Opp(d Dir) Dir { return (d + 2) & 3 }

// Tile is a *variant* (rotation/flip already baked in if you want variants).
// For early WFC prototyping, storing per-side "socket IDs" is a clean approach.
type Tile struct {
	Name    string
	Sockets [4]uint16 // N,E,S,W socket ids
	Weight  float32
	Color   color.RGBA
}

// TileSet contains precomputed compatibility masks.
type TileSet struct {
	Tiles     []Tile
	Compat    [4][]uint64 // Compat[dir][t] => bitmask of neighbor tiles allowed in dir from tile t
	AllMask   uint64
	TileCount int
}

// MakeDemoTiles: a simple "roads" like socket demo.
// Socket IDs: 0 = empty, 1 = road
// Rules: sockets must match to connect.
func MakeDemoTiles() TileSet {
	tiles := []Tile{
		// Name: city_plot_4, Sockets: [18]{1, 0, 1, 0, 0 ...}, Weight: 8, InstanceLayout: ?{}?
		/*
			If Sockets becomes nested arrays, we can ensure we have a match for each of them.
			Hmm...
		*/
		{Name: "Empty", Sockets: [4]uint16{0, 0, 0, 0}, Weight: 4, Color: color.RGBA{255, 240, 210, 255}},
		{Name: "Road_NS", Sockets: [4]uint16{1, 0, 1, 0}, Weight: 2, Color: color.RGBA{220, 220, 220, 255}},
		{Name: "Road_EW", Sockets: [4]uint16{0, 1, 0, 1}, Weight: 3, Color: color.RGBA{220, 220, 220, 255}},
		{Name: "Corner_NE", Sockets: [4]uint16{1, 1, 0, 0}, Weight: 3, Color: color.RGBA{200, 200, 200, 255}},
		{Name: "Corner_ES", Sockets: [4]uint16{0, 1, 1, 0}, Weight: 1, Color: color.RGBA{200, 200, 200, 255}},
		{Name: "Corner_SW", Sockets: [4]uint16{0, 0, 1, 1}, Weight: 1, Color: color.RGBA{200, 200, 200, 255}},
		{Name: "Corner_WN", Sockets: [4]uint16{1, 0, 0, 1}, Weight: 3, Color: color.RGBA{200, 200, 200, 255}},
	}

	ts := NewTileSetFromSockets(tiles)
	return ts
}

func NewTileSetFromSockets(tiles []Tile) TileSet {
	if len(tiles) > 64 {
		panic("this skeleton uses uint64 bitmasks; expand to []uint64 chunks if you need >64 tiles")
	}

	var ts TileSet
	ts.Tiles = tiles
	ts.TileCount = len(tiles)
	ts.AllMask = (uint64(1) << uint(ts.TileCount)) - 1 // 0111111

	// Precompute compat masks by socket matching.
	// For a tile t and direction dir, neighbor u is allowed if:
	// t.Sockets[dir] == u.Sockets[Opp(dir)]
	for dir := 0; dir < 4; dir++ {
		ts.Compat[dir] = make([]uint64, ts.TileCount)
	}

	/* Note- This is a very simple compatibility sequence just to hammer in the understanding */

	// Create a mask for this tile, denoting all of the other tiles that align with each of its directions
	for t := 0; t < ts.TileCount; t++ {
		for dir := Dir(0); dir < 4; dir++ {
			var mask uint64
			want := ts.Tiles[t].Sockets[dir]	// 1 or 0 (or whatever our convention might be)
			od := Opp(dir)						// Index to the opposing direction of a Socket from `want`
			for u := 0; u < ts.TileCount; u++ {	// Re-iterate all tiles in the set
				if ts.Tiles[u].Sockets[od] == want {	// The opposite direction matched in value at tile `u`
					mask |= (uint64(1) << uint(u))		// Set bit `u` to 1 in our mask
				}
			}
			ts.Compat[dir][t] = mask					// Store the bitmask representing the tiles compatible with `t` on its `dir` facing side
		}
	}

	return ts
}
