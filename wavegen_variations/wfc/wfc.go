package wfc

import (
	"fmt"
	"image/color"
	"math"
	"math/bits"
	"math/rand"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

const (
	GridW      = 40
	GridH      = 24
	CellSizePx = 24
)

var CollapsedDomains []int

type Solver struct {
	W, H int

	TS  TileSet
	Rng *rand.Rand

	Domains []uint64 // per cell domain bitset
	Collapsed []int
	Steps   int

	contradiction bool
	done          bool
}

func NewSolver(w, h int, ts TileSet, rng *rand.Rand) *Solver {
	return &Solver{
		W:       w, // 40
		H:       h, // 24
		TS:      ts,
		Rng:     rng,
		// Collapsed: make([]int, w*h),
		Collapsed: make([]int, 0, w*h),
		Domains: make([]uint64, w*h),
	}
}

func (s *Solver) Reset() {
	for i := range s.Domains {
		s.Domains[i] = s.TS.AllMask
	}
	s.Steps = 0
	s.contradiction = false
	s.done = false
}

func (s *Solver) Done() bool           { return s.done }
func (s *Solver) Contradiction() bool  { return s.contradiction }
func (s *Solver) idx(x, y int) int     { return y*s.W + x }
func (s *Solver) DomainAt(x, y int) uint64 { return s.Domains[s.idx(x, y)] }

func (s *Solver) CollapsedCount() int {
	n := 0
	for _, d := range s.Domains {
		if bits.OnesCount64(d) == 1 {
			n++
		}
	}
	return n
}

func (s *Solver) Step() {
	if s.done || s.contradiction {
		return
	}
	s.Steps++

	// Get cell with least entropy
	// Determine which tile it should be based on its domain
	// Alter the domains of neighboring tiles

	// 1) get the index of the cell (Domains[c]) with the lowest entropy (>1 possibilities)
	c := s.pickLowestEntropyCell() 
	if c < 0 {
		s.done = true
		return
	}

	// 2) collapse it (weighted random among its options)
	d := s.Domains[c]
	choice := s.weightedPick(d) // Pick which tile it should be.
	if choice < 0 {
		s.contradiction = true
		return
	}
	// Domain for cell `c` now consists of only 1 set bit, the chosen tile represented by its index in Tiles[]
	s.Domains[c] = (uint64(1) << uint(choice))
	s.Collapsed = append(s.Collapsed, c)

	// 3) propagate constraints outward - altering neighbor domains
	s.propagateFrom(c)
}

func (s *Solver) ForceRandomTile(x, y int) {
	if s.done || s.contradiction {
		return
	}
	i := s.idx(x, y)
	d := s.Domains[i]
	if d == 0 {
		s.contradiction = true
		return
	}
	t := s.weightedPick(d)
	if t < 0 {
		s.contradiction = true
		return
	}
	s.Domains[i] = (uint64(1) << uint(t))
	s.propagateFrom(i)
}

func (s *Solver) pickLowestEntropyCell() int {
	best := -1
	bestCount := math.MaxInt32

	// scan-based for skeleton simplicity; swap to heap later
	for i, d := range s.Domains {
		c := bits.OnesCount64(d) // number of possibilities
		if c == 0 {
			s.contradiction = true
			return -1
		}
		if c == 1 {
			continue
		}
		if c < bestCount {
			bestCount = c
			best = i // take the cell's index
		} else if c == bestCount && best >= 0 {
			// random tie-break to avoid patterns
			if s.Rng.Intn(2) == 0 {
				best = i
			}
		}
	}
	return best // index of CELL with lowest entropy (least number of 1's in its domain)
}

/* 	pick the index based upon the bits of the (already acquired) lowest entropy 
	domain, representing the tile we're going to select */
func (s *Solver) weightedPick(domain uint64) int {
	// For prototype: iterate set bits and do a weighted roulette.
	// (You can replace with alias tables per domain, etc.)
	var total float64
	type opt struct {
		tile int
		w    float64
	}
	opts := make([]opt, 0, bits.OnesCount64(domain))

	m := domain
	for m != 0 { // This implies that the domain has at least one valid state
		b := bits.TrailingZeros64(m) // the index of a possible tile `t` from the s.TS.Tiles array
		t := b
		m &^= (uint64(1) << uint(b)) // Clear the domain bit at `b`
		w := float64(s.TS.Tiles[t].Weight) // Read the weight assigned to this tile
		if w <= 0 {
			continue	// ignore
		}
		total += w
		opts = append(opts, opt{tile: t, w: w})
	}
	if total <= 0 || len(opts) == 0 {
		return -1
	}

	/* Super rudimentary weighted selection but still kinda sensible as long as weights are drastically different */
	r := s.Rng.Float64() * total // reduce the total based on a random coefficient from [0.0, 1.0)
	for _, o := range opts {
		r -= o.w
		if r <= 0 {
			return o.tile
		}
	}
	return opts[len(opts)-1].tile
}

/* Beginning from cell at Domain[start] */ 
func (s *Solver) propagateFrom(start int) {
	// Fast queue (slice head index)
	q := make([]int, 0, 64)
	q = append(q, start) // First index is our starting point
	head := 0

	for head < len(q) {
		c := q[head] // Retrive the next cell from front to back
		head++

		/*
		Pro Tip:
			To get zero-based index of cell at (x, y) in a 1d array:
				y * width + x = idx

			To get the (x, y) of a cell based on the cells index in the 1D array:
				x = idx % width
				y = idx / width
		
		*/

		// Get the x, y of C in cell coordinates based on its flat index
		x := c % s.W 
		y := c / s.W 

		curDomain := s.Domains[c]
		if curDomain == 0 { // The domain represents no tiles. Oh no
			s.contradiction = true
			return
		}

		/* This loop implies that 4 is a hard coded constant of the number of borders that dictate where a neighbor is */
		// For each neighbor, constrain neighbor domain based on current domain.
		for dir := Dir(0); dir < 4; dir++ {
			nx, ny := x, y
			switch dir {
			case N:
				ny--
			case E:
				nx++
			case S:
				ny++
			case W:
				nx--
			}
			if nx < 0 || nx >= s.W || ny < 0 || ny >= s.H {
				continue // Not in our area of interest
			}
			// The index of this neighbor
			ni := s.idx(nx, ny)
			// cache it
			oldN := s.Domains[ni]
			if oldN == 0 { // This neighbor's domain has a mask with no valid possibilities.
				s.contradiction = true
				return
			}

			// Using our current domain, find the mask of compatible domains for neighbor at `dir`
			allowed := s.allowedByNeighborRule(curDomain, dir)
			// & the mask to the neighbor's domain to effectively collapse the neighbor to its possible choices
			newN := oldN & allowed

			if newN != oldN {
				// Update the neighbor's domain
				s.Domains[ni] = newN
				// Neighbor can't have 0 possible states
				if newN == 0 {
					s.contradiction = true
					return
				}
				// only append if the `Domains[ni]` was shrunk
				q = append(q, ni)
			}
		}
	}
}

/* Iterate over a `domain` and construct a mask that represents every compatible index from a given direction */
func (s *Solver) allowedByNeighborRule(domain uint64, dir Dir) uint64 {
	// allowed = union_{t in domain} compat[dir][t]
	var allowed uint64
	m := domain
	for m != 0 {	// Acquire every compatible index
		b := bits.TrailingZeros64(m) 	// The first compatible Tiles[] index based on our domain arg
		m &^= (uint64(1) << uint(b))	// set << b to zero - so we can continue to iterate over the domain
		allowed |= s.TS.Compat[dir][b] 	// There's one `b` for every tile in the tile set
	}
	return allowed // Returns a mask representing the compatible indices from Tiles[]
}

func (s *Solver) CellColor(domain uint64) color.Color {
	c := bits.OnesCount64(domain)
	if c == 0 {
		return color.RGBA{255, 0, 0, 255} // contradiction marker
	}
	if c == 1 {
		t := bits.TrailingZeros64(domain)
		return s.TS.Tiles[t].Color
	}
	// Uncollapsed: visualize entropy (more options = brighter)
	// Map [2..tileCount] -> [40..200]
	tcount := s.TS.TileCount
	v := 40 + int(float64(c-2)/float64(max(1, tcount-2))*160.0)
	if v > 220 {
		v = 220
	}
	return color.RGBA{uint8(v), uint8(v), uint8(v), 255}
}

func (s *Solver) RoadColor(x float32, y float32) color.Color {

	return color.RGBA{0, 0, 255, 255}
}

func (s *Solver) RoadLines(screen *ebiten.Image) {

	var screenXr1 float32
	var screenXr2 float32
	var screenYr1 float32
	var screenYr2 float32

	var height1 float32
	var height2 float32
	var width1 float32
	var width2 float32

	for _, val := range s.Collapsed {

		if bits.OnesCount64(s.Domains[val]) > 1 {
			fmt.Println("Skipping:\t", s.Domains[val])
			return
		}

		cMask := s.Domains[val]

		t := bits.TrailingZeros(uint(cMask)) // Type of tile
		x := float32((val % s.W) * CellSizePx)	// screen x
		y := float32((val / s.W) * CellSizePx)	// screen y

		//tile_x := (val % s.W)
		//tile_y := (val / s.W)

		midpoint := float32(CellSizePx / 2)

		//NS := false
		//EW := false
		// S := false
		// E := false
		// W := false

		switch t {
		case 0:
			continue
		case 1: // NS
			
			height1 = midpoint
			height2 = midpoint
			width1 = 3
			width2 = 3

			// N
			screenXr1 = (x + midpoint - 2)
			screenYr1 = y
			// S
			screenXr2 = (x + midpoint - 2)
			screenYr2 = (y + midpoint)

		case 2: // EW

			height1 = 3
			height2 = 3 	
			width1 = midpoint
			width2 = midpoint

			// E
			screenXr1 = midpoint + x
			screenYr1 = (y + midpoint - 1)

			// W
			screenXr2 = x
			screenYr2 = (y + midpoint - 1)

		case 3: // Corner_NE

			height1 = midpoint
			height2 = 3
			width1 = 3
			width2 = midpoint

			// N
			screenXr1 = (x + midpoint - 2)
			screenYr1 = y

			// E
			screenXr2 = midpoint + x
			screenYr2 = (y + midpoint - 1)

		case 4: // Corner_ES

			height1 = 3
			height2 = midpoint
			width1 = midpoint
			width2 = 3

			// E
			screenXr1 = midpoint + x
			screenYr1 = (y + midpoint - 1)
			// S
			screenXr2 = (x + midpoint - 2)
			screenYr2 = (y + midpoint)
		case 5: // Corner_SW

			height1 = midpoint
			height2 = 3
			width1 = 3
			width2 = midpoint

			// S
			screenXr1 = (x + midpoint - 2)
			screenYr1 = (y + midpoint)
			// W
			screenXr2 = x
			screenYr2 = (y + midpoint - 1)
		case 6: // Corner_WN

			// W
			screenXr2 = x
			screenYr2 = (y + midpoint - 1)
			height2 = 3
			width2 = midpoint
		
			// N
			screenXr1 = (x + midpoint - 2)
			screenYr1 = y
			height1 = midpoint
			width1 = 3
			
		}


		// r1
		vector.FillRect(screen,
			screenXr1,
			screenYr1,
			width1,
			height1,
			s.RoadColor(x, y),
			false,
		)

		// r2
		vector.FillRect(screen,
			screenXr2,
			screenYr2,
			width2,
			height2,
			s.RoadColor(x, y),
			false,
		)
		
	}

}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
