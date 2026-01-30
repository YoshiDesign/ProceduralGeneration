package eros

/** Some Notes:

### Direction Calculation
- In [Vol] a random direction is chosen if the new directions magnitude is below a
threshold. This happens most often when a drop is initialized inside a flat cell. My
algorithm only chooses a random direction if the direction otherwise would be 0. The
direction of motion is then normalized.

This direction is a blended value between g and dirold defined by a parameter
pinertia with a value between 0 and 1. 1 means g is not taken into account and the
direction never changes and 0 means the previous direction dirold is ignored and the
new direction is along the negative gradient.

dirNew = dirOld · pinertia − g · (1 − pinertia)
posNew += dirNew


*/

import "procedural_generation/terrain_generation/duals2/core"

type ErosionManager struct {
	Config ErosionConfig
}

type ErosionConfig struct {
	Pinertia float64
	MinSlopeThreshold float64
	Gravity float64
}

type Droplet struct {
	Pos core.ChunkCoord
	Vel core.Vec2
	Dir core.Vec2 // According to gradient
	Water float64
	Sediment float64
}

func MoveDroplet(droplet *Droplet, dt float64) {

}

func Erode(site *core.Site, droplet *Droplet, idx core.SiteIndex, heights []float64) {
	//siteHeight := heights[idx]
}

func Deposit(site *core.Site, droplet *Droplet, idx core.SiteIndex, heights []float64) {
	//siteHeight := heights[idx]
}