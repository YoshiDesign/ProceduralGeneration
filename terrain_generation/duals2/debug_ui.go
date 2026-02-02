package duals2

import (
	"bytes"
	"fmt"
	"image/color"
	"procedural_generation/terrain_generation/duals2/core"
	"procedural_generation/terrain_generation/duals2/hydro"
	"strconv"
	"time"

	"github.com/ebitenui/ebitenui"
	"github.com/ebitenui/ebitenui/image"
	"github.com/ebitenui/ebitenui/widget"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/text/v2"
	"golang.org/x/image/font/gofont/goregular"
)

// DebugUI provides a debug panel for tuning terrain generation parameters.
type DebugUI struct {
	ui       *ebitenui.UI
	visible  bool
	fontFace text.Face

	// Current parameter values (editable copies)
	NoiseParams core.NoiseParams
	HydroConfig hydro.HydroConfig

	// Callback for regeneration
	OnRegenerate func(core.NoiseParams, hydro.HydroConfig)

	// Value display labels (for updating when sliders change)
	noiseLabels map[string]*widget.Text
	hydroLabels map[string]*widget.Text

	// Debounce state for auto-regeneration
	dirty          bool
	lastChangeTime time.Time
	debounceDelay  time.Duration
}

// NewDebugUI creates a new debug UI panel.
func NewDebugUI(initialNoise core.NoiseParams, initialHydro hydro.HydroConfig, onRegenerate func(core.NoiseParams, hydro.HydroConfig)) *DebugUI {
	d := &DebugUI{
		NoiseParams:   initialNoise,
		HydroConfig:   initialHydro,
		OnRegenerate:  onRegenerate,
		visible:       false,
		noiseLabels:   make(map[string]*widget.Text),
		hydroLabels:   make(map[string]*widget.Text),
		debounceDelay: 150 * time.Millisecond,
	}

	d.fontFace = d.loadFont()
	d.ui = d.buildUI()

	return d
}

// loadFont loads a basic font for the UI.
func (d *DebugUI) loadFont() text.Face {
	source, err := text.NewGoTextFaceSource(bytes.NewReader(goregular.TTF))
	if err != nil {
		panic(err)
	}
	return &text.GoTextFace{
		Source: source,
		Size:   14,
	}
}

// buildUI constructs the ebitenui interface.
func (d *DebugUI) buildUI() *ebitenui.UI {
	// Create a root container
	rootContainer := widget.NewContainer(
		widget.ContainerOpts.Layout(widget.NewAnchorLayout()),
	)

	// Create the main panel container (positioned top-left)
	panelContainer := widget.NewContainer(
		widget.ContainerOpts.Layout(widget.NewRowLayout(
			widget.RowLayoutOpts.Direction(widget.DirectionVertical),
			widget.RowLayoutOpts.Padding(widget.NewInsetsSimple(10)),
			widget.RowLayoutOpts.Spacing(8),
		)),
		widget.ContainerOpts.BackgroundImage(d.createPanelBackground()),
		widget.ContainerOpts.WidgetOpts(
			widget.WidgetOpts.LayoutData(widget.AnchorLayoutData{
				HorizontalPosition: widget.AnchorLayoutPositionStart,
				VerticalPosition:   widget.AnchorLayoutPositionStart,
				Padding:            widget.NewInsetsSimple(10),
			}),
			widget.WidgetOpts.MinSize(320, 0),
		),
	)

	// Title
	panelContainer.AddChild(d.createLabel("DEBUG PARAMETERS", color.RGBA{255, 220, 100, 255}))

	// Noise Parameters Section
	panelContainer.AddChild(d.createLabel("-- Noise --", color.RGBA{180, 180, 255, 255}))
	panelContainer.AddChild(d.createIntSlider("Octaves", &d.NoiseParams.Octaves, 1, 12, "octaves"))
	panelContainer.AddChild(d.createFloatSlider("Frequency", &d.NoiseParams.Frequency, 0.001, 0.05, "frequency"))
	panelContainer.AddChild(d.createFloatSlider("Amplitude", &d.NoiseParams.Amplitude, 1.0, 50.0, "amplitude"))
	panelContainer.AddChild(d.createFloatSlider("Persistence", &d.NoiseParams.Persistence, 0.1, 1.0, "persistence"))
	panelContainer.AddChild(d.createFloatSlider("Lacunarity", &d.NoiseParams.Lacunarity, 1.0, 4.0, "lacunarity"))

	// Hydrology Parameters Section
	panelContainer.AddChild(d.createLabel("-- Hydrology --", color.RGBA{180, 180, 255, 255}))
	panelContainer.AddChild(d.createFloatSlider("Sea Level", &d.HydroConfig.SeaLevel, -10.0, 10.0, "seaLevel"))
	panelContainer.AddChild(d.createFloatSlider("River Src Elev", &d.HydroConfig.RiverSourceMinElev, 0.0, 20.0, "riverSrcElev"))
	panelContainer.AddChild(d.createFloatSlider("River Max Elev", &d.HydroConfig.RiverSourceMaxElev, 0.0, 20.0, "riverMaxElev"))
	panelContainer.AddChild(d.createFloatSlider("River Prob", &d.HydroConfig.RiverSourceProbability, 0.0, 1.0, "riverProb"))
	panelContainer.AddChild(d.createFloatSlider("Lake Prob", &d.HydroConfig.LakeProbability, 0.0, 1.0, "lakeProb"))
	panelContainer.AddChild(d.createFloatSlider("Lake Min Depth", &d.HydroConfig.LakeMinDepth, 0.0, 10.0, "lakeMinDepth"))
	panelContainer.AddChild(d.createFloatSlider("Ocean Prob", &d.HydroConfig.OceanProbability, 0.0, 1.0, "oceanProb"))

	// Instructions
	panelContainer.AddChild(d.createLabel("Changes apply automatically", color.RGBA{128, 128, 128, 255}))
	panelContainer.AddChild(d.createLabel("Press D to toggle panel", color.RGBA{128, 128, 128, 255}))

	rootContainer.AddChild(panelContainer)

	return &ebitenui.UI{Container: rootContainer}
}

// createPanelBackground creates a semi-transparent background for the panel.
func (d *DebugUI) createPanelBackground() *image.NineSlice {
	// Create a simple colored background
	img := ebiten.NewImage(1, 1)
	img.Fill(color.RGBA{30, 35, 45, 230})
	return image.NewNineSliceSimple(img, 0, 0)
}

// createLabel creates a text label.
func (d *DebugUI) createLabel(text string, clr color.Color) *widget.Text {
	return widget.NewText(
		widget.TextOpts.Text(text, &d.fontFace, clr),
		widget.TextOpts.WidgetOpts(
			widget.WidgetOpts.LayoutData(widget.RowLayoutData{
				Stretch: true,
			}),
		),
	)
}

// createIntSlider creates a labeled slider for an integer value.
func (d *DebugUI) createIntSlider(label string, value *int, min, max int, key string) *widget.Container {
	container := widget.NewContainer(
		widget.ContainerOpts.Layout(widget.NewRowLayout(
			widget.RowLayoutOpts.Direction(widget.DirectionHorizontal),
			widget.RowLayoutOpts.Spacing(10),
		)),
	)

	// Label
	labelWidget := widget.NewText(
		widget.TextOpts.Text(label, &d.fontFace, color.RGBA{200, 200, 200, 255}),
		widget.TextOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(100, 0),
		),
	)
	container.AddChild(labelWidget)

	// Value display
	valueLabel := widget.NewText(
		widget.TextOpts.Text(fmt.Sprintf("%d", *value), &d.fontFace, color.RGBA{255, 255, 255, 255}),
		widget.TextOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(40, 0),
		),
	)
	d.noiseLabels[key] = valueLabel

	// Slider
	slider := widget.NewSlider(
		widget.SliderOpts.Direction(widget.DirectionHorizontal),
		widget.SliderOpts.MinMax(min, max),
		widget.SliderOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(120, 24),
			widget.WidgetOpts.LayoutData(widget.RowLayoutData{
				Position: widget.RowLayoutPositionCenter,
			}),
		),
		widget.SliderOpts.Images(d.createSliderImages(), d.createSliderHandleImages()),
		widget.SliderOpts.PageSizeFunc(func() int {
			return 1
		}),
		widget.SliderOpts.ChangedHandler(func(args *widget.SliderChangedEventArgs) {
			*value = args.Current
			valueLabel.Label = fmt.Sprintf("%d", *value)
			d.markDirty()
		}),
	)
	slider.Current = *value

	container.AddChild(slider)
	container.AddChild(valueLabel)

	return container
}

// createFloatSlider creates a labeled slider for a float64 value.
func (d *DebugUI) createFloatSlider(label string, value *float64, min, max float64, key string) *widget.Container {
	container := widget.NewContainer(
		widget.ContainerOpts.Layout(widget.NewRowLayout(
			widget.RowLayoutOpts.Direction(widget.DirectionHorizontal),
			widget.RowLayoutOpts.Spacing(10),
		)),
	)

	// Label
	labelWidget := widget.NewText(
		widget.TextOpts.Text(label, &d.fontFace, color.RGBA{200, 200, 200, 255}),
		widget.TextOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(100, 0),
		),
	)
	container.AddChild(labelWidget)

	// Value display
	valueLabel := widget.NewText(
		widget.TextOpts.Text(formatFloat(*value), &d.fontFace, color.RGBA{255, 255, 255, 255}),
		widget.TextOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(50, 0),
		),
	)
	if key == "seaLevel" || key == "riverSrcElev" || key == "riverProb" || key == "lakeProb" || key == "oceanProb" {
		d.hydroLabels[key] = valueLabel
	} else {
		d.noiseLabels[key] = valueLabel
	}

	// Convert to int range (0-100) for slider granularity
	sliderMin := 0
	sliderMax := 100

	slider := widget.NewSlider(
		widget.SliderOpts.Direction(widget.DirectionHorizontal),
		widget.SliderOpts.MinMax(sliderMin, sliderMax),
		widget.SliderOpts.WidgetOpts(
			widget.WidgetOpts.MinSize(120, 24),
			widget.WidgetOpts.LayoutData(widget.RowLayoutData{
				Position: widget.RowLayoutPositionCenter,
			}),
		),
		widget.SliderOpts.Images(d.createSliderImages(), d.createSliderHandleImages()),
		widget.SliderOpts.PageSizeFunc(func() int {
			return 1
		}),
		widget.SliderOpts.ChangedHandler(func(args *widget.SliderChangedEventArgs) {
			// Map slider position (0-100) to actual value range
			t := float64(args.Current-sliderMin) / float64(sliderMax-sliderMin)
			*value = min + t*(max-min)
			valueLabel.Label = formatFloat(*value)
			d.markDirty()
		}),
	)

	// Set initial slider position
	t := (*value - min) / (max - min)
	slider.Current = sliderMin + int(t*float64(sliderMax-sliderMin))

	container.AddChild(slider)
	container.AddChild(valueLabel)

	return container
}

// createSliderImages creates the track images for a slider.
func (d *DebugUI) createSliderImages() *widget.SliderTrackImage {
	idle := ebiten.NewImage(32, 8)
	idle.Fill(color.RGBA{80, 80, 100, 255})

	hover := ebiten.NewImage(32, 8)
	hover.Fill(color.RGBA{100, 100, 120, 255})

	return &widget.SliderTrackImage{
		Idle:  image.NewNineSliceSimple(idle, 4, 4),
		Hover: image.NewNineSliceSimple(hover, 4, 4),
	}
}

// createSliderHandleImages creates the handle images for a slider.
func (d *DebugUI) createSliderHandleImages() *widget.ButtonImage {
	idle := ebiten.NewImage(20, 20)
	idle.Fill(color.RGBA{150, 150, 180, 255})

	hover := ebiten.NewImage(20, 20)
	hover.Fill(color.RGBA{180, 180, 220, 255})

	pressed := ebiten.NewImage(20, 20)
	pressed.Fill(color.RGBA{200, 200, 255, 255})

	return &widget.ButtonImage{
		Idle:    image.NewNineSliceSimple(idle, 4, 4),
		Hover:   image.NewNineSliceSimple(hover, 4, 4),
		Pressed: image.NewNineSliceSimple(pressed, 4, 4),
	}
}

// Toggle toggles the visibility of the debug panel.
func (d *DebugUI) Toggle() {
	d.visible = !d.visible
}

// markDirty marks the parameters as changed, triggering debounced regeneration.
func (d *DebugUI) markDirty() {
	d.dirty = true
	d.lastChangeTime = time.Now()
}

// IsVisible returns whether the debug panel is visible.
func (d *DebugUI) IsVisible() bool {
	return d.visible
}

// Update updates the UI state and handles debounced regeneration.
func (d *DebugUI) Update() {
	if d.visible {
		d.ui.Update()
	}

	// Check for debounced regeneration
	if d.dirty && time.Since(d.lastChangeTime) >= d.debounceDelay {
		d.dirty = false
		if d.OnRegenerate != nil {
			d.OnRegenerate(d.NoiseParams, d.HydroConfig)
		}
	}
}

// Draw draws the UI if visible.
func (d *DebugUI) Draw(screen *ebiten.Image) {
	if d.visible {
		d.ui.Draw(screen)
	}
}

// formatFloat formats a float for display.
func formatFloat(v float64) string {
	s := strconv.FormatFloat(v, 'f', 3, 64)
	// Trim trailing zeros but keep at least one decimal
	for len(s) > 1 && s[len(s)-1] == '0' && s[len(s)-2] != '.' {
		s = s[:len(s)-1]
	}
	return s
}
