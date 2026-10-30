package main

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/muesli/termenv"
)

const (
	colorBackdrop = "#0f1117"
	colorPanel    = "#161a23"
	colorStripe   = "#1b2130" // the shade of every other strip along the bottom
	colorLine     = "#2a3040"
	colorDim      = "#5b6472"
	colorMuted    = "#8b93a1"
	colorText     = "#d7dce5"
	colorBright   = "#f2f5fa"
	colorAccent   = "#7aa2f7"
	colorWarn     = "#e0af68"
	colorAlert    = "#f7768e"
	colorGood     = "#9ece6a"
)

// palette is the set of colors every coloring scheme draws from. A palette
// changes the shades, never their meaning: a mismatch still takes its base's
// color, and a stop codon still takes the alarm color.
type palette struct {
	name                     string
	bases                    map[byte]string
	forward, reverse         string
	dim, cool, loud          string
	good, warn, alert, other string
}

var palettes = []palette{
	{
		name:    "Default",
		bases:   map[byte]string{'A': "#7ee787", 'C': "#79c0ff", 'G': "#ffb454", 'T': "#ff7b72", 'N': "#6e7681", '-': "#454c5c"},
		forward: "#5b8fd6", reverse: "#c98fb0",
		dim: "#41475a", cool: "#9dc9f0", loud: "#f7768e",
		good: "#9ece6a", warn: "#e0af68", alert: "#f7768e", other: "#bb9af7",
	},
	{
		name:    "Color-blind safe",
		bases:   map[byte]string{'A': "#009e73", 'C': "#56b4e9", 'G': "#e69f00", 'T': "#d55e00", 'N': "#9a9a9a", '-': "#4b4b4b"},
		forward: "#0072b2", reverse: "#cc79a7",
		dim: "#3f4653", cool: "#56b4e9", loud: "#d55e00",
		good: "#009e73", warn: "#e69f00", alert: "#d55e00", other: "#cc79a7",
	},
	{
		name:    "High contrast",
		bases:   map[byte]string{'A': "#00ff87", 'C': "#00d7ff", 'G': "#ffd700", 'T': "#ff4d4d", 'N': "#c0c0c0", '-': "#6a6a6a"},
		forward: "#4d8dff", reverse: "#ff6ec7",
		dim: "#2b3140", cool: "#00e5ff", loud: "#ff2d55",
		good: "#00ff87", warn: "#ffb300", alert: "#ff2d55", other: "#c77dff",
	},
	{
		name:    "Sunset",
		bases:   map[byte]string{'A': "#f6c177", 'C': "#9ccfd8", 'G': "#c4a7e7", 'T': "#eb6f92", 'N': "#908caa", '-': "#524f67"},
		forward: "#3e8fb0", reverse: "#eb6f92",
		dim: "#2a273f", cool: "#9ccfd8", loud: "#eb6f92",
		good: "#a3be8c", warn: "#f6c177", alert: "#eb6f92", other: "#c4a7e7",
	},
	{
		name:    "Meadow",
		bases:   map[byte]string{'A': "#a3be8c", 'C': "#88c0d0", 'G': "#ebcb8b", 'T': "#bf616a", 'N': "#8f9aa8", '-': "#4c566a"},
		forward: "#81a1c1", reverse: "#b48ead",
		dim: "#3b4252", cool: "#88c0d0", loud: "#bf616a",
		good: "#a3be8c", warn: "#ebcb8b", alert: "#bf616a", other: "#b48ead",
	},
}

// chosenPalette is the palette in force. The view is drawn once, so there is
// one palette: choosing another repaints every read on the next frame.
var chosenPalette = palettes[0]

// cellStyle is a comparable description of one screen cell, so that styles can
// be cached and equal neighbours merged into a single escape sequence.
type cellStyle struct {
	fg, bg string
	bold   bool
	faint  bool
}

type cell struct {
	r rune
	s cellStyle
}

var profile = termenv.ColorProfile()

var escapes = map[cellStyle]string{}

// escape is the sequence that switches the terminal to a style. Terminals that
// cannot show 24-bit color get the nearest color they can.
func (c cellStyle) escape() string {
	if e, ok := escapes[c]; ok {
		return e
	}
	if profile == termenv.Ascii {
		escapes[c] = "" // a terminal that cannot color gets plain text
		return ""
	}
	parts := []string{"0"}
	if c.bold {
		parts = append(parts, "1")
	}
	if c.faint {
		parts = append(parts, "2")
	}
	if c.fg != "" {
		parts = append(parts, profile.Color(c.fg).Sequence(false))
	}
	if c.bg != "" {
		parts = append(parts, profile.Color(c.bg).Sequence(true))
	}
	escapes[c] = "\x1b[" + strings.Join(parts, ";") + "m"
	return escapes[c]
}

// renderRow turns cells into a line of text, changing style only where the
// style changes.
func renderRow(row []cell) string {
	var out strings.Builder
	var current cellStyle
	for i, c := range row {
		if i == 0 || c.s != current {
			out.WriteString(c.s.escape())
			current = c.s
		}
		out.WriteRune(c.r)
	}
	out.WriteString(cellStyle{}.escape()) // back to plain text at the end of the line
	return out.String()
}

// highlights pick out the bases that differ from something. They form one
// section of the coloring panel, and can be left off.
var highlights = []string{
	"Mismatch (from consensus)",
	"Mismatch (from reference)",
}

// colors give every base the scheme's own color. They form the other section
// of the panel, and can be left off too.
var colors = []string{
	"Nucleotides",
	"Read strand",
	"Mapping quality",
	"Base quality",
	"Mismatch count",
	"Distance between mates",
	"Reads with indels",
	"Read length",
	"Protein change",
}

// namedIndex finds a choice by the name the panel and --color use, or -1.
func namedIndex(list []string, name string) int {
	for i, choice := range list {
		if strings.EqualFold(choice, name) {
			return i
		}
	}
	return -1
}

// quietGrey is the neutral shade of a base that has no color scheme of its own.
const quietGrey = "#575e6d"

// readColor is the color of a base that is not a difference: the chosen
// scheme's own color, or quiet grey when no scheme is chosen. typical is the
// usual read length here, which the length scheme shades against.
func readColor(color int, a *Alignment, pos, typical int) string {
	if color < 0 {
		return quietGrey
	}
	switch colors[color] {
	case "Nucleotides":
		if a.Bases == nil {
			return quietGrey // zoomed out, where one cell covers many bases
		}
		return chosenPalette.bases[a.BaseAt(pos)]
	case "Read strand":
		if a.Reverse() {
			return chosenPalette.reverse
		}
		return chosenPalette.forward
	case "Mapping quality":
		return confidence(float64(min(int(a.MapQ), 60)) / 60)
	case "Base quality":
		if a.Quals == nil {
			return quietGrey // zoomed out, where one cell covers many bases
		}
		q := a.QualAt(pos)
		if q == 0xff {
			return quietGrey // the file gives this base no score
		}
		return confidence(float64(min(int(q), 40)) / 40)
	case "Mismatch count":
		return alarm(float64(min(a.Mismatches, 10)) / 10)
	case "Distance between mates":
		if !a.Paired() || a.InsertSize == 0 {
			return quietGrey
		}
		return alarm(float64(min(abs(a.InsertSize), 1000)) / 1000)
	case "Reads with indels":
		if len(a.Insertions) > 0 || a.HasIndel {
			return colorWarn
		}
		return quietGrey
	case "Read length":
		return confidence(float64(min(len(a.Sequence), 2*typical)) / float64(max(2*typical, 1)))
	}
	return quietGrey
}

// confidence shades from the palette's dim shade at 0 to its cool one at 1.
func confidence(t float64) string { return mix(chosenPalette.dim, chosenPalette.cool, t) }

// alarm shades from quiet grey at 0 to the palette's loud shade at 1.
func alarm(t float64) string { return mix(quietGrey, chosenPalette.loud, t) }

// mix blends two colors, remembering shades it has already worked out because
// a screen full of reads asks for the same ones over and over.
var blends = map[blend]string{}

type blend struct {
	from, to string
	step     int // t, in sixty-fourths
}

func mix(from, to string, t float64) string {
	key := blend{from, to, int(t * 63)}
	if color, ok := blends[key]; ok {
		return color
	}
	fr, fg, fb := hex(from)
	tr, tg, tb := hex(to)
	at := func(a, b int) int { return a + (b-a)*key.step/63 }
	blends[key] = fmt.Sprintf("#%02x%02x%02x", at(fr, tr), at(fg, tg), at(fb, tb))
	return blends[key]
}

func hex(s string) (int, int, int) {
	n, _ := strconv.ParseUint(strings.TrimPrefix(s, "#"), 16, 32)
	return int(n >> 16), int(n>>8) & 0xff, int(n) & 0xff
}

func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}
