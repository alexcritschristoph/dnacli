package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"unicode"
)

// The exported picture is the frame already on screen: the same grid of cells,
// drawn again as text and shapes so that it stays sharp at any size. Cell sizes
// follow the font, so every row and column lands where the terminal had it.
const (
	svgFont     = 13.0            // the size letters are drawn at, in pixels
	svgCellW    = 7.8             // pixels across one character cell
	svgCellH    = 16.0            // pixels down one row
	svgPad      = 12.0            // margin around the picture
	svgBaseline = svgCellH * 0.78 // where the bottom of a letter sits in its cell
)

// svgTheme is the paper the picture is drawn on. The light theme keeps every
// hue but darkens the letters and tints the panels, so that it prints.
type svgTheme struct {
	light  bool
	canvas string
}

var (
	darkSVG  = svgTheme{canvas: colorBackdrop}
	lightSVG = svgTheme{light: true, canvas: "#ffffff"}
)

// ink turns a color of the screen into one that reads on white: the same hue,
// darkened until it contrasts with the page. The palettes keep their character
// rather than fading to greys.
func (t svgTheme) ink(color string) string {
	if !t.light || color == "" {
		return color
	}
	if shade, ok := lightInks[color]; ok {
		return shade
	}
	h, s, l := hsl(color)
	for l > 0.16 && contrastWhite(hslColor(h, s, l)) < 4.5 {
		l -= 0.03
	}
	lightInks[color] = hslColor(h, s, l)
	return lightInks[color]
}

// paper tints a background just enough to keep the panels and strips apart
// from the page, keeping the hue of the color it came from.
func (t svgTheme) paper(color string) string {
	if !t.light || color == "" || color == colorBackdrop {
		return t.canvas
	}
	if shade, ok := lightPapers[color]; ok {
		return shade
	}
	h, s, l := hsl(color)
	shade := hslColor(h, s*0.5, 0.84+0.6*l)
	lightPapers[color] = shade
	return shade
}

var (
	lightInks   = map[string]string{}
	lightPapers = map[string]string{}
)

// contrastWhite is how strongly a color stands out on a white page, as the
// ratio the Web Content Accessibility Guidelines use.
func contrastWhite(color string) float64 {
	r, g, b := hex(color)
	lum := 0.2126*linear(float64(r)/255) + 0.7152*linear(float64(g)/255) + 0.0722*linear(float64(b)/255)
	return 1.05 / (lum + 0.05)
}

func linear(v float64) float64 {
	if v <= 0.03928 {
		return v / 12.92
	}
	return math.Pow((v+0.055)/1.055, 2.4)
}

// hsl splits a color into hue, saturation and lightness, each from 0 to 1.
func hsl(color string) (h, s, l float64) {
	r, g, b := hex(color)
	hi, lo := float64(max(r, max(g, b)))/255, float64(min(r, min(g, b)))/255
	l = (hi + lo) / 2
	if hi == lo {
		return 0, 0, l
	}
	d := hi - lo
	if l > 0.5 {
		s = d / (2 - hi - lo)
	} else {
		s = d / (hi + lo)
	}
	switch hi {
	case float64(r) / 255:
		h = (float64(g) - float64(b)) / 255 / d
		if g < b {
			h += 6
		}
	case float64(g) / 255:
		h = (float64(b)-float64(r))/255/d + 2
	default:
		h = (float64(r)-float64(g))/255/d + 4
	}
	return h / 6, s, l
}

// hslColor turns a hue, saturation and lightness back into a hex color.
func hslColor(h, s, l float64) string {
	if s == 0 {
		v := int(l*255 + 0.5)
		return fmt.Sprintf("#%02x%02x%02x", v, v, v)
	}
	q := l * (1 + s)
	if l >= 0.5 {
		q = l + s - l*s
	}
	p := 2*l - q
	to := func(t float64) int {
		if t < 0 {
			t++
		}
		if t > 1 {
			t--
		}
		switch {
		case t < 1.0/6:
			return int((p+(q-p)*6*t)*255 + 0.5)
		case t < 0.5:
			return int(q*255 + 0.5)
		case t < 2.0/3:
			return int((p+(q-p)*(2.0/3-t)*6)*255 + 0.5)
		}
		return int(p*255 + 0.5)
	}
	return fmt.Sprintf("#%02x%02x%02x", to(h+1.0/3), to(h), to(h-1.0/3))
}

// runExport saves the view as an SVG. The word light or dark anywhere in the
// argument picks the paper; whatever is left is the file to write.
func (m *model) runExport(arg string) {
	theme := lightSVG
	var words []string
	for _, word := range strings.Fields(arg) {
		switch strings.ToLower(word) {
		case "light":
			theme = lightSVG
		case "dark":
			theme = darkSVG
		default:
			words = append(words, word)
		}
	}
	path := expandPath(strings.Join(words, " "))
	if path == "" {
		path = m.exportName()
	}
	if m.width < 40 || m.height < 18 {
		m.notify("The window is too small to export.")
		return
	}
	if err := m.writeSVG(path, theme); err != nil {
		m.notify("Cannot write %s: %v", filepath.Base(path), err)
		return
	}
	paper := "dark"
	if theme.light {
		paper = "light"
	}
	m.notify("Saved %s, %s background.", filepath.Base(path), paper)
}

// exportName is what an export is called when no file is given: the contig and
// the visible range, with anything a file name cannot hold replaced.
func (m *model) exportName() string {
	c := m.contigs()[m.contig]
	name := strings.Map(func(r rune) rune {
		if unicode.IsLetter(r) || unicode.IsDigit(r) || strings.ContainsRune("-_.", r) {
			return r
		}
		return '-'
	}, c.Name)
	return fmt.Sprintf("dnacli-%s-%d-%d.svg", name, m.view.start+1, min(c.Length, m.view.start+m.span()))
}

// contentRect is the part of the frame a picture keeps: the ruler, genes,
// variants, reference, protein, depth and reads, without the bars, tabs and
// panels around them, and without the empty rows under the last read.
func (m *model) contentRect(l layout) (top, bottom, right int) {
	top, bottom, right = trackTop, l.pileupTop+l.pileupRows, gutterWidth+l.columns
	for bottom > top && m.blankRow(bottom-1, right) {
		bottom--
	}
	return top, bottom, right
}

func (m *model) blankRow(row, right int) bool {
	for col := 0; col < right; col++ {
		if m.frame[row][col].r != ' ' {
			return false
		}
	}
	return true
}

func (m *model) writeSVG(path string, theme svgTheme) error {
	m.reload() // make sure the reads and columns match the view being exported
	m.exporting = true
	m.View() // the frame on screen is the frame that gets exported
	m.exporting = false
	l := m.layout()
	top, bottom, right := m.contentRect(l)
	rows, cols := bottom-top, right
	contig := m.contigs()[m.contig]
	width := 2*svgPad + float64(cols)*svgCellW
	height := 2*svgPad + float64(rows)*svgCellH

	var b strings.Builder
	b.WriteString(`<?xml version="1.0" encoding="UTF-8"?>` + "\n")
	fmt.Fprintf(&b, `<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f"`+
		` font-family="'DejaVu Sans Mono', 'Liberation Mono', Menlo, Consolas, monospace" font-size="%.0f">`+"\n",
		width, height, width, height, svgFont)
	fmt.Fprintf(&b, "<title>dnacli %s:%d-%d</title>\n",
		xmlText(contig.Name), m.view.start+1, min(contig.Length, m.view.start+m.span()))
	fmt.Fprintf(&b, `<rect width="100%%" height="100%%" fill="%s"/>`+"\n", theme.canvas)

	b.WriteString(`<g id="background" shape-rendering="crispEdges">` + "\n")
	for i := 0; i < rows; i++ {
		row := top + i
		for col := 0; col < cols; {
			bg := theme.paper(m.frame[row][col].s.bg)
			start := col
			for col < cols && theme.paper(m.frame[row][col].s.bg) == bg {
				col++
			}
			if bg != theme.canvas {
				b.WriteString(rect(svgPad+float64(start)*svgCellW, svgPad+float64(i)*svgCellH,
					float64(col-start)*svgCellW, svgCellH, bg) + "\n")
			}
		}
	}
	b.WriteString("</g>\n")

	b.WriteString(`<g id="symbols" shape-rendering="crispEdges">` + "\n")
	for i := 0; i < rows; i++ {
		row := top + i
		for col := 0; col < cols; {
			cell := m.frame[row][col]
			if full, isBlock := blockIndex(cell.r); isBlock {
				run := 1
				for col+run < cols && m.frame[row][col+run].r == cell.r && m.frame[row][col+run].s.fg == cell.s.fg {
					run++
				}
				h := svgCellH * float64(full) / 8
				b.WriteString(rect(svgPad+float64(col)*svgCellW, svgPad+float64(i)*svgCellH+svgCellH-h,
					float64(run)*svgCellW, h, theme.ink(cell.s.fg)) + "\n")
				col += run
				continue
			}
			if shape := glyphShape(cell.r, col, i, theme.ink(cell.s.fg)); shape != "" {
				b.WriteString(shape + "\n")
			}
			col++
		}
	}
	b.WriteString("</g>\n")

	b.WriteString(`<g id="text" xml:space="preserve">` + "\n")
	for i := 0; i < rows; i++ {
		row := top + i
		for col := 0; col < cols; {
			first := m.frame[row][col]
			if !drawsText(first.r) {
				col++
				continue
			}
			start, run := col, []rune{}
			for col < cols {
				next := m.frame[row][col]
				if next.s != first.s || !drawsText(next.r) {
					break
				}
				run = append(run, next.r)
				col++
			}
			fmt.Fprintf(&b, `<text x="%.1f" y="%.1f" fill="%s"`, svgPad+float64(start)*svgCellW,
				svgPad+float64(i)*svgCellH+svgBaseline, theme.ink(first.s.fg))
			if first.s.bold {
				b.WriteString(` font-weight="bold"`)
			}
			if first.s.faint {
				b.WriteString(` fill-opacity="0.65"`)
			}
			if len(run) > 1 {
				fmt.Fprintf(&b, ` textLength="%.1f" lengthAdjust="spacingAndGlyphs"`, float64(len(run))*svgCellW)
			}
			fmt.Fprintf(&b, ">%s</text>\n", xmlText(string(run)))
		}
	}
	b.WriteString("</g>\n</svg>\n")

	return os.WriteFile(path, []byte(b.String()), 0o644)
}

// drawsText is true for a cell that is drawn as a letter rather than a shape or
// a patch of background.
func drawsText(r rune) bool {
	return r != 0 && r != ' ' && unicode.IsPrint(r) && glyphShape(r, 0, 0, "#000000") == ""
}

// blockIndex says how full a block character is: 1 for the lowest eighth, 8 for
// a full cell.
func blockIndex(r rune) (int, bool) {
	for i, block := range []rune("▁▂▃▄▅▆▇█") {
		if block == r {
			return i + 1, true
		}
	}
	return 0, false
}

// glyphShape draws the block, line, arrow and marker characters as shapes, so
// that the picture never depends on the fonts of whoever opens it.
func glyphShape(r rune, col, row int, fill string) string {
	x, y := svgPad+float64(col)*svgCellW, svgPad+float64(row)*svgCellH
	if i, isBlock := blockIndex(r); isBlock {
		h := svgCellH * float64(i) / 8
		return rect(x, y+svgCellH-h, svgCellW, h, fill)
	}
	switch r {
	case '─':
		return rect(x, y+svgCellH/2-0.6, svgCellW, 1.2, fill)
	case '━':
		return rect(x, y+svgCellH/2-1.2, svgCellW, 2.4, fill)
	case '│':
		return rect(x+svgCellW/2-0.6, y, 1.2, svgCellH, fill)
	case '·':
		return fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="1.4" fill="%s"/>`,
			x+svgCellW/2, y+svgCellH/2, fill)
	case '▪':
		return rect(x+svgCellW/2-1.8, y+svgCellH/2-1.8, 3.6, 3.6, fill)
	case '▼':
		return polygon(fill, x+svgCellW*0.15, y+svgCellH*0.3, x+svgCellW*0.85, y+svgCellH*0.3,
			x+svgCellW*0.5, y+svgCellH*0.78)
	case '▶':
		return polygon(fill, x+svgCellW*0.2, y+svgCellH*0.22, x+svgCellW*0.2, y+svgCellH*0.78,
			x+svgCellW*0.85, y+svgCellH*0.5)
	case '◀':
		return polygon(fill, x+svgCellW*0.8, y+svgCellH*0.22, x+svgCellW*0.8, y+svgCellH*0.78,
			x+svgCellW*0.15, y+svgCellH*0.5)
	case '▸':
		return polygon(fill, x+svgCellW*0.3, y+svgCellH*0.32, x+svgCellW*0.3, y+svgCellH*0.68,
			x+svgCellW*0.7, y+svgCellH*0.5)
	}
	return ""
}

func rect(x, y, w, h float64, fill string) string {
	return fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>`, x, y, w, h, fill)
}

func polygon(fill string, points ...float64) string {
	var b strings.Builder
	fmt.Fprintf(&b, `<polygon fill="%s" points="`, fill)
	for i := 0; i+1 < len(points); i += 2 {
		if i > 0 {
			b.WriteString(" ")
		}
		fmt.Fprintf(&b, "%.1f,%.1f", points[i], points[i+1])
	}
	b.WriteString(`"/>`)
	return b.String()
}

var xmlSafe = strings.NewReplacer("&", "&amp;", "<", "&lt;", ">", "&gt;")

func xmlText(s string) string { return xmlSafe.Replace(s) }
