package main

import (
	"fmt"
	"strings"
)

var sparkChars = []rune("▁▂▃▄▅▆▇█")

// coloringLine is one line of the coloring panel: a section title, an entry's
// name, or a blank line between the sections. entry is the choice's place in
// both lists together, and -1 on anything but a name.
type coloringLine struct {
	title string
	name  string
	entry int
}

// coloringPanel lays the coloring panel out once, so that drawing a frame and
// clicking a row always agree on which line is which.
var coloringPanel = coloringLayout()

func coloringLayout() []coloringLine {
	var lines []coloringLine
	entry := 0
	paletteNames := make([]string, len(palettes))
	for i, p := range palettes {
		paletteNames[i] = p.name
	}
	sections := []struct {
		title string
		list  []string
	}{{"SNP highlighting", highlights}, {"Color by", colors}, {"Color palette", paletteNames}}
	for i, section := range sections {
		if i > 0 {
			lines = append(lines, coloringLine{entry: -1}) // a blank line between the sections
		}
		lines = append(lines, coloringLine{title: section.title, entry: -1})
		for _, name := range section.list {
			lines = append(lines, coloringLine{name: name, entry: entry})
			entry++
		}
	}
	return lines
}

func sparkline(profile []float64, width int) string {
	if len(profile) == 0 {
		return strings.Repeat("·", width)
	}
	peak := 0.0
	for _, v := range profile {
		peak = max(peak, v)
	}
	var b strings.Builder
	for i := 0; i < width; i++ {
		v := profile[i*len(profile)/width]
		switch {
		case peak == 0:
			b.WriteRune('·')
		case v == 0:
			b.WriteRune('·')
		default:
			b.WriteRune(sparkChars[min(int(v/peak*float64(len(sparkChars))), len(sparkChars)-1)])
		}
	}
	return b.String()
}

func (m *model) drawPanel(g grid, l layout) {
	left := m.width - l.panel
	back := cellStyle{fg: colorText, bg: colorPanel}
	title := cellStyle{fg: colorAccent, bg: colorPanel, bold: true}
	for row := trackTop; row < m.height-l.barRows; row++ { // between the command bar and the strips
		g.fill(row, left, m.width, ' ', back)
		g.draw(row, left, "│", cellStyle{fg: colorLine, bg: colorPanel})
	}
	inner := left + 2
	width := l.panel - 3
	top, bottom := l.listTop, l.listTop+l.listRows

	switch m.panel {
	case "variants":
		g.draw(trackTop, inner, fmt.Sprintf("Variants  %d", len(m.variants)), title)
		m.drawVariantList(g, inner, width, top, bottom)
	case "contigs":
		g.draw(trackTop, inner, fmt.Sprintf("Contigs  %d", len(m.contigs())), title)
		m.drawContigList(g, inner, width, top, bottom)
	case "colors":
		g.draw(trackTop, inner, "Read coloring", title)
		offset := m.coloringOffset(l)
		for i, line := range coloringPanel[offset:] {
			row := top + i
			if row >= bottom { // the rest do not fit in the window
				break
			}
			switch {
			case line.title != "":
				g.draw(row, inner, line.title, cellStyle{fg: colorMuted, bg: colorPanel, bold: true})
			case line.entry < 0: // a blank line between the sections
			default:
				cursor := line.entry == m.schemeRow
				style := cellStyle{fg: colorText, bg: colorPanel}
				switch {
				case m.chosen(line.entry):
					style = cellStyle{fg: colorBright, bg: colorLine, bold: true}
				case cursor: // the cursor is here, but this entry is not the chosen one
					style = cellStyle{fg: colorText, bg: colorStripe}
				}
				if style.bg != colorPanel {
					g.fill(row, inner-1, m.width-1, ' ', style)
				}
				mark := "  "
				if cursor {
					mark = "▸ "
				}
				g.draw(row, inner, mark+line.name, style)
			}
		}
		hint := "↑ ↓ move   enter choose or clear"
		if !m.coloringLocked {
			hint = "↑ ↓ preview   enter locks"
		}
		g.draw(bottom, inner, hint, cellStyle{fg: colorDim, bg: colorPanel})
	case "filters":
		g.draw(trackTop, inner, "Filters", title)
		for i, r := range filterRows {
			style := back
			if m.panelFocus && i == m.filterIndex {
				style = cellStyle{fg: colorBright, bg: colorLine}
				g.fill(top+i*2, inner-1, m.width-1, ' ', style)
			}
			value := r.value(m.filter)
			g.draw(top+i*2, inner, r.name, style)
			g.draw(top+i*2, inner+width-len(value), value, cellStyle{fg: colorBright, bg: style.bg, bold: true})
		}
		g.draw(top+len(filterRows)*2+1, inner, "↑ ↓ choose   ← → change", cellStyle{fg: colorDim, bg: colorPanel})
	}
}

// filterRows are the read filters as the panel lists them, and how the left and
// right keys change each one.
var filterRows = []struct {
	name  string
	value func(Filter) string
	move  func(*Filter, int)
}{
	{"Mapping quality at least", func(f Filter) string { return fmt.Sprint(f.MinMapQ) },
		func(f *Filter, d int) { f.MinMapQ = clamp(f.MinMapQ+5*d, 0, 60) }},
	{"Average base quality at least", func(f Filter) string { return fmt.Sprint(f.MinMeanQual) },
		func(f *Filter, d int) { f.MinMeanQual = clamp(f.MinMeanQual+5*d, 0, 40) }},
	{"Mismatches at most", mismatchLimit,
		func(f *Filter, d int) { f.MaxMismatches = clamp(f.MaxMismatches+d, -1, 50) }},
	{"Secondary alignments", func(f Filter) string { return shown(!f.HideSecondary) },
		func(f *Filter, _ int) { f.HideSecondary = !f.HideSecondary }},
	{"Duplicates", func(f Filter) string { return shown(!f.HideDuplicates) },
		func(f *Filter, _ int) { f.HideDuplicates = !f.HideDuplicates }},
	{"Reads not in a proper pair", func(f Filter) string { return shown(!f.HideUnpaired) },
		func(f *Filter, _ int) { f.HideUnpaired = !f.HideUnpaired }},
	{"Read length at least", func(f Filter) string { return comma(f.MinLength) },
		func(f *Filter, d int) { f.MinLength = along(lengthLadder, f.MinLength, d) }},
	{"Soft-clipped bases", func(f Filter) string { return shown(f.SoftClips) },
		func(f *Filter, _ int) { f.SoftClips = !f.SoftClips }},
}

// lengthLadder is what the arrow keys walk through for the read length filter,
// which has to cover both short reads and long ones.
var lengthLadder = []int{0, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000}

// along moves one place through a ladder of values.
func along(ladder []int, value, delta int) int {
	at := 0
	for i, v := range ladder {
		if v <= value {
			at = i
		}
	}
	return ladder[clamp(at+delta, 0, len(ladder)-1)]
}

func shown(on bool) string {
	if on {
		return "shown"
	}
	return "hidden"
}

func (m *model) drawVariantList(g grid, left, width, top, bottom int) {
	if len(m.variants) == 0 {
		g.draw(top, left, "Load a VCF file to list variants.", cellStyle{fg: colorDim, bg: colorPanel})
		return
	}
	visible := bottom - top
	offset := max(0, min(m.variantIndex-visible/2, len(m.variants)-visible))
	for i := 0; i < visible && offset+i < len(m.variants); i++ {
		v := m.variants[offset+i]
		row := top + i
		style := cellStyle{fg: colorText, bg: colorPanel}
		if offset+i == m.variantIndex {
			style = cellStyle{fg: colorBright, bg: colorLine, bold: true}
			g.fill(row, left-1, m.width-1, ' ', style)
		}
		mark := cellStyle{fg: impactColor(v.Impact), bg: style.bg, bold: true}
		g.draw(row, left, v.Impact.Short(), mark)
		g.draw(row, left+2, comma(v.Start+1), style)
		g.draw(row, left+12, shorten(v.Label(), 9), style)
		kind := v.Impact.Kind
		if v.Filter != "" && v.Filter != "PASS" && v.Filter != "." {
			kind = v.Filter
		}
		g.draw(row, left+22, shorten(kind, width-22), cellStyle{fg: colorMuted, bg: style.bg})
	}
	g.draw(bottom, left, "↑ ↓ to jump between them", cellStyle{fg: colorDim, bg: colorPanel})
}

func (m *model) drawContigList(g grid, left, width, top, bottom int) {
	contigs := m.contigs()
	visible := (bottom - top) / 3
	offset := max(0, min(m.contigIndex-visible/2, len(contigs)-visible))
	for i := 0; i < visible && offset+i < len(contigs); i++ {
		c := contigs[offset+i]
		row := top + i*3
		style := cellStyle{fg: colorText, bg: colorPanel}
		if offset+i == m.contigIndex {
			style = cellStyle{fg: colorBright, bg: colorLine, bold: true}
			g.fill(row, left-1, m.width-1, ' ', style)
			g.fill(row+1, left-1, m.width-1, ' ', style)
		}
		g.draw(row, left, shorten(c.Name, width-14), style)
		count := fmt.Sprintf("%s reads", comma(int(c.Reads)))
		g.draw(row, left+width-len(count), count, cellStyle{fg: colorMuted, bg: style.bg})
		g.draw(row+1, left, sparkline(c.Profile, width-10), cellStyle{fg: "#5b8fd6", bg: style.bg})
		size := comma(c.Length) + " bp"
		g.draw(row+1, left+width-len(size), size, cellStyle{fg: colorDim, bg: style.bg})
	}
	g.draw(bottom, left, "↑ ↓ to open one", cellStyle{fg: colorDim, bg: colorPanel})
}

var helpLines = [][2]string{
	{"Command bar keywords (just type)", ""},
	{"go 1000", "jump to position on contig"},
	{"gene [name]", "jump to a gene by name, if a Genbank is provided"},
	{"open [file]", "open another BAM, GenBank or VCF file"},
	{"export [file]", "save this view as an SVG"},
	{"contigs", "view and switch between contigs"},
	{"colors", "change read coloring and highlighting"},
	{"filters", "filter which reads are shown"},
	{"variants", "scan by variant, if a VCF is provided"},
	{"", ""},
	{"Navigation", ""},
	{"← →", "pan left and right"},
	{"t", "switch to the next open BAM"},
	{"↑ ↓  page up/down", "scroll reads vertically"},
	{"home  end", "jump to the ends of the contig"},
	{"+ -  or the wheel", "zoom in and out"},
	{"left click", "left clicking on a read copies it"},
	{"0", "fit the whole contig on screen"},
	{"esc", "dismiss an open menu or panel"},
	{"m", "turn mouse on/off, makes text selectable"},
	{"?", "this help"},
	{"q", "quit"},
}

func (m *model) drawHelp(g grid) {
	width, height := 0, len(helpLines)+4
	for _, line := range helpLines {
		width = max(width, 28+len([]rune(line[1])))
	}
	width = min(width, m.width-4)
	left, top := max(0, (m.width-width)/2), max(2, (m.height-height)/2)
	back := cellStyle{fg: colorText, bg: colorPanel}
	for row := top; row < min(top+height, m.height-2); row++ {
		g.fill(row, left, left+width, ' ', back)
	}
	g.draw(top+1, left+3, "dnacli - commands", cellStyle{fg: colorAccent, bg: colorPanel, bold: true})
	for i, line := range helpLines {
		row := top + 3 + i
		if row >= m.height-2 {
			break
		}
		if line[1] == "" {
			g.draw(row, left+3, line[0], cellStyle{fg: colorBright, bg: colorPanel, bold: true})
			continue
		}
		g.draw(row, left+4, line[0], cellStyle{fg: colorAccent, bg: colorPanel})
		g.draw(row, left+24, line[1], back)
	}
}
