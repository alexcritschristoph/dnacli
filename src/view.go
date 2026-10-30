package main

import (
	"fmt"
	"path/filepath"
	"sort"
	"strings"
)

const gutterWidth = 12

// The command bar sits under the title bar, and the tracks and the side panel
// start under that. Padding rows either side of the row typed into are had by
// raising commandRows.
const commandTop, commandRows = 1, 1
const commandRow = commandTop + commandRows/2 // the middle row, the one typed into
const trackTop = commandTop + commandRows

// layout is where everything sits on screen for the current window size.
type layout struct {
	panel      int // width of the side panel, 0 when closed
	columns    int // columns of reference shown
	geneRows   int
	aminoRow   int
	refRow     int
	barRows    int // everything below the pileup: the strips and the keys
	listTop    int // first row of a panel's list
	listRows   int
	pileupTop  int
	pileupRows int
}

func (m *model) layout() layout {
	l := layout{}
	if m.panel != "" && m.panel != "help" && m.width > 90 {
		l.panel = 38
	}
	l.columns = max(8, m.width-gutterWidth-l.panel-1)
	l.geneRows = len(m.geneTrack(l.columns))
	if m.view.bpp == 1 && len(m.codingTrack(l.columns)) > 0 {
		l.aminoRow = 1
	}
	l.refRow = trackTop + 1 + l.geneRows + 1    // title bar, command bar, ruler, genes, variants
	l.pileupTop = l.refRow + 1 + l.aminoRow + 1 // reference, aminos, coverage
	l.barRows = 2                               // a blank row above the strips, and the keys below them
	for _, b := range m.bars {
		l.barRows += len(b.lines)
	}
	l.pileupRows = max(0, m.height-l.pileupTop-l.barRows)
	l.listTop = trackTop + 2 // below the panel's title
	l.listRows = m.height - l.barRows - 1 - l.listTop
	return l
}

// columnSpec says what one screen column shows: a reference base, or one of the
// bases that reads have inserted after it.
type columnSpec struct {
	pos    int
	insert int // 0 for the reference base itself, otherwise the 1st, 2nd… inserted base
}

// maxInsertColumns caps how much room one insertion may take, so that a long
// one cannot push the rest of the view off the screen.
const maxInsertColumns = 10

// buildColumns decides what each column shows, giving the insertions in view
// columns of their own so that the bases either side of them still line up.
func (m *model) buildColumns(columns int) {
	m.cols, m.inserts = make([]columnSpec, 0, columns), nil
	if m.view.bpp > 1 {
		for col := 0; col < columns; col++ {
			m.cols = append(m.cols, columnSpec{m.view.start + col*m.view.bpp, 0})
		}
		return
	}
	m.inserts = map[int]int{}
	for _, a := range m.reads {
		for pos, seq := range a.Insertions {
			if pos >= m.view.start && pos < m.view.start+columns {
				m.inserts[pos] = min(max(m.inserts[pos], len(seq)), maxInsertColumns)
			}
		}
	}
	for pos := m.view.start; len(m.cols) < columns; pos++ {
		m.cols = append(m.cols, columnSpec{pos, 0})
		for i := 1; i <= m.inserts[pos] && len(m.cols) < columns; i++ {
			m.cols = append(m.cols, columnSpec{pos, i})
		}
	}
}

// column returns the screen column showing a reference position, which is
// negative for positions off the left edge.
func (m *model) column(pos int) int {
	if pos < m.view.start {
		return (pos - m.view.start - m.view.bpp + 1) / m.view.bpp
	}
	col := (pos - m.view.start) / m.view.bpp
	for at, width := range m.inserts {
		if at < pos {
			col += width
		}
	}
	return col
}

// grid is the whole screen as styled cells.
type grid [][]cell

var blank = cell{' ', cellStyle{bg: colorBackdrop}}

func newGrid(width, height int) grid {
	cells := make([]cell, width*height)
	for i := range cells {
		cells[i] = blank
	}
	g := make(grid, height)
	for i := range g {
		g[i] = cells[i*width : (i+1)*width]
	}
	return g
}

// blankFrame clears the grid we drew the last frame into and hands it back, so
// that a mouse move does not allocate a screenful of cells every time.
func (m *model) blankFrame() grid {
	if len(m.frame) != m.height || len(m.frame[0]) != m.width {
		m.frame = newGrid(m.width, m.height)
		return m.frame
	}
	for _, row := range m.frame {
		for i := range row {
			row[i] = blank
		}
	}
	return m.frame
}

func (g grid) String() string {
	rows := make([]string, len(g))
	for i, row := range g {
		rows[i] = renderRow(row)
	}
	return strings.Join(rows, "\n")
}

// set paints one cell, ignoring anything off the edge of the screen.
func (g grid) set(row, col int, r rune, s cellStyle) {
	if row >= 0 && row < len(g) && col >= 0 && col < len(g[row]) {
		g[row][col] = cell{r, s}
	}
}

func (g grid) draw(row, col int, text string, s cellStyle) {
	for _, r := range text {
		g.set(row, col, r, s)
		col++
	}
}

func (g grid) fill(row, from, to int, r rune, s cellStyle) {
	for col := from; col < to; col++ {
		g.set(row, col, r, s)
	}
}

func (m *model) View() string {
	if m.width < 40 || m.height < 18 {
		return "Error: dnacli needs a larger window."
	}
	m.refreshBars()
	l := m.layout()
	if len(m.cols) != l.columns {
		m.buildColumns(l.columns) // the window has changed size since the last load
	}
	g := m.blankFrame()
	m.drawHeader(g, l)
	m.drawTracks(g, l)
	m.drawPileup(g, l)
	m.markColumn(g, l)
	if l.panel > 0 {
		m.drawPanel(g, l)
	}
	m.drawBars(g, l) // last, so the strips span the whole width
	if m.panel == "help" && !m.exporting {
		m.drawHelp(g)
	}
	if !m.exporting {
		m.drawSuggestions(g)
	}
	return g.String()
}

func (m *model) drawHeader(g grid, l layout) {
	head := cellStyle{fg: colorText, bg: colorPanel}
	g.fill(0, 0, m.width, ' ', head)
	badge := " dnacli "
	g.draw(0, 0, badge, cellStyle{fg: colorBackdrop, bg: colorAccent, bold: true})

	c := m.contigs()[m.contig]
	span := l.columns * m.view.bpp
	right := fmt.Sprintf("%s   %s–%s of %s ",
		c.Name, comma(m.view.start+1), comma(min(m.view.start+span, c.Length)), comma(c.Length))
	rightWidth := len([]rune(right))
	g.draw(0, m.width-rightWidth, right, cellStyle{fg: colorBright, bg: colorPanel, bold: true})

	col := m.drawTabs(g, len([]rune(badge))+1, m.width-rightWidth-11)
	if m.ref != nil {
		label := "reference " + filepath.Base(m.refPath)
		if space := m.width - rightWidth - col - 3; space > 14 {
			g.draw(0, col+2, shorten(label, space), cellStyle{fg: colorDim, bg: colorPanel})
		}
	}
	m.drawCommandBar(g)
}

// drawTabs draws one tab for each alignment file, the open one standing out.
func (m *model) drawTabs(g grid, left, room int) int {
	m.tabHits = m.tabHits[:0]
	if len(m.tabs) == 0 || room < 4 {
		return left
	}
	width := 24
	if n := len(m.tabs); n > 1 {
		width = clamp((room-3*n+1)/n, 3, 24)
	}
	fit := max(1, (room+1)/(width+3))
	first := clamp(m.tab-fit+1, 0, max(0, len(m.tabs)-fit))
	col := left
	for i := first; i < len(m.tabs) && i < first+fit; i++ {
		name := shorten(filepath.Base(m.tabs[i].Path), width)
		style := cellStyle{fg: colorDim, bg: colorPanel}
		if i == m.tab {
			style = cellStyle{fg: colorBright, bg: colorLine, bold: true}
		}
		if i > first {
			col++
		}
		g.draw(0, col, " "+name+" ", style)
		m.tabHits = append(m.tabHits, tabHit{col, col + len([]rune(name)) + 1, i})
		col += len([]rune(name)) + 2
	}
	return col
}

// drawCommandBar draws the line everything is typed into, with what has been
// typed so far and where it could go.
func (m *model) drawCommandBar(g grid) {
	bar := cellStyle{fg: colorBright, bg: colorLine}
	for row := commandTop; row < commandTop+commandRows; row++ {
		g.fill(row, 0, m.width, ' ', bar)
	}
	row := commandRow
	g.draw(row, 2, "›", cellStyle{fg: colorAccent, bg: colorLine, bold: true})
	switch {
	case m.message != "":
		note := shorten(m.message, m.width-8) + " "
		g.draw(row, m.width-len([]rune(note))-2, note, cellStyle{fg: colorGood, bg: colorLine})
	case !m.mouseOn:
		g.draw(row, m.width-26, " mouse off: select text", cellStyle{fg: colorWarn, bg: colorLine})
	}
	if m.command == "" {
		if m.message == "" { // a message gets the whole bar, so it is not hidden behind the hint
			g.draw(row, 4, "type a command: contigs, colors, filters, variants, go, gene, export", cellStyle{fg: colorDim, bg: colorLine})
		}
		return
	}
	g.draw(row, 4, m.command+"▌", bar)
}

// drawSuggestions lists what the command could become, under the bar.
func (m *model) drawSuggestions(g grid) {
	list := m.suggestions()
	if m.command == "" || len(list) == 0 {
		return
	}
	width := min(64, m.width-2)
	for i, s := range list {
		style := cellStyle{fg: colorText, bg: colorPanel}
		if i == clamp(m.choice, 0, len(list)-1) {
			style = cellStyle{fg: colorBright, bg: colorLine, bold: true}
		}
		row := commandTop + commandRows + i
		g.fill(row, 0, m.width, ' ', cellStyle{bg: colorBackdrop}) // clear the strip behind the list
		g.fill(row, 0, 1+width, ' ', style)
		g.draw(row, 4, s.label, style)
		g.draw(row, 1+width-len([]rune(s.help))-2, shorten(s.help, width-len([]rune(s.label))-8), cellStyle{fg: colorDim, bg: style.bg})
	}
}

func mismatchLimit(f Filter) string {
	if f.MaxMismatches < 0 {
		return "any"
	}
	return fmt.Sprintf("≤ %d", f.MaxMismatches)
}

// gene track ---------------------------------------------------------------

type geneBar struct {
	feature Feature
	from    int
	to      int
}

// geneTrack lays the annotations in view onto as few rows as possible.
func (m *model) geneTrack(columns int) [][]geneBar {
	if m.ref == nil {
		return nil
	}
	contig := m.contigs()[m.contig].Name
	end := m.view.start + columns*m.view.bpp
	inView := m.ref.featuresIn(contig, m.view.start, end)
	coding := map[string]bool{}
	for _, f := range inView {
		if f.Kind == "CDS" {
			coding[f.Name] = true
		}
	}
	var rows [][]geneBar
	for _, f := range inView {
		if f.End() <= m.view.start || f.Start() >= end {
			continue
		}
		if f.End()-f.Start() < m.view.bpp { // thinner than a column, so not worth a mark
			continue
		}
		if f.Kind == "gene" && coding[f.Name] {
			continue // the CDS says the same thing
		}
		bar := geneBar{f, max(m.column(f.Start()), 0), min(m.column(f.End()-1)+1, columns)}
		if bar.to <= bar.from {
			bar.to = bar.from + 1
		}
		placed := false
		for r := range rows {
			last := rows[r][len(rows[r])-1]
			if bar.from > last.to+1 {
				rows[r], placed = append(rows[r], bar), true
				break
			}
		}
		if !placed && len(rows) < 3 {
			rows = append(rows, []geneBar{bar})
		}
	}
	return rows
}

// codingTrack returns the features in view whose protein is drawn: every CDS,
// and any other feature that carries a translation of its own.
func (m *model) codingTrack(columns int) []Feature {
	if m.ref == nil {
		return nil
	}
	contig := m.contigs()[m.contig].Name
	end := m.view.start + columns*m.view.bpp
	var out []Feature
	for _, f := range m.ref.featuresIn(contig, m.view.start, end) {
		if (f.Kind == "CDS" || f.Protein != "") && f.End() > m.view.start && f.Start() < end {
			out = append(out, f)
		}
	}
	return out
}

func (m *model) drawTracks(g grid, l layout) {
	back := cellStyle{bg: colorBackdrop}
	label := cellStyle{fg: colorDim, bg: colorBackdrop}
	row := trackTop
	contig := m.contigs()[m.contig].Name

	// ruler
	g.draw(row, gutterWidth-2-len("position"), "position", label)
	step := rulerStep(m.view.bpp)
	for pos := (m.view.start/step + 1) * step; pos < m.view.start+l.columns*m.view.bpp; pos += step {
		g.draw(row, gutterWidth+m.column(pos), "│"+comma(pos+1), cellStyle{fg: colorMuted, bg: colorBackdrop})
	}
	row++

	for _, bars := range m.geneTrack(l.columns) {
		g.draw(row, gutterWidth-2-len("genes"), "genes", label)
		for _, b := range bars {
			style := cellStyle{fg: "#7dcfff", bg: colorBackdrop}
			if b.feature.Kind != "CDS" {
				style.fg = "#9aa5b8"
			}
			g.fill(row, gutterWidth+b.from, gutterWidth+b.to, '─', style)
			if b.feature.Reverse {
				g.draw(row, gutterWidth+b.from, "◀", style)
			} else {
				g.draw(row, gutterWidth+b.to-1, "▶", style)
			}
			name := geneName(b.feature)
			if width := b.to - b.from; len(name)+2 <= width {
				at := gutterWidth + b.from + (width-len(name)-2)/2
				g.draw(row, at, " "+name+" ", cellStyle{fg: colorBright, bg: colorBackdrop, bold: true})
			}
		}
		row++
	}

	// variant markers, from the VCF where we have one and from the reads otherwise
	g.draw(row, gutterWidth-2-len("variants"), "variants", label)
	g.fill(row, gutterWidth, gutterWidth+l.columns, ' ', back)
	for pos, s := range m.sites {
		if col := m.column(pos); col >= 0 && col < l.columns {
			shade := colorMuted // 5% and up
			switch {
			case s.Ratio >= 0.5:
				shade = colorAlert
			case s.Ratio >= 0.3:
				shade = colorWarn
			}
			g.draw(row, gutterWidth+col, "▪", cellStyle{fg: shade, bg: colorBackdrop})
		}
	}
	for i, v := range m.variants {
		if v.Contig != contig {
			continue
		}
		if col := m.column(v.Start); col >= 0 && col < l.columns {
			style := cellStyle{fg: impactColor(v.Impact), bg: colorBackdrop, bold: true}
			if i == m.variantIndex {
				style.bg = colorLine
			}
			g.draw(row, gutterWidth+col, "▼", style)
		}
	}
	row++

	// reference sequence
	g.draw(row, gutterWidth-2-len("reference"), "reference", label)
	if seq := m.ref.sequence(contig); m.view.bpp == 1 && seq != nil {
		for col, spec := range m.cols {
			if col >= l.columns || spec.pos >= len(seq) {
				break
			}
			if spec.insert > 0 { // the reference has nothing here
				g.set(row, gutterWidth+col, '-', cellStyle{fg: colorLine, bg: colorBackdrop})
				continue
			}
			base := upper(seq[spec.pos])
			g.set(row, gutterWidth+col, rune(base), cellStyle{fg: chosenPalette.bases[base], bg: colorBackdrop, bold: true})
		}
	} else {
		g.fill(row, gutterWidth, gutterWidth+l.columns, '─', cellStyle{fg: colorLine, bg: colorBackdrop})
	}
	row++

	// amino acids of the coding sequences in view
	if l.aminoRow == 1 {
		g.draw(row, gutterWidth-2-len("protein"), "protein", label)
		for _, f := range m.codingTrack(l.columns) {
			m.drawAminoAcids(g, row, l, f)
		}
		row++
	}

	// coverage
	g.draw(row, gutterWidth-2-len("depth"), "depth", label)
	if m.covPeak > 0 {
		depth := shorten(fmt.Sprintf("depth %d", m.covPeak), 9)
		if m.truncated {
			depth = "sample"
		}
		g.draw(row, gutterWidth-2-len(depth), depth, label)
	}
	blocks := []rune("▁▂▃▄▅▆▇█")
	for col, n := range m.cov {
		if col >= l.columns {
			break
		}
		if n == 0 {
			g.set(row, gutterWidth+col, '·', cellStyle{fg: colorLine, bg: colorBackdrop})
			continue
		}
		level := min((n*len(blocks)-1)/max(m.covPeak, 1), len(blocks)-1)
		shade := mix("#2f4858", "#7dcfff", float64(n)/float64(max(m.covPeak, 1)))
		g.set(row, gutterWidth+col, blocks[level], cellStyle{fg: shade, bg: colorBackdrop})
	}
}

// markedPosition is the base the mark and the position line describe: the one
// under the mouse, or the selected variant when the mouse is elsewhere. It is
// -1 when a column is several bases wide, where there is no one base to point
// at.
func (m *model) markedPosition() int {
	if m.view.bpp != 1 {
		return -1
	}
	if m.hoverPos >= 0 {
		return min(m.hoverPos, m.contigs()[m.contig].Length-1)
	}
	if m.marked >= 0 {
		return min(m.marked, m.contigs()[m.contig].Length-1)
	}
	if v := m.selectedVariant(); v != nil && v.Contig == m.contigs()[m.contig].Name {
		return v.Start
	}
	return -1
}

// markColumn lights up one column of bases: the reference base and the read
// bases under it, at the position the mouse is over or at the selected variant
// when the mouse is elsewhere. Empty cells are left alone, and nothing is
// marked once a column holds more than one base.
func (m *model) markColumn(g grid, l layout) {
	from := m.markedPosition()
	if from < 0 {
		return
	}
	to := from + 1
	if m.hoverPos < 0 && m.marked < 0 { // the mark came from a variant, so its whole span lights up
		to = m.selectedVariant().End()
	}
	light := func(row, col int) {
		if c := g[row][gutterWidth+col]; c.r != ' ' {
			c.s.bg, c.s.bold = colorLine, true
			g[row][gutterWidth+col] = c
		}
	}
	for pos := from; pos < to; pos++ {
		col := m.column(pos)
		if col < 0 || col >= l.columns {
			continue
		}
		light(l.refRow, col)
		for row := l.pileupTop; row < l.pileupTop+l.pileupRows; row++ {
			light(row, col)
		}
	}
}

func (m *model) drawAminoAcids(g grid, row int, l layout, f Feature) {
	contig := m.contigs()[m.contig].Name
	positions := codingPositions(f)
	protein := m.ref.givenProtein(f, contig, positions)
	if protein == "" && f.Kind != "CDS" {
		return // a feature with no translation of its own, and none to work out
	}
	for i := 0; i+2 < len(positions); i += 3 {
		col := m.column(positions[i+1])
		if col < 0 || col >= l.columns {
			continue
		}
		var amino byte
		switch {
		case i/3 < len(protein):
			amino = protein[i/3]
		case f.Kind == "CDS": // a CDS keeps going past a translation that leaves the stop out
			amino = aminoAt(f.Table, codonAt(m.ref, f, contig, positions[i:]), i/3+1)
		default:
			continue
		}
		style := cellStyle{fg: colorMuted, bg: colorBackdrop}
		if (i/3)%2 == 0 {
			style.bg = "#191d27"
		}
		switch amino {
		case '*':
			style.fg, style.bold = chosenPalette.alert, true
		case 'M':
			style.fg, style.bold = chosenPalette.good, true
		}
		for k := 0; k < 3; k++ {
			if c := m.column(positions[i+k]); c >= 0 && c < l.columns {
				g.draw(row, gutterWidth+c, " ", style)
			}
		}
		g.draw(row, gutterWidth+col, string(amino), style)
	}
}

func rulerStep(bpp int) int {
	for _, step := range []int{10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, 100000, 500000, 1000000} {
		if step/bpp >= 12 {
			return step
		}
	}
	return 5000000
}

// pileup -------------------------------------------------------------------

func (m *model) drawPileup(g grid, l layout) {
	contig := m.contigs()[m.contig].Name
	for _, a := range m.reads {
		row := l.pileupTop + a.Row - m.rowOffset
		if a.Row < m.rowOffset || row >= l.pileupTop+l.pileupRows {
			continue
		}
		m.drawRead(g, row, l, a, contig)
	}
	if !m.showReads() && m.contigs()[m.contig].Profile == nil {
		g.draw(l.pileupTop+1, gutterWidth+2, "Measuring coverage across the contig…", cellStyle{fg: colorMuted, bg: colorBackdrop})
	} else if !m.showReads() {
		g.draw(l.pileupTop+1, gutterWidth+2, "Zoomed out: the depth track covers the whole contig.", cellStyle{fg: colorMuted, bg: colorBackdrop})
		g.draw(l.pileupTop+2, gutterWidth+2, "Press + to zoom to reads.", cellStyle{fg: colorDim, bg: colorBackdrop})
	} else if len(m.reads) == 0 {
		g.draw(l.pileupTop+1, gutterWidth+2, "No reads here that pass the filters.", cellStyle{fg: colorDim, bg: colorBackdrop})
	}
	if hidden := m.rowCount - m.rowOffset - l.pileupRows; hidden > 0 {
		note := fmt.Sprintf(" %d more rows below ", hidden)
		g.draw(l.pileupTop+l.pileupRows-1, m.width-l.panel-len(note)-1, note, cellStyle{fg: colorDim, bg: colorPanel})
	}
}

func (m *model) drawRead(g grid, row int, l layout, a *Alignment, contig string) {
	paired, flash := a.PairedWith(m.hover), a.PairedWith(m.flash)
	first, last := m.column(a.Start), m.column(a.End-1)+m.inserts[a.End-1]
	for col := max(first, 0); col <= last && col < min(l.columns, len(m.cols)); col++ {
		spec := m.cols[col]
		style := cellStyle{fg: readColor(m.color, a, spec.pos, m.readLength), bg: colorBackdrop, faint: a.SoftClipped(spec.pos)}
		if paired {
			style.bg = "#242a38"
		}
		if flash {
			style.bg, style.bold = colorAccent, true
		}
		r := ' '
		switch {
		case spec.insert > 0:
			r, style.fg = '-', colorLine
			if seq := a.Insertions[spec.pos]; spec.insert <= len(seq) {
				r, style.fg, style.bg = rune(seq[spec.insert-1]), colorBright, "#4a2f6d"
			}
		case m.view.bpp == 1:
			base := a.BaseAt(spec.pos)
			if base == 0 {
				continue
			}
			r = rune(base)
			style.fg, style.bold = m.baseColor(a, spec.pos, base, m.ref.Base(contig, spec.pos))
			if base == '-' {
				style.fg = chosenPalette.bases['-']
			}
		default:
			r = '─'
			if at := m.deviationInColumn(a, spec.pos); at >= 0 {
				r = '━'
				style.fg, style.bold = m.baseColor(a, at, a.BaseAt(at), m.ref.Base(contig, at))
			}
			if col == first && a.Reverse() {
				r = '◀'
			}
			if col == last && !a.Reverse() {
				r = '▶'
			}
		}
		g.set(row, gutterWidth+col, r, style)
	}
}

// deviationInColumn returns the first position in a column where a read
// differs from what the coloring compares against, or -1. The places were
// found as the region loaded: comparing every base again on every frame is
// what made highlighting lag when zoomed out.
func (m *model) deviationInColumn(a *Alignment, pos int) int {
	if m.highlight < 0 {
		return -1 // nothing is being picked out
	}
	places := a.MismatchAt
	if m.fromConsensus() {
		if !a.consensusDone { // worked out the first time the read is drawn, not for the whole screen
			a.ConsensusAt, a.consensusDone = consensusDeviations(a, m.consensus, m.view.start), true
		}
		places = a.ConsensusAt
	}
	if i := sort.SearchInts(places, pos); i < len(places) && places[i] < pos+m.view.bpp {
		return places[i]
	}
	return -1
}

// bars ---------------------------------------------------------------------

// bar is one labelled strip along the bottom of the screen.
type bar struct {
	label string
	tint  string // the color of its label, and of the tick down its left edge
	lines []string
}

// barIndent is how far the text of a strip sits from the left edge: past its
// label, and lined up with the reads above.
const barIndent = gutterWidth

// refreshBars works out the strips below the pileup: the base being pointed
// at, the read under the pointer, and how the reads are drawn. A strip with
// nothing to say is left out.
func (m *model) refreshBars() {
	gap := strings.Repeat(" ", clamp(2+m.width/100, 2, 6)) // a wider window, a little more air
	m.bars = nil
	add := func(label, tint string, fields []string) {
		if len(fields) > 0 {
			m.bars = append(m.bars, bar{label, tint, wrapFields(fields, m.width-barIndent-2, gap)})
		}
	}
	position := m.positionSummary()
	if position == nil { // the strip stays, so that the reads above it do not shift
		position = []string{""}
	}
	add("Position", colorAccent, position)
	if m.hover != nil {
		add("Read", "#7dcfff", m.readSummary(m.hover))
	}
	add("Filters", colorMuted, m.filterSummary())
}

// drawBars draws the strips below the pileup, each labelled along its left
// edge, and the keys along the very bottom.
func (m *model) drawBars(g grid, l layout) {
	row := m.height - l.barRows + 1 // the row above the strips stays empty
	for i, b := range m.bars {
		shade := []string{colorPanel, colorStripe}[i%2]
		for j, line := range b.lines {
			g.fill(row+j, 0, m.width, ' ', cellStyle{bg: shade})
			g.draw(row+j, barIndent, line, cellStyle{fg: colorText, bg: shade})
		}
		g.draw(row, 1, " "+b.label+" ", cellStyle{fg: b.tint, bg: colorLine, bold: true})
		row += len(b.lines)
	}
	keys := " ←→ pan   t tab   ↑↓ rows   +/- zoom   m mouse   ? help   q quit"
	g.draw(m.height-1, 0, keys, cellStyle{fg: colorDim, bg: colorBackdrop})
}

// filterSummary says how the reads are drawn and which of them are, value by
// value.
func (m *model) filterSummary() []string {
	drawn := fmt.Sprintf("Visualizing %d reads", len(m.reads))
	switch {
	case !m.showReads():
		drawn = "Zoom in to view reads."
	case m.filter != defaultFilter() && m.readsSeen > 0: // once a filter is on, say what it left out
		drawn = fmt.Sprintf("Visualizing %d of %d (%d%%) reads",
			len(m.reads), m.readsSeen, 100*len(m.reads)/m.readsSeen)
	}
	if m.truncated {
		drawn += ", a sample of this region"
	}
	fields := []string{m.coloringSummary(), drawn,
		fmt.Sprintf("Mapping quality ≥ %d", m.filter.MinMapQ), "Mismatches " + mismatchLimit(m.filter)}
	if m.filter.HideUnpaired {
		fields = append(fields, "Proper Pairs")
	}
	if m.filter.MinLength > 0 {
		fields = append(fields, "Read length ≥ "+comma(m.filter.MinLength))
	}
	if m.filter.MinMeanQual > 0 {
		fields = append(fields, fmt.Sprintf("Average quality ≥ %d", m.filter.MinMeanQual))
	}
	return fields
}

// wrapFields spaces values out along a line, carrying on to a second line when
// they do not fit, and breaking between values rather than inside them.
func wrapFields(fields []string, width int, gap string) []string {
	lines := []string{""}
	for _, field := range fields {
		last := len(lines) - 1
		switch {
		case lines[last] == "":
			lines[last] = field
		case len([]rune(lines[last]))+len(gap)+len([]rune(field)) <= width:
			lines[last] += gap + field
		case last == 0:
			lines = append(lines, "  "+field)
		default:
			lines[last] += gap + field
		}
	}
	for i, line := range lines {
		lines[i] = shorten(line, width)
	}
	return lines
}

// positionSummary describes the marked base, value by value: what the reads
// show there, and what a difference from the reference means for the protein.
func (m *model) positionSummary() []string {
	pos := m.markedPosition()
	if pos < 0 {
		return nil
	}
	contig := m.contigs()[m.contig].Name
	counts := tally(m.reads, pos)
	values := []string{contig + ":" + comma(pos+1)}
	if counts.Depth == 0 {
		return append(values, "no reads here")
	}
	share := func(base byte) string {
		return fmt.Sprintf("%c (%.0f%%)", base, 100*float64(counts.Count(base))/float64(counts.Depth))
	}
	ranked := counts.ranked()
	var alleles []string
	for _, base := range ranked {
		if counts.Count(base) > 0 {
			alleles = append(alleles, share(base))
		}
	}
	values = append(values, fmt.Sprintf("Depth: %dx", counts.Depth), "Alleles: "+strings.Join(alleles, " "))

	v := m.variantAt(pos)
	if v != nil && m.hover == nil { // the VCF's own call, which a pileup cannot show for an indel
		values = append(values, v.Label()+" "+v.Impact.String())
	} else if c := m.changeAt(pos); c.to > 0 {
		values = append(values, m.spell(c))
	}
	if v != nil {
		values = append(values, v.Details()...)
	}
	return values
}

// change is one difference from the reference: the bases from..to give way to
// alt. The zero value means no difference at all.
type change struct {
	from, to int
	alt      string
}

// changeAt is the difference the position line describes: the one in the read
// being pointed at, and otherwise the one most of the reads agree on. Taking
// the commonest keeps the line steady as the pointer moves along a deletion
// that a stray read spells differently.
func (m *model) changeAt(pos int) change {
	if m.hover != nil {
		if c := m.changeIn(m.hover, pos); c.to > 0 {
			return c
		}
	}
	votes, best := map[change]int{}, change{}
	for _, a := range m.reads {
		c := m.changeIn(a, pos)
		if c.to == 0 {
			continue // this read agrees with the reference
		}
		if votes[c]++; votes[c] > votes[best] {
			best = c
		}
	}
	return best
}

// changeIn is what one read shows at a position: a different base, an
// insertion, or the whole deletion it falls inside.
func (m *model) changeIn(a *Alignment, pos int) change {
	contig := m.contigs()[m.contig].Name
	refBase, base := m.ref.Base(contig, pos), a.BaseAt(pos)
	switch {
	case refBase == 0: // without a reference there is nothing to compare against
	case a.Insertions[pos] != "":
		return change{pos, pos + 1, string(refBase) + a.Insertions[pos]}
	case base == '-':
		from, to := a.deletion(pos)
		return change{from - 1, to, m.ref.Span(contig, from-1, from)}
	case base != 0 && base != ' ' && base != refBase:
		return change{pos, pos + 1, string(base)}
	}
	return change{}
}

// spell writes a change out with what it does to the protein.
func (m *model) spell(c change) string {
	contig := m.contigs()[m.contig].Name
	ref := m.ref.Span(contig, c.from, c.to)
	return label(ref, c.alt) + " " + m.ref.CodingImpact(contig, c.from, ref, c.alt).String()
}

func (m *model) variantAt(pos int) *Variant {
	for i, v := range m.variants {
		if v.Contig == m.contigs()[m.contig].Name && v.Start == pos {
			return &m.variants[i]
		}
	}
	return nil
}

// readSummary describes the read being pointed at, value by value.
func (m *model) readSummary(a *Alignment) []string {
	strand := "forward"
	if a.Reverse() {
		strand = "reverse"
	}
	mismatches := fmt.Sprintf("%d mismatches", a.Mismatches)
	if a.Mismatches == 1 {
		mismatches = "1 mismatch"
	}
	parts := []string{
		a.Name,
		fmt.Sprintf("%s–%s (%s bp)", comma(a.Start+1), comma(a.End), comma(a.End-a.Start)),
		strand,
		fmt.Sprintf("mapping quality %d", a.MapQ),
		fmt.Sprintf("average base quality %d", a.MeanQual),
		mismatches,
	}
	for pos, seq := range a.Insertions {
		parts = append(parts, fmt.Sprintf("insertion %s at %s", seq, comma(pos+1)))
		break
	}
	switch {
	case !a.Paired():
		parts = append(parts, "unpaired")
	case a.MateContig != "" && a.MateContig != m.contigs()[m.contig].Name:
		parts = append(parts, "paired, mate on "+a.MateContig)
	case a.InsertSize != 0:
		parts = append(parts, fmt.Sprintf("paired, pair distance: %s bp", comma(abs(a.InsertSize))))
	default:
		parts = append(parts, "paired")
	}
	return parts
}

func impactColor(i Impact) string {
	switch i.Kind {
	case "synonymous":
		return chosenPalette.good
	case "missense", "complex substitution":
		return chosenPalette.warn
	case "nonsense", "start loss", "stop loss", "frameshift insertion", "frameshift deletion":
		return chosenPalette.alert
	case "in-frame insertion", "in-frame deletion", "insertion", "deletion":
		return chosenPalette.other
	}
	return colorMuted
}

func shorten(s string, width int) string {
	r := []rune(s)
	if width <= 1 || len(r) <= width {
		return s
	}
	return string(r[:width-1]) + "…"
}

func comma(n int) string {
	s := fmt.Sprint(n)
	for i := len(s) - 3; i > 0; i -= 3 {
		s = s[:i] + "," + s[i:]
	}
	return s
}
