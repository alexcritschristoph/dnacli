package main

import (
	"encoding/base64"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	tea "github.com/charmbracelet/bubbletea"
)

// press sends a key the way the terminal would: named keys have their own type,
// everything else arrives as runes.
func press(m *model, key string) {
	named := map[string]tea.KeyType{
		"left": tea.KeyLeft, "right": tea.KeyRight, "up": tea.KeyUp, "down": tea.KeyDown,
		"home": tea.KeyHome, "end": tea.KeyEnd, "enter": tea.KeyEnter, "esc": tea.KeyEsc,
		"tab": tea.KeyTab, "backspace": tea.KeyBackspace, "pgup": tea.KeyPgUp, "pgdown": tea.KeyPgDown,
		" ": tea.KeySpace,
	}
	msg := tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune(key)}
	if kind, ok := named[key]; ok {
		msg = tea.KeyMsg{Type: kind}
	}
	m.handleKey(msg)
}

// run types a command into the bar and presses enter, the way a user would.
func run(m *model, text string) {
	for _, r := range text {
		press(m, string(r))
	}
	press(m, "enter")
}

func mouse(m *model, x, y int, action tea.MouseAction, button tea.MouseButton) {
	m.handleMouse(tea.MouseMsg{X: x, Y: y, Action: action, Button: button})
}

func ready(t *testing.T, stem string) *model {
	m := testModel(t, stem)
	m.width, m.height = 150, 44
	m.reload()
	return m
}

func TestScrollingStopsAtBothEndsOfTheContig(t *testing.T) {
	m := ready(t, "viral_variants")
	for i := 0; i < 40; i++ {
		press(m, "left")
	}
	if m.view.start != 0 {
		t.Fatalf("panning left ended at %d, want 0", m.view.start)
	}
	length := m.contigs()[m.contig].Length
	for i := 0; i < 1000; i++ {
		press(m, "right")
	}
	if want := length - m.span(); m.view.start != want {
		t.Fatalf("panning right ended at %d, want %d", m.view.start, want)
	}
	press(m, "home")
	if m.view.start != 0 {
		t.Fatalf("home went to %d, want 0", m.view.start)
	}
	press(m, "end")
	if m.view.start+m.span() != length {
		t.Fatalf("end left %d bases past the contig", m.view.start+m.span()-length)
	}
}

func TestZoomStaysCentredAndStopsAtWholeContig(t *testing.T) {
	m := ready(t, "viral_variants")
	m.goTo(4000)
	centre := m.view.start + m.span()/2
	for i := 0; i < 3; i++ {
		press(m, "-")
	}
	if m.view.bpp == 1 {
		t.Fatal("zooming out did not change the scale")
	}
	if off := abs(m.view.start + m.span()/2 - centre); off > m.view.bpp*2 {
		t.Fatalf("zoom moved the centre by %d bases", off)
	}
	for i := 0; i < 30; i++ {
		press(m, "-")
	}
	if m.span() > m.contigs()[m.contig].Length+m.layout().columns {
		t.Fatalf("zoomed out past the whole contig: %d bases on screen", m.span())
	}
	for i := 0; i < 30; i++ {
		press(m, "+")
	}
	if m.view.bpp != 1 {
		t.Fatalf("zooming all the way in gave %d bases per column", m.view.bpp)
	}
}

func TestPairsShareARowAndHoverLightsBoth(t *testing.T) {
	m := ready(t, "paired_reads")
	m.view = viewState{0, 10} // wide enough for both mates to be on screen
	m.reload()
	mates := map[string][]*Alignment{}
	for _, a := range m.reads {
		mates[a.Name] = append(mates[a.Name], a)
	}
	pairs := 0
	for name, both := range mates {
		if len(both) != 2 {
			continue
		}
		pairs++
		if both[0].Row != both[1].Row && both[1].Start-both[0].End < m.span() {
			t.Fatalf("the mates of %s landed on rows %d and %d", name, both[0].Row, both[1].Row)
		}
	}
	if pairs < 5 {
		t.Fatalf("only %d pairs are on screen", pairs)
	}
	if m.rowCount > len(m.reads)/2 {
		t.Fatalf("%d reads needed %d rows; pairs are not sharing rows", len(m.reads), m.rowCount)
	}

	read := m.reads[0]
	l := m.layout()
	mouse(m, gutterWidth+max(m.column(read.Start), 0), l.pileupTop+read.Row-m.rowOffset, tea.MouseActionMotion, tea.MouseButtonNone)
	if m.hover != read {
		t.Fatalf("hovering the read gave %v", m.hover)
	}
	lit := 0
	for _, a := range m.reads {
		if a.PairedWith(m.hover) {
			lit++
		}
	}
	if lit != 2 {
		t.Fatalf("hovering lit %d reads, want the read and its mate", lit)
	}
	if summary := strings.Join(m.readSummary(m.hover), " "); !strings.Contains(summary, "paired, pair distance:") {
		t.Fatalf("the read details say nothing about the pair: %s", summary)
	}
}

// copied returns what a mouse press sends to the terminal, which is where the
// OSC 52 escape that carries the clipboard goes.
func copied(t *testing.T, m *model, x, y int, button tea.MouseButton) string {
	t.Helper()
	terminal, w, _ := os.Pipe()
	stderr := os.Stderr
	os.Stderr = w
	mouse(m, x, y, tea.MouseActionPress, button)
	w.Close()
	os.Stderr = stderr
	sent, _ := io.ReadAll(terminal)
	_, payload, _ := strings.Cut(string(sent), ";c;")
	text, err := base64.StdEncoding.DecodeString(strings.TrimSuffix(payload, "\a"))
	if err != nil {
		t.Fatalf("the terminal was sent %q", sent)
	}
	return string(text)
}

func TestClickingAReadCopiesItAsFasta(t *testing.T) {
	m := ready(t, "viral_variants")
	read := m.reads[0]
	before := m.view
	text := copied(t, m, gutterWidth+max(m.column(read.Start), 0), m.layout().pileupTop+read.Row, tea.MouseButtonLeft)
	if want := ">" + read.Name + "\n" + read.Sequence; text != want {
		t.Fatalf("clicking a read copied %q, want %q", text, want)
	}
	if len(read.Sequence) < 50 || strings.ContainsAny(read.Sequence, "- ") {
		t.Fatalf("the copied sequence is not a plain read: %q", read.Sequence)
	}
	if m.view != before {
		t.Fatalf("clicking a read moved the view to %v", m.view)
	}
	if !strings.Contains(m.message, "read "+read.Name+" copied") {
		t.Fatalf("clicking a read said %q", m.message)
	}
}

func TestClickingAReadFlashesIt(t *testing.T) {
	m := ready(t, "viral_variants")
	read := m.reads[0]
	x, y := gutterWidth+max(m.column(read.Start), 0), m.layout().pileupTop+read.Row
	if text := copied(t, m, x, y, tea.MouseButtonRight); text != "" || m.flash != nil {
		t.Fatalf("right-clicking a read copied %q", text)
	}
	copied(t, m, x, y, tea.MouseButtonLeft)
	if m.flash != read {
		t.Fatal("the read did not light up")
	}
	m.View()
	if m.frame[y][x].s.bg != colorAccent {
		t.Fatal("the flash is not drawn")
	}

	m.flashUntil = time.Now().Add(-time.Second) // once the flash has had its moment
	m.Update(fadeMsg{})
	if m.flash != nil {
		t.Fatal("the flash stayed on")
	}
	if m.message == "" {
		t.Fatal("the message went away with the flash")
	}
}

func TestJumpingBetweenVariants(t *testing.T) {
	m := ready(t, "viral_variants")
	run(m, "variants")
	if m.panel != "variants" || !m.panelFocus {
		t.Fatalf("the variants command left the panel %q, focus %v", m.panel, m.panelFocus)
	}
	press(m, "down")
	first := m.selectedVariant()
	if first == nil || m.view.start > first.Start || m.view.start+m.span() < first.Start {
		t.Fatalf("the first variant is not on screen: %v", first)
	}
	press(m, "down")
	if second := m.selectedVariant(); second == nil || second.Start <= first.Start {
		t.Fatal("down did not move to a later variant")
	}
	press(m, "up")
	if m.selectedVariant().Start != first.Start {
		t.Fatal("up did not move back")
	}
	if !strings.Contains(strings.Join(m.positionSummary(), " "), "protein:") {
		t.Fatalf("the position line has no protein change: %s", m.positionSummary())
	}

	mouse(m, 10, m.layout().pileupTop, tea.MouseActionPress, tea.MouseButtonLeft)
	if m.panelFocus {
		t.Fatal("clicking the reads did not take the arrow keys back")
	}
	chosen := m.variantIndex
	press(m, "down")
	if m.variantIndex != chosen {
		t.Fatal("down still moved the variant after the reads were clicked")
	}
}

func TestGoToAcceptsPositionsAndRanges(t *testing.T) {
	m := ready(t, "multi_contig_uneven")
	m.parseGoto("CHR2:3,000")
	if m.contigs()[m.contig].Name != "CHR2" {
		t.Fatalf("landed on contig %s", m.contigs()[m.contig].Name)
	}
	if off := abs(m.view.start + m.span()/2 - 2999); off > m.view.bpp {
		t.Fatalf("centred on %d, want 2999", m.view.start+m.span()/2)
	}
	if m.markedPosition() != 2999 { // the spot gone to is marked, as a chosen variant is
		t.Fatalf("go left the mark at %d, want 2999", m.markedPosition())
	}
	m.parseGoto("CHR1:1000-2000")
	if m.span() < 1000 || m.span() > 1400 {
		t.Fatalf("a 1000 base range gave a %d base view", m.span())
	}
	m.parseGoto("NOPE:1")
	if !strings.Contains(m.message, "NOPE") {
		t.Fatalf("no complaint about an unknown contig: %q", m.message)
	}
}

func TestFilterCommandsChangeWhatIsDrawn(t *testing.T) {
	m := ready(t, "supertest")
	m.view = viewState{0, 20}
	m.reload()
	before := len(m.reads)
	run(m, "filters mapq 30")
	if m.filter.MinMapQ != 30 {
		t.Fatalf("mapping quality is %d, want 30", m.filter.MinMapQ)
	}
	if m.panel != "filters" {
		t.Fatal("a filter command should show the filters")
	}
	if len(m.reads) >= before {
		t.Fatalf("raising the mapping quality kept %d of %d reads", len(m.reads), before)
	}
	run(m, "filters reset")
	if len(m.reads) != before {
		t.Fatalf("reset gave %d reads, want %d", len(m.reads), before)
	}
	run(m, "filters mismatches 0")
	if m.filter.MaxMismatches != 0 {
		t.Fatalf("mismatch limit is %d, want 0", m.filter.MaxMismatches)
	}
	for _, a := range m.reads {
		if a.Mismatches > 0 {
			t.Fatalf("read %s has %d mismatches", a.Name, a.Mismatches)
		}
	}

	run(m, "filters reset") // the arrow keys change the chosen filter too
	for filterRows[m.filterIndex].name != "Mismatches at most" {
		press(m, "down")
	}
	press(m, "right")
	if m.filter.MaxMismatches != 0 {
		t.Fatalf("the arrows left filter %d at %v", m.filterIndex, m.filter)
	}
	for m.filterIndex > 0 {
		press(m, "up")
	}
	for i := 0; i < 6; i++ {
		press(m, "right")
	}
	if m.filter.MinMapQ != 30 {
		t.Fatalf("the arrows left the mapping quality at %d", m.filter.MinMapQ)
	}
}

func TestColorSchemesAndPanels(t *testing.T) {
	m := ready(t, "viral_variants")
	if m.color != -1 {
		t.Fatalf("a color scheme is on by default: %s", colors[m.color])
	}
	run(m, "colors")
	if m.panel != "colors" {
		t.Fatalf("panel is %q", m.panel)
	}
	press(m, "down")
	press(m, "down")
	press(m, "down")
	if m.color < 0 || colors[m.color] != "Mapping quality" {
		t.Fatalf("three presses down previewed the color %d", m.color)
	}
	if m.highlight < 0 || highlights[m.highlight] != "Mismatch (from reference)" {
		t.Fatal("moving through the colors cleared the SNP highlighting")
	}
	press(m, "enter") // the preview becomes the locked choice
	if !m.coloringLocked {
		t.Fatal("enter did not lock the previewed choice")
	}
	for i := range colors {
		m.color = i
		if color := readColor(i, m.reads[0], m.reads[0].Start, m.readLength); !strings.HasPrefix(color, "#") {
			t.Fatalf("scheme %s gave the color %q", colors[i], color)
		}
		screen(t, m, 150, 44) // every scheme has to draw without trouble
	}
	press(m, "esc")
	if m.panel != "" {
		t.Fatal("esc did not close the panel")
	}
}

func TestColoringPanelPreviewsThenLocks(t *testing.T) {
	m := ready(t, "viral_variants")
	run(m, "colors")
	press(m, "down") // before anything is locked the cursor and the highlight move together
	if m.color != namedIndex(colors, "Nucleotides") {
		t.Fatalf("moving the cursor previewed the color %d", m.color)
	}
	if m.highlight < 0 {
		t.Fatal("moving the cursor cleared the SNP highlighting")
	}
	press(m, "enter") // enter locks the previewed choice
	if !m.coloringLocked || m.color != namedIndex(colors, "Nucleotides") {
		t.Fatal("enter did not lock the previewed color")
	}
	press(m, "down") // now only the cursor moves
	if m.color != namedIndex(colors, "Nucleotides") {
		t.Fatal("moving the cursor changed the locked color")
	}
	if highlights[m.highlight] != "Mismatch (from reference)" {
		t.Fatal("moving through the colors changed the SNP highlighting")
	}
	press(m, "up")
	m.View() // both sections have a choice, and both have to show it
	for _, name := range []string{"Mismatch (from reference)", "Nucleotides"} {
		for i, line := range coloringPanel {
			if line.name != name {
				continue
			}
			if bg := m.frame[m.layout().listTop+i][m.width-3].s.bg; bg != colorLine {
				t.Fatalf("%s is not drawn as chosen: %q", name, bg)
			}
		}
	}
	press(m, "enter") // on the locked choice, enter resets it to nothing
	if m.color != -1 {
		t.Fatalf("enter left the locked color at %d", m.color)
	}
	if m.highlight < 0 {
		t.Fatal("clearing the color cleared the SNP highlighting too")
	}
	if color, loud := m.baseColor(m.reads[0], m.reads[0].Start, m.reads[0].BaseAt(m.reads[0].Start), 0); loud {
		t.Fatalf("with no color scheme a matching base still stands out in %s", color)
	}
	press(m, "down")
	if m.color != -1 {
		t.Fatal("moving the cursor after a reset chose a color")
	}
	press(m, "up")
	press(m, "enter") // and again brings it back
	if m.color != namedIndex(colors, "Nucleotides") {
		t.Fatalf("enter again gave the color %d", m.color)
	}
	press(m, "up") // onto the locked SNP highlighting
	press(m, "enter")
	if m.highlight != -1 {
		t.Fatalf("enter left the SNP highlighting at %d", m.highlight)
	}
	for _, a := range m.reads {
		base := a.BaseAt(a.Start)
		if _, loud := m.baseColor(a, a.Start, base, m.ref.Base(m.contigs()[m.contig].Name, a.Start)); loud {
			t.Fatal("with the SNP highlighting off, a difference still stands out")
		}
	}
	click := func(name string) {
		for i, line := range coloringPanel {
			if line.name == name {
				mouse(m, m.width-4, m.layout().listTop+i, tea.MouseActionPress, tea.MouseButtonLeft)
				return
			}
		}
		t.Fatalf("the panel has no row called %s", name)
	}
	click("Read strand")
	if m.color != namedIndex(colors, "Read strand") {
		t.Fatal("clicking a color did not choose it")
	}
	click("Read strand")
	if m.color != -1 {
		t.Fatal("clicking the locked color again did not reset it")
	}
	press(m, "esc")
	run(m, "colors") // opening the panel again starts previewing afresh
	if m.coloringLocked {
		t.Fatal("a freshly opened panel is still locked")
	}
	press(m, "up")
	if m.color != namedIndex(colors, "Nucleotides") {
		t.Fatalf("previewing after reopening gave the color %d", m.color)
	}
}

func TestMismatchColoringUsesTheReferenceOrTheConsensus(t *testing.T) {
	m := ready(t, "supertest")
	contig := m.contigs()[m.contig].Name
	if name := highlights[m.highlight]; name != "Mismatch (from reference)" {
		t.Fatalf("a reference is loaded, so the coloring should compare against it, not %s", name)
	}
	if m.color != -1 {
		t.Fatalf("a color scheme is on by default: %s", colors[m.color])
	}
	pos := 56 // the one place in this file where the reads and the reference disagree
	refBase := m.ref.Base(contig, pos)
	if refBase != 'T' || m.consensusBaseAt(pos) != 'A' {
		t.Fatalf("testdata changed: the reference gives %c, the consensus %c", refBase, m.consensusBaseAt(pos))
	}
	cellOf := func(base byte) cell {
		for _, a := range m.reads {
			if a.BaseAt(pos) == base {
				return m.frame[m.layout().pileupTop+a.Row-m.rowOffset][gutterWidth+m.column(pos)]
			}
		}
		t.Fatalf("no read carries %c at %d", base, pos)
		return cell{}
	}
	m.View()
	if c := cellOf('A'); !c.s.bold || c.s.fg != chosenPalette.bases['A'] {
		t.Fatalf("against the reference the differing A is drawn %q", c.s.fg)
	}
	if c := cellOf('T'); c.s.bold || c.s.fg != quietGrey {
		t.Fatalf("against the reference the matching T is drawn %q", c.s.fg)
	}
	m.setHighlight(namedIndex(highlights, "Mismatch (from consensus)"))
	m.View()
	if c := cellOf('A'); c.s.bold || c.s.fg != quietGrey {
		t.Fatalf("against the consensus the matching A is drawn %q", c.s.fg)
	}
	if c := cellOf('T'); !c.s.bold || c.s.fg != chosenPalette.bases['T'] {
		t.Fatalf("against the consensus the differing T is drawn %q", c.s.fg)
	}
}

func TestGeneSearchWithoutAReference(t *testing.T) {
	src, err := OpenAlignments("testdata/viral_variants.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	m.width, m.height = 150, 44
	m.reload()
	m.editCommand(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("gene ")})
	if list := m.suggestions(); len(list) != 0 {
		t.Fatalf("no reference is loaded, so gene suggested %v", list)
	}
	m.View()
}

func TestConsensusHighlightingWithoutAReference(t *testing.T) {
	src, err := OpenAlignments("testdata/viral_variants.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	if name := highlights[m.highlight]; name != "Mismatch (from consensus)" {
		t.Fatalf("no reference is loaded, so the coloring should compare against the consensus, not %s", name)
	}
	if m.color != -1 {
		t.Fatalf("a color scheme is on by default: %s", colors[m.color])
	}
	m.width, m.height = 150, 44
	m.reload()
	m.View()
	found := false
	for _, a := range m.reads {
		for pos := a.Start; pos < a.End; pos++ {
			base, consensus := a.BaseAt(pos), m.consensusBaseAt(pos)
			if consensus == 0 || baseIndex(base) < 0 || base == consensus {
				continue
			}
			found = true
			if c := m.frame[m.layout().pileupTop+a.Row-m.rowOffset][gutterWidth+m.column(pos)]; !c.s.bold || c.s.fg != chosenPalette.bases[base] {
				t.Fatalf("%c at %d differs from the consensus %c but is drawn %q", base, pos, consensus, c.s.fg)
			}
		}
	}
	if !found {
		t.Fatal("no read differed from the consensus, so nothing was tested")
	}

	m.view.bpp = 20 // zoomed out, where a column holds many bases
	m.reload()
	m.View()
	for _, a := range m.reads {
		if len(a.ConsensusAt) == 0 {
			continue
		}
		at := a.ConsensusAt[0]
		c := m.frame[m.layout().pileupTop+a.Row-m.rowOffset][gutterWidth+m.column(at)]
		if !c.s.bold || c.s.fg != chosenPalette.bases[a.BaseAt(at)] {
			t.Fatalf("%c at %d differs from the consensus but its column is drawn %q", a.BaseAt(at), at, c.s.fg)
		}
		return
	}
	t.Fatal("zoomed out, no read was marked as differing from the consensus")
}

func TestContigListShowsCoverageSketches(t *testing.T) {
	m := ready(t, "multi_contig_uneven")
	m.Update(sketchCoverage(m.src)())
	run(m, "contigs")
	press(m, "down")
	if m.contigs()[m.contig].Name != "CHR2" {
		t.Fatalf("opened contig %s", m.contigs()[m.contig].Name)
	}
	if m.view.start != 0 || m.view.bpp != 1 {
		t.Fatalf("a new contig opened at %v", m.view)
	}
	line := sparkline(m.contigs()[0].Profile, 25)
	if strings.Count(line, "·") > 5 {
		t.Fatalf("the coverage sketch is mostly empty: %q", line)
	}
}

func TestEveryTestFileDraws(t *testing.T) {
	for _, stem := range []string{"viral_variants", "plasmid_even", "mito_lowcov", "supertest", "multi_contig_uneven", "paired_reads"} {
		m := testModel(t, stem)
		for _, size := range [][2]int{{150, 44}, {80, 24}, {200, 60}} {
			for _, bpp := range []int{1, 10, 200} {
				m.view = viewState{0, bpp}
				lines := screen(t, m, size[0], size[1])
				if !strings.Contains(lines[0], "dnacli") {
					t.Fatalf("%s at %v: header is missing", stem, size)
				}
			}
		}
	}
}

// The header names the program, without a release number, at every width.
func TestHeaderNamesTheProgram(t *testing.T) {
	m := testModel(t, "viral_variants")
	for _, width := range []int{40, 80, 150} {
		lines := screen(t, m, width, 30)
		if !strings.Contains(lines[0], "dnacli") {
			t.Fatalf("at %d columns the header says %q", width, lines[0])
		}
		if strings.Contains(lines[0], version) {
			t.Fatalf("at %d columns the header shows the version: %q", width, lines[0])
		}
	}
}

func TestAminoRowUsesTheFilesTranslations(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{0, 1}
	l := m.layout()
	draw := func(f Feature) string {
		g := newGrid(m.width, 2)
		m.drawAminoAcids(g, 0, l, f)
		var out []rune
		for _, c := range g[0] {
			if c.r != ' ' {
				out = append(out, c.r)
			}
		}
		return string(out)
	}
	if got := draw(Feature{Kind: "CDS", Protein: "QQQ", Ranges: []Range{{0, 9}}}); got != "QQQ" {
		t.Fatalf("the CDS drew %q, want the file's QQQ", got)
	}
	if got := draw(Feature{Kind: "misc_feature", Protein: "WWW", Ranges: []Range{{0, 9}}}); got != "WWW" {
		t.Fatalf("a misc_feature with a translation drew %q", got)
	}
	if got := draw(Feature{Kind: "misc_feature", Ranges: []Range{{0, 9}}}); got != "" {
		t.Fatalf("a misc_feature without a translation drew %q", got)
	}
	plain := Feature{Kind: "CDS", Ranges: []Range{{0, 9}}}
	want := ""
	for i := 0; i < 9; i += 3 {
		codon := codonAt(m.ref, plain, "VIRUS1", []int{i, i + 1, i + 2})
		want += string(aminoAt(plain.Table, codon, i/3+1))
	}
	if got := draw(plain); got != want {
		t.Fatalf("a CDS without a translation drew %q, want the sequence's %q", got, want)
	}
}

func TestVariantBasesStandOutFromMatchingOnes(t *testing.T) {
	m := ready(t, "viral_variants")
	m.showVariant(0)
	v := m.variants[0]
	refBase := m.ref.Base(v.Contig, v.Start)
	carrying := 0
	for _, a := range m.reads {
		base := a.BaseAt(v.Start)
		color, loud := m.baseColor(a, v.Start, base, refBase)
		switch base {
		case 'C':
			carrying++
			if !loud || color != chosenPalette.bases['C'] {
				t.Fatalf("the variant base is drawn %s, loud=%v", color, loud)
			}
		case refBase:
			if loud {
				t.Fatal("a matching base is drawn as loudly as a variant")
			}
		}
	}
	if carrying == 0 {
		t.Fatal("no read carries the variant")
	}
	for i, s := range colors {
		if s == "Protein change" {
			m.color = i
		}
	}
	for _, a := range m.reads {
		if a.BaseAt(v.Start) == 'C' {
			if color, _ := m.baseColor(a, v.Start, 'C', refBase); color != impactColor(v.Impact) {
				t.Fatalf("protein coloring drew a %s change as %s", v.Impact.Kind, color)
			}
			break
		}
	}
}

func TestTheSelectedColumnHighlightsItsBases(t *testing.T) {
	m := ready(t, "viral_variants")
	l := m.layout()
	read := m.reads[0]
	col := max(m.column(read.Start), 0) + 2
	mouse(m, gutterWidth+col, l.pileupTop+read.Row, tea.MouseActionMotion, tea.MouseButtonNone)
	m.View()
	lit := func(row, col int) bool {
		return m.frame[row][gutterWidth+col].s.bg == colorLine
	}
	if !lit(l.refRow, col) {
		t.Fatal("the reference base in the column is not lit")
	}
	if !lit(l.pileupTop+read.Row, col) {
		t.Fatal("the read base in the column is not lit")
	}
	if lit(l.refRow, col+1) {
		t.Fatal("the next column along is lit too")
	}
	for _, row := range []int{trackTop, l.refRow - 1, l.pileupTop - 1} { // ruler, variants, depth
		if lit(row, col) {
			t.Fatalf("row %d is lit; only bases should be", row)
		}
	}
	if lit(l.pileupTop+l.pileupRows-1, col) {
		t.Fatal("an empty row below the reads is lit")
	}

	m.view.bpp = 10 // a column is no longer one base, so nothing should be marked
	m.reload()
	mouse(m, gutterWidth+col, l.pileupTop, tea.MouseActionMotion, tea.MouseButtonNone)
	m.View()
	for row := trackTop; row < l.pileupTop+l.pileupRows; row++ {
		if lit(row, col) {
			t.Fatalf("row %d is still marked at 10 bases per column", row)
		}
	}
}

func TestMessagesSitInTheCommandBarAndFade(t *testing.T) {
	m := ready(t, "supertest") // no VCF here, so n has something to complain about
	read := m.reads[0]
	mouse(m, gutterWidth+max(m.column(read.Start), 0), m.layout().pileupTop+read.Row, tea.MouseActionMotion, tea.MouseButtonNone)
	rows := m.layout().pileupRows
	run(m, "gene nothing")
	lines := strings.Split(ansi.ReplaceAllString(m.View(), ""), "\n")
	if !strings.Contains(lines[commandRow], "No gene") {
		t.Fatalf("the message is not in the command bar: %q", lines[commandRow])
	}
	if !strings.Contains(strings.Join(lines, "\n"), read.Name) {
		t.Fatalf("the read strip is gone: %q", lines)
	}
	if m.layout().pileupRows != rows {
		t.Fatal("the message took a row from the pileup")
	}

	m.messageUntil = time.Now().Add(-time.Second)
	m.Update(fadeMsg{})
	if m.message != "" {
		t.Fatalf("the message stayed on screen: %q", m.message)
	}
	press(m, "m")
	m.messageUntil = time.Now().Add(-time.Second)
	m.Update(fadeMsg{})
	lines = strings.Split(ansi.ReplaceAllString(m.View(), ""), "\n")
	if !strings.Contains(lines[commandRow], "mouse off") { // the standing notice has the spot back
		t.Fatalf("the mouse notice is missing: %q", lines[commandRow])
	}
}

func TestThePositionLineDescribesTheColumn(t *testing.T) {
	m := ready(t, "viral_variants")
	m.width = 200    // wide enough for the VCF's own fields at the end of the line
	m.showVariant(0) // A→C at 1,201, in orf1a
	line := strings.Join(m.positionSummary(), "  ")
	for _, want := range []string{"VIRUS1:1,201", "Depth: 17x", "Alleles: C (100%)",
		"A→C orf1a:902 protein: N301T missense", "quality: 100", "PASS"} {
		if !strings.Contains(line, want) {
			t.Fatalf("the position line is missing %q: %s", want, line)
		}
	}
	if silent := m.ref.CodingImpact("VIRUS1", 1201, "C", "T"); strings.Contains(silent.String(), "protein") {
		t.Fatalf("a synonymous change should not spell out a protein change: %s", silent)
	}

	m.width = 80 // too narrow for one line, so it wraps between fields
	m.refreshBars()
	wrapped := m.bars[0].lines
	if len(wrapped) != 2 || strings.Contains(wrapped[0], "…") {
		t.Fatalf("the position strip did not wrap: %q", wrapped)
	}
	bare := func(s string) string { return strings.ReplaceAll(s, " ", "") }
	if bare(strings.Join(wrapped, "")) != bare(line) {
		t.Fatalf("wrapping changed the text: %q", wrapped)
	}

	m.view.bpp = 10 // no single base to describe out here
	m.reload()
	if m.positionSummary() != nil {
		t.Fatalf("the position line still describes a base at 10 bases per column: %s", m.positionSummary())
	}
}

func TestTheReadLineOnlyAppearsWhilePointingAtARead(t *testing.T) {
	m := ready(t, "paired_reads")
	m.refreshBars()
	rows, strips := m.layout().pileupRows, len(m.bars)
	read := m.reads[0]
	mouse(m, gutterWidth+max(m.column(read.Start), 0), m.layout().pileupTop+read.Row, tea.MouseActionMotion, tea.MouseButtonNone)
	if len(m.bars) != strips+1 || m.layout().pileupRows != rows-1 {
		t.Fatal("pointing at a read did not open a strip for it")
	}
	line := strings.Join(m.readSummary(read), " ")
	if !strings.Contains(line, "paired, pair distance: 960 bp") || strings.Contains(line, read.Cigar) {
		t.Fatalf("the read line should carry the distance between the mates and no CIGAR: %s", line)
	}
	for _, a := range m.reads {
		if !a.Paired() {
			if !strings.Contains(strings.Join(m.readSummary(a), " "), "unpaired") {
				t.Fatalf("a single-end read is not called unpaired: %s", m.readSummary(a))
			}
			break
		}
	}
}

func TestHoveringAnIndelSaysSo(t *testing.T) {
	m := ready(t, "viral_variants")
	deletion, insertion := m.variants[2], m.variants[3] // AAAAA…→A at 3,100 and C→CTTT at 4,300
	m.showVariant(2)
	hover := func(pos int, a *Alignment) string {
		mouse(m, gutterWidth+m.column(pos), m.layout().pileupTop+a.Row, tea.MouseActionMotion, tea.MouseButtonNone)
		return strings.Join(m.positionSummary(), "  ")
	}
	for _, a := range m.reads {
		if a.BaseAt(deletion.Start+1) != '-' {
			continue
		}
		line := hover(deletion.Start+1, a)
		if !strings.Contains(line, "deletion") || !strings.Contains(line, "→A") {
			t.Fatalf("pointing at a deleted base does not call it a deletion: %s", line)
		}
		if !strings.Contains(line, "Alleles: - (") {
			t.Fatalf("deleted bases are missing from the alleles: %s", line)
		}
		break
	}
	for pos := deletion.Start + 1; pos < deletion.Start+4; pos++ { // one call for the whole run
		if m.spell(m.changeAt(pos)) != m.spell(m.changeAt(deletion.Start+1)) {
			t.Fatalf("the same deletion is called two things: %s at %d", m.spell(m.changeAt(pos)), pos+1)
		}
	}

	m.showVariant(3)
	for _, a := range m.reads {
		if a.Insertions[insertion.Start] == "" {
			continue
		}
		if line := hover(insertion.Start, a); !strings.Contains(line, "C→CTTT") || !strings.Contains(line, "insertion") {
			t.Fatalf("pointing at an inserted base does not call it an insertion: %s", line)
		}
		return
	}
	t.Fatal("no read carries the insertion")
}

func TestInsertionsGetColumnsOfTheirOwn(t *testing.T) {
	m := ready(t, "supertest")
	m.parseGoto("SUPER_A:1480-1530")
	m.view.bpp = 1
	m.reload()
	at, width := 1496, 3 // three bases inserted after this position, in some of the reads
	if m.inserts[at] != width {
		t.Fatalf("the insertion at %d got %d columns, want %d", at+1, m.inserts[at], width)
	}
	if m.column(at+1) != m.column(at)+width+1 {
		t.Fatal("the bases after the insertion do not line up again")
	}
	m.View()
	shown := func(row int) string {
		out := ""
		for col := m.column(at) + 1; col <= m.column(at)+width; col++ {
			out += string(m.frame[row][gutterWidth+col].r)
		}
		return out
	}
	l := m.layout()
	if got := shown(l.refRow); got != "---" {
		t.Fatalf("the reference should be a gap across the insertion, not %q", got)
	}
	with, without := "", ""
	for _, a := range m.reads {
		row := l.pileupTop + a.Row - m.rowOffset
		if a.Start >= at || a.End <= at+1 || row >= l.pileupTop+l.pileupRows {
			continue
		}
		if a.Insertions[at] != "" {
			with = shown(row)
		} else {
			without = shown(row)
		}
	}
	if with != "GGC" {
		t.Fatalf("a read carrying the insertion shows %q, want its inserted bases", with)
	}
	if without != "---" {
		t.Fatalf("a read without the insertion shows %q, want a gap", without)
	}
}

func TestSoftClipsCanBeTurnedOff(t *testing.T) {
	m := ready(t, "paired_reads")
	if !m.filter.SoftClips {
		t.Fatal("soft clips should be shown to begin with")
	}
	clipped := func() *Alignment {
		for _, a := range m.reads {
			if a.Name == "pair000" {
				return a
			}
		}
		return nil
	}
	on := clipped()
	if on.LeadClip != 8 || on.Start != 12 || on.Mismatches != 0 {
		t.Fatalf("the clipped read came out as %d-%d, %d clipped, %d mismatches", on.Start, on.End, on.LeadClip, on.Mismatches)
	}
	if !on.SoftClipped(on.Start) || on.SoftClipped(on.Start+on.LeadClip) {
		t.Fatal("the clipped bases are not marked as clipped")
	}
	if tally(m.reads, on.Start).Depth != 0 {
		t.Fatal("clipped bases are being counted as evidence")
	}

	run(m, "filters clips")
	off := clipped()
	if off.LeadClip != 0 || off.Start != on.Start+on.LeadClip {
		t.Fatalf("with clips off the read still starts at %d", off.Start)
	}
	if off.Sequence != on.Sequence {
		t.Fatal("the read's own sequence should not change with the setting")
	}
}

func TestTheCommandBarTakesWhatIsTyped(t *testing.T) {
	m := ready(t, "viral_variants")
	press(m, "c") // no key opens a panel any more; it goes into the bar
	if m.panel == "colors" || m.command != "c" {
		t.Fatalf("typing c gave panel %q and command %q", m.panel, m.command)
	}
	if list := m.suggestions(); len(list) != 2 || list[0].label != "contigs" || list[1].label != "colors" {
		t.Fatalf("c suggested %v", list)
	}
	press(m, "down") // the second line, then tab to fill it in
	press(m, "tab")
	if m.command != "colors" {
		t.Fatalf("tab left the bar at %q", m.command)
	}
	press(m, "enter")
	if m.panel != "colors" || !m.panelFocus || m.command != "" {
		t.Fatalf("enter left panel %q, focus %v, bar %q", m.panel, m.panelFocus, m.command)
	}

	m.editCommand(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("gene spi")}) // typed in one burst
	if list := m.suggestions(); len(list) != 1 || list[0].label != "spike" {
		t.Fatalf("gene spi suggested %v", list)
	}
	press(m, "enter")
	spike := m.ref.Features["VIRUS1"]
	for _, f := range spike {
		if geneName(f) == "spike" && f.Kind == "CDS" {
			if m.view.start > f.Start() || m.view.start+m.span() < f.End() {
				t.Fatalf("the spike gene is not on screen: view %v", m.view)
			}
		}
	}
	press(m, "esc")

	run(m, "colornonsense")
	if !strings.Contains(m.message, "colornonsense") {
		t.Fatalf("an unknown command said %q", m.message)
	}
	run(m, "go 4,300")
	if off := abs(m.view.start + m.span()/2 - 4299); off > m.view.bpp {
		t.Fatalf("4,300 landed at %d", m.view.start+m.span()/2)
	}
	for _, junk := range []string{"35", "20M", "M", "80;20"} { // what a busy terminal sends
		m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune(junk)})
		if m.command != "" {
			t.Fatalf("%q opened the bar", junk)
		}
	}
}

func TestTheWheelZooms(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view.bpp = 10
	m.reload()
	wheel := func(button tea.MouseButton, ctrl, shift bool) {
		m.wheelTurned = time.Now().Add(-time.Second) // one click of a mouse wheel, long after the last
		m.handleMouse(tea.MouseMsg{X: 40, Y: m.layout().pileupTop, Button: button, Ctrl: ctrl, Shift: shift})
	}
	wheel(tea.MouseButtonWheelUp, false, false)
	if m.view.bpp != 5 {
		t.Fatalf("one click of the wheel left the zoom at %d bases per column, want 5", m.view.bpp)
	}
	for i := 0; i < 20; i++ { // a track pad sends turns far faster than a hand can click
		m.handleMouse(tea.MouseMsg{X: 40, Y: m.layout().pileupTop, Button: tea.MouseButtonWheelUp})
	}
	if m.view.bpp != 5 {
		t.Fatalf("a track pad's stream of turns ran the zoom on to %d bases per column", m.view.bpp)
	}
	for i := 0; i < 2; i++ { // a fast hand on the wheel is not held back
		m.wheelTurned = time.Now().Add(-2 * mouseClick)
		m.handleMouse(tea.MouseMsg{X: 40, Y: m.layout().pileupTop, Button: tea.MouseButtonWheelDown})
	}
	if m.view.bpp != 20 {
		t.Fatalf("two quick clicks of the wheel left the zoom at %d bases per column, want 20", m.view.bpp)
	}
	wheel(tea.MouseButtonWheelUp, false, false)
	if m.view.bpp != 10 {
		t.Fatalf("the wheel back up left %d bases per column", m.view.bpp)
	}
	for i := 0; i < 3; i++ { // zooming closes in on the pointer, not the middle
		m.wheelTurned = time.Now().Add(-time.Second)
		m.handleMouse(tea.MouseMsg{X: gutterWidth + 10, Y: m.layout().pileupTop, Button: tea.MouseButtonWheelUp})
	}
	if under := m.cols[10].pos; under != m.view.start+10*m.view.bpp {
		t.Fatalf("the base under the pointer moved to column %d", (under-m.view.start)/m.view.bpp)
	}

	start := m.view.start
	wheel(tea.MouseButtonWheelDown, false, true) // shift still pans
	if m.view.start <= start {
		t.Fatal("shift and the wheel did not pan")
	}
	m.rowOffset, m.rowCount = 0, 100
	wheel(tea.MouseButtonWheelDown, true, false) // control still scrolls the reads
	if m.rowOffset != 3 {
		t.Fatalf("control and the wheel left the pileup at row %d", m.rowOffset)
	}
}

func TestAverageQualityFilter(t *testing.T) {
	m := ready(t, "supertest")
	m.view = viewState{0, 20}
	m.reload()
	before := len(m.reads)
	sample := m.reads[0]
	for _, a := range m.reads {
		if a.LeadClip > 0 || a.TailClip > 0 || len(a.Insertions) > 0 || len(a.Quals) != len(a.Bases) {
			continue
		}
		sum := 0
		for _, q := range a.Quals {
			sum += int(q)
		}
		if mean := sum / len(a.Quals); a.MeanQual != mean {
			t.Fatalf("read %s has average quality %d, want %d", a.Name, a.MeanQual, mean)
		}
	}

	run(m, "filters quality 30")
	if m.filter.MinMeanQual != 30 {
		t.Fatalf("the quality filter is %d, want 30", m.filter.MinMeanQual)
	}
	if len(m.reads) >= before {
		t.Fatalf("asking for quality 30 kept %d of %d reads", len(m.reads), before)
	}
	for _, a := range m.reads {
		if a.MeanQual < 30 {
			t.Fatalf("read %s has average quality %d", a.Name, a.MeanQual)
		}
	}
	if summary := strings.Join(m.filterSummary(), " "); !strings.Contains(summary, "Average quality ≥ 30") {
		t.Fatalf("the filter line does not say quality: %s", summary)
	}
	run(m, "filters reset")
	if len(m.reads) != before {
		t.Fatalf("reset gave %d reads, want %d", len(m.reads), before)
	}

	for filterRows[m.filterIndex].name != "Average base quality at least" {
		press(m, "down")
	}
	press(m, "right")
	press(m, "right")
	if m.filter.MinMeanQual != 10 {
		t.Fatalf("two presses right left the quality filter at %d", m.filter.MinMeanQual)
	}
	run(m, "filters quality 90")
	if m.filter.MinMeanQual != 40 {
		t.Fatalf("a quality of 90 was held at %d, want 40", m.filter.MinMeanQual)
	}
	if details := strings.Join(m.readSummary(sample), " "); !strings.Contains(details, "average base quality") {
		t.Fatalf("the read details do not give the average quality: %s", details)
	}
}

func TestOpenAddsTabsAndSwitches(t *testing.T) {
	m := ready(t, "multi_contig_uneven")
	m.view = viewState{0, 1}
	second := filepath.Join(t.TempDir(), "second.sam")
	sam := "@HD\tVN:1.6\n@SQ\tSN:CHR2\tLN:8000\nr1\t0\tCHR2\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
	if err := os.WriteFile(second, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}
	run(m, "open "+second)
	if len(m.tabs) != 2 || m.tab != 1 || m.src != m.tabs[1] {
		t.Fatalf("opening a BAM left %d tabs, on tab %d", len(m.tabs), m.tab)
	}
	if m.src.Path != second {
		t.Fatalf("the second tab is %q", m.src.Path)
	}
	if m.view.start != 0 || m.view.bpp != 1 {
		t.Fatalf("a new tab changed the shared position to %v", m.view)
	}
	top := screen(t, m, 150, 44)[0]
	if !strings.Contains(top, "multi_contig_uneven.bam") || !strings.Contains(top, "second.sam") {
		t.Fatalf("the top bar does not show both files: %q", top)
	}
	press(m, "t")
	if m.tab != 0 {
		t.Fatalf("t left the view on %d", m.tab)
	}
	press(m, "t")
	if m.tab != 1 {
		t.Fatalf("t again left the view on %d", m.tab)
	}

	// The tabs along the top take a click.
	mouse(m, m.tabHits[0].from+1, 0, tea.MouseActionPress, tea.MouseButtonLeft)
	if m.tab != 0 {
		t.Fatalf("clicking the first tab left the view on %d", m.tab)
	}
	mouse(m, m.tabHits[1].from+1, 0, tea.MouseActionPress, tea.MouseButtonLeft)
	if m.tab != 1 {
		t.Fatalf("clicking the second tab left the view on %d", m.tab)
	}

	src, err := OpenAlignments(second)
	if err != nil {
		t.Fatal(err)
	}
	contig := src.Contigs[0]
	src.Close()
	m.Update(profileMsg{m.tabs[1], map[string][]float64{contig.Name: make([]float64, profileBins)}})
	if m.tabs[1].Contigs[0].Profile == nil {
		t.Fatal("a coverage sketch did not reach the tab it was measured for")
	}
}

func TestOpenShowsReferenceAndVariants(t *testing.T) {
	m := ready(t, "viral_variants")
	m.ref, m.refPath, m.variants = nil, "", nil
	m.highlight = namedIndex(highlights, "Mismatch (from consensus)")
	m.schemeRow = m.highlight
	m.panel = ""
	m.reload()

	// A reference opened where there was none brings its features with it.
	run(m, "open testdata/viral_variants.gb")
	if m.ref == nil || !strings.HasSuffix(m.refPath, "viral_variants.gb") {
		t.Fatalf("opening a reference left %q", m.refPath)
	}
	m.view = viewState{3199, 1}
	m.reload()
	if m.layout().geneRows == 0 {
		t.Fatal("the opened reference has no genes on screen")
	}
	if !strings.Contains(screen(t, m, 150, 44)[0], "reference viral_variants.gb") {
		t.Fatal("the top bar does not name the reference")
	}

	// A VCF opened where there was none brings its variants with it.
	run(m, "open testdata/viral_variants.vcf")
	if len(m.variants) != 9 || m.panel != "variants" {
		t.Fatalf("opening a VCF gave %d variants, panel %q", len(m.variants), m.panel)
	}
	if m.variants[0].Impact.Gene == "" {
		t.Fatal("the opened variants were not worked out against the reference")
	}
}

func TestOpenOffersExamplesAndIsInTheHelp(t *testing.T) {
	m := ready(t, "viral_variants")
	for _, r := range "open " {
		press(m, string(r))
	}
	list := m.suggestions()
	if len(list) != 3 {
		t.Fatalf("open offered %v", list)
	}
	for i, want := range []string{"reads.bam", "genes.gbk", "calls.vcf"} {
		if list[i].label != want {
			t.Fatalf("open example %d is %q, want %q", i, list[i].label, want)
		}
	}
	if list[0].help != "another BAM file in a new tab" {
		t.Fatalf("the BAM example says %q", list[0].help)
	}
	found := false
	for _, line := range helpLines {
		if line[0] == "open [file]" {
			found = true
		}
	}
	if !found {
		t.Fatal("the help does not mention open")
	}
}

func TestOpenRefusesFilesThatDoNotMatch(t *testing.T) {
	m := ready(t, "supertest")
	ref := m.ref
	run(m, "open testdata/bf_fragilis.gb")
	if m.ref != ref || !strings.Contains(m.message, "does not match") {
		t.Fatalf("a foreign reference was taken: %q", m.message)
	}
	run(m, "open testdata/bf_fragilis.bam")
	if len(m.tabs) != 1 || !strings.Contains(m.message, "Cannot open bf_fragilis.bam: does not share reference contigs with supertest.bam") {
		t.Fatalf("a foreign BAM was taken: %d tabs, %q", len(m.tabs), m.message)
	}
	before := len(m.variants)
	run(m, "open testdata/bf_fragilis.vcf")
	if len(m.variants) != before || !strings.Contains(m.message, "does not match") {
		t.Fatalf("a foreign VCF was taken: %q", m.message)
	}
	run(m, "open ")
	if !strings.Contains(m.message, "Type a file") {
		t.Fatalf("open with no file said %q", m.message)
	}
	run(m, "open testdata/missing.bam")
	if !strings.Contains(m.message, "Cannot read") {
		t.Fatalf("a missing file said %q", m.message)
	}
	run(m, "open notes.txt")
	if m.message != "Error: wrong filetype. Please open a BAM, Genbank, or VCF." {
		t.Fatalf("an unknown filetype said %q", m.message)
	}

	// A file with a contig the first one does not have cannot join it.
	other := ready(t, "multi_contig_uneven")
	run(other, "open testdata/supertest.bam")
	if len(other.tabs) != 1 || !strings.Contains(other.message,
		"Cannot open supertest.bam: does not share reference contigs with multi_contig_uneven.bam") {
		t.Fatalf("a BAM with an extra contig was taken: %d tabs, %q", len(other.tabs), other.message)
	}
}

func TestCommandLineTakesSeveralAlignments(t *testing.T) {
	in, err := readInputs([]string{"a.bam", "b.sam", "genes.gbk", "calls.vcf"})
	if err != nil {
		t.Fatal(err)
	}
	if len(in.alignments) != 2 || in.alignments[0] != "a.bam" || in.alignments[1] != "b.sam" {
		t.Fatalf("the alignments came out %v", in.alignments)
	}
	if in.reference != "genes.gbk" || in.variants != "calls.vcf" {
		t.Fatalf("the reference and VCF came out %q and %q", in.reference, in.variants)
	}
	for _, args := range [][]string{
		{"a.gbk", "b.fa"},  // two references
		{"a.vcf", "b.vcf"}, // two VCFs
		{"reads.bam", "notes.txt"},
	} {
		if _, err := readInputs(args); err == nil {
			t.Fatalf("%v was accepted", args)
		}
	}
}

func TestSeveralBamsShareOneReference(t *testing.T) {
	src, err := OpenAlignments("testdata/shared_a.bam")
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(src.Close)
	ref, err := LoadReference("testdata/shared_ref.gb")
	if err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs(src.Contigs)
	m := newModel(src, ref, "testdata/shared_ref.gb", nil)
	m.width, m.height = 150, 44
	m.view = viewState{1200, 1}
	m.reload()
	start := time.Now()
	for _, path := range []string{"testdata/shared_b.bam", "testdata/shared_c.bam"} {
		m.openAlignments(path)
	}
	t.Logf("two more BAMs opened in %s", time.Since(start).Round(time.Microsecond))
	if len(m.tabs) != 3 {
		t.Fatalf("three BAMs made %d tabs: %q", len(m.tabs), m.message)
	}
	for _, tab := range m.tabs {
		if len(tab.Contigs) != 3 {
			t.Fatalf("%s has %d contigs, want 3", tab.Path, len(tab.Contigs))
		}
	}
	top := screen(t, m, 150, 44)[0]
	for _, name := range []string{"shared_a.bam", "shared_b.bam", "shared_c.bam", "reference shared_ref.gb"} {
		if !strings.Contains(top, name) {
			t.Fatalf("the top bar does not show %q: %q", name, top)
		}
	}

	m.switchTab(0)
	first := len(m.reads)
	m.switchTab(1)
	second := len(m.reads)
	m.switchTab(2)
	third := len(m.reads)
	if first == 0 || second == 0 || third == 0 || first == third {
		t.Fatalf("the three tabs hold %d, %d and %d reads", first, second, third)
	}
	if m.view.start != 1200 {
		t.Fatalf("switching tabs moved the shared position to %d", m.view.start)
	}
	if m.layout().geneRows == 0 {
		t.Fatal("the shared reference has no genes on screen")
	}
}

func TestReadLengthFiltersAndColors(t *testing.T) {
	m := ready(t, "supertest")
	m.view = viewState{0, 20}
	m.reload()
	before, longest := len(m.reads), 0
	for _, a := range m.reads {
		longest = max(longest, len(a.Sequence))
	}
	run(m, "filters length 200")
	if m.filter.MinLength != 200 {
		t.Fatalf("the length filter is %d, want 200", m.filter.MinLength)
	}
	if len(m.reads) >= before {
		t.Fatalf("asking for 200 bases kept %d of %d reads", len(m.reads), before)
	}
	for _, a := range m.reads {
		if len(a.Sequence) < 200 {
			t.Fatalf("read %s is only %d bases", a.Name, len(a.Sequence))
		}
	}
	run(m, "filters reset")
	if len(m.reads) != before {
		t.Fatal("reset did not bring the short reads back")
	}

	press(m, "down") // the arrows walk the same filter through its ladder
	for m.filterIndex < len(filterRows)-1 && filterRows[m.filterIndex].name != "Read length at least" {
		press(m, "down")
	}
	press(m, "right")
	press(m, "right")
	if m.filter.MinLength != 50 {
		t.Fatalf("two presses right left the length filter at %d, want 50", m.filter.MinLength)
	}

	scheme := 0
	for i, s := range colors {
		if s == "Read length" {
			scheme = i
		}
	}
	short, long := m.reads[0], m.reads[0]
	for _, a := range m.reads {
		if len(a.Sequence) < len(short.Sequence) {
			short = a
		}
		if len(a.Sequence) > len(long.Sequence) {
			long = a
		}
	}
	if readColor(scheme, short, short.Start, m.readLength) == readColor(scheme, long, long.Start, m.readLength) {
		t.Fatalf("a %d base read and a %d base one came out the same color", len(short.Sequence), len(long.Sequence))
	}
}

func TestKeysOutsideTheCommandBarLeaveIt(t *testing.T) {
	m := ready(t, "viral_variants")
	press(m, "g")
	press(m, "e")
	start := m.view.start
	press(m, "right") // an arrow key belongs to the view, not the bar
	if m.command != "" {
		t.Errorf("the bar still holds %q", m.command)
	}
	if m.view.start == start {
		t.Error("the arrow key did not pan")
	}
	press(m, "g")
	press(m, "?") // so does the help key
	if m.command != "" || m.panel != "help" {
		t.Errorf("bar %q, panel %q", m.command, m.panel)
	}
	press(m, "g")
	m.handleMouse(tea.MouseMsg{Action: tea.MouseActionPress, Button: tea.MouseButtonLeft, X: 10, Y: 20})
	if m.command != "" {
		t.Errorf("clicking left the bar holding %q", m.command)
	}
	press(m, "g")
	m.handleMouse(tea.MouseMsg{Action: tea.MouseActionPress, Button: tea.MouseButtonWheelDown, X: 10, Y: 20})
	if m.command != "" {
		t.Errorf("the wheel left the bar holding %q", m.command)
	}
	m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("\x1b[<35;80;20M")}) // a sequence the terminal sent
	if m.command != "" {
		t.Errorf("the bar took %q", m.command)
	}
}

// The viral reads rd00001 to rd00017 are unpaired reads that happen to share
// a name with another read far away. A shared name is only a pair when the
// reads say they are paired, so these must be packed and lit up on their own.
func TestUnpairedReadsWithTheSameNameStayApart(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{0, 50}
	m.reload()
	var rowEnd []int // the packer's own bookkeeping, replayed to check every read took the first free row
	for _, a := range m.reads {
		for row := 0; row < a.Row; row++ {
			if row < len(rowEnd) && a.Start >= rowEnd[row] {
				t.Fatalf("%s at %d sits on row %d though row %d was free", a.Name, a.Start, a.Row, row)
			}
		}
		for len(rowEnd) <= a.Row {
			rowEnd = append(rowEnd, 0)
		}
		rowEnd[a.Row] = a.End + m.view.bpp
	}
	read := m.reads[0]
	l := m.layout()
	mouse(m, gutterWidth+max(m.column(read.Start), 0), l.pileupTop+read.Row-m.rowOffset, tea.MouseActionMotion, tea.MouseButtonNone)
	m.View()
	lit := 0
	for _, a := range m.reads {
		if a.Row == read.Row && m.frame[l.pileupTop+read.Row][gutterWidth+max(m.column(a.Start), 0)].s.bg != colorBackdrop {
			lit++
		}
	}
	if lit != 1 {
		t.Fatalf("hovering %s lit %d reads on its row, want just itself", read.Name, lit)
	}
}

// A region with more reads than fit is loaded slice by slice. A read spanning
// two slices must be loaded once, not once for each slice it touches.
func TestSampledRegionsLoadEachReadOnce(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{0, 20}
	m.reload()
	whole := len(m.reads)
	m.src.Contigs[0].Reads = 1 << 30 // pretend the contig is enormous so the region is sampled
	m.reload()
	seen := map[string]bool{}
	for _, a := range m.reads {
		if at := a.Name + ":" + comma(a.Start); seen[at] {
			t.Fatalf("%s was loaded twice", at)
		} else {
			seen[at] = true
		}
	}
	if !m.truncated || len(m.reads) != whole {
		t.Fatalf("sampling gave %d reads, sampled %v, want all %d", len(m.reads), m.truncated, whole)
	}
}

// Choosing a command that takes more fills it in and lists what can follow,
// instead of running it with nothing; the list is a menu two levels deep.
func TestChoosingACommandListsWhatItTakes(t *testing.T) {
	m := ready(t, "viral_variants")
	press(m, "f")
	press(m, "enter") // filters says something on its own, so it runs and opens its panel
	if m.command != "" || m.panel != "filters" {
		t.Fatalf("enter on filters left the bar at %q and the panel %q", m.command, m.panel)
	}
	for _, r := range "gene" {
		press(m, string(r))
	}
	press(m, "enter") // gene wants a name, so it is filled in and its list opens
	if m.command != "gene " || len(m.suggestions()) == 0 {
		t.Fatalf("enter on gene left the bar at %q with %d names", m.command, len(m.suggestions()))
	}
	press(m, "esc")
	for _, r := range "filters " {
		press(m, string(r))
	}
	if list := m.suggestions(); len(list) != 9 || list[0].label != "mapq 30" {
		t.Fatalf("filters listed %v", list)
	}
	press(m, "enter") // the first filter, which wants a number
	if m.command != "filters mapq " {
		t.Fatalf("enter on mapq left the bar at %q", m.command)
	}
	run(m, "30")
	if m.filter.MinMapQ != 30 || m.panel != "filters" || m.command != "" {
		t.Fatalf("mapq 30 gave %d, panel %q, bar %q", m.filter.MinMapQ, m.panel, m.command)
	}
	for _, r := range "filters mapq 40" {
		press(m, string(r))
	}
	press(m, "tab") // tab must not throw away the number
	if m.command != "filters mapq 40" {
		t.Fatalf("tab left the bar at %q", m.command)
	}
	press(m, "enter")
	if m.filter.MinMapQ != 40 {
		t.Fatalf("mapq 40 gave %d", m.filter.MinMapQ)
	}
	press(m, "g")
	press(m, "e")
	press(m, "enter")
	if list := m.suggestions(); m.command != "gene " || len(list) == 0 {
		t.Fatalf("enter on gene left %q listing %v", m.command, list)
	}
	press(m, "esc")
	press(m, "g")
	press(m, "enter")
	if list := m.suggestions(); m.command != "go " || len(list) != 3 {
		t.Fatalf("enter on go left %q listing %v", m.command, list)
	}
	press(m, "enter") // with nothing typed there is nowhere to go, and the bar says so
	if m.command != "" || !strings.Contains(m.message, "position") {
		t.Fatalf("enter on an empty go left %q saying %q", m.command, m.message)
	}
}

// TestUsageListsOnlyFileArguments pins the command line text: files of each
// kind, in the order they were written, and no options at all.
func TestUsageListsOnlyFileArguments(t *testing.T) {
	for _, line := range []string{
		"dnacli: a modern terminal read alignment viewer.",
		"version " + version,
		"  dnacli reads.bam [reference.gb|reference.fa]",
		"  dnacli reads.bam [reference.gb|reference.fa] [variants.vcf]",
		"  dnacli reads.bam [more.bam ...]",
		"  dnacli sample_a.bam sample_b.bam genome.gbk snvs.vcf",
		"  dnacli samples*.bam",
		"Files are recognised by their extensions, and can be in any order.",
	} {
		if !strings.Contains(usage, line) {
			t.Fatalf("the usage does not say %q", line)
		}
	}
	if strings.Contains(usage, "Options:") || strings.Contains(usage, "--") {
		t.Fatalf("the usage lists more than the files:\n%s", usage)
	}
}
