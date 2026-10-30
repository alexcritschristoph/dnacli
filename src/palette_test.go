package main

import (
	"strings"
	"testing"
	"time"

	tea "github.com/charmbracelet/bubbletea"
)

// TestColorPaletteRepaintsEveryScheme checks the palette section of the
// coloring panel: choosing a palette leaves the highlighting and the color
// scheme alone, and every scheme then draws with the palette's colors.
func TestColorPaletteRepaintsEveryScheme(t *testing.T) {
	m := ready(t, "supertest")
	run(m, "colors")
	first := len(highlights) + len(colors)

	// The section lists every palette, in order, below the other two.
	var names []string
	for i, line := range coloringPanel {
		if line.title == "Color palette" && i != len(coloringPanel)-len(palettes)-1 {
			t.Fatalf("the palette title is on line %d", i)
		}
		if line.entry >= first {
			names = append(names, line.name)
		}
	}
	if strings.Join(names, ",") != "Default,Color-blind safe,High contrast,Sunset,Meadow" {
		t.Fatalf("the palette section lists %v", names)
	}

	// Moving onto another palette previews it, without touching anything else.
	beforeColor, beforeHighlight := m.color, m.highlight
	m.focusColoring(first + 1)
	if m.palette != 1 || chosenPalette.name != "Color-blind safe" {
		t.Fatalf("previewing gave palette %d (%s)", m.palette, chosenPalette.name)
	}
	if m.color != beforeColor || m.highlight != beforeHighlight {
		t.Fatal("previewing a palette changed the color scheme or the highlighting")
	}
	m.pressColoring() // and enter locks it
	if !m.coloringLocked || chosenPalette.name != "Color-blind safe" {
		t.Fatal("enter did not lock the previewed palette")
	}
	m.pressColoring() // a palette is always in force, so pressing it again keeps it
	if m.palette != 1 {
		t.Fatal("pressing the palette again cleared it")
	}
	m.focusColoring(first + 2) // locked, so the cursor alone moves
	if m.palette != 1 || chosenPalette.name != "Color-blind safe" {
		t.Fatal("moving the cursor changed a locked palette")
	}

	// Every scheme draws from the palette now.
	contig := m.contigs()[m.contig].Name
	a := m.reads[0]
	for pos := a.Start; pos < a.End; pos++ {
		base := a.BaseAt(pos)
		if base != 0 && readColor(namedIndex(colors, "Nucleotides"), a, pos, m.readLength) != chosenPalette.bases[base] {
			t.Fatalf("%c was drawn in the default color", base)
		}
	}
	forward, reverse := chosenPalette.forward, chosenPalette.reverse
	if readColor(namedIndex(colors, "Read strand"), a, a.Start, m.readLength) != forward {
		t.Fatal("a forward read did not take the palette's forward strand color")
	}
	b := m.reads[0]
	for _, r := range m.reads {
		if r.Reverse() {
			b = r
			break
		}
	}
	if readColor(namedIndex(colors, "Read strand"), b, b.Start, m.readLength) != reverse {
		t.Fatal("a reverse read did not take the palette's reverse strand color")
	}
	if got := impactColor(Impact{Kind: "synonymous"}); got != chosenPalette.good {
		t.Fatalf("a synonymous change is drawn %s, want the palette's %s", got, chosenPalette.good)
	}
	if got := impactColor(Impact{Kind: "nonsense"}); got != chosenPalette.alert {
		t.Fatalf("a nonsense change is drawn %s, want the palette's %s", got, chosenPalette.alert)
	}

	// The reference row and a mismatching read base are drawn in the palette
	// too, and the panel marks the palette as chosen.
	m.view = viewState{0, 1}
	m.reload()
	m.View()
	l := m.layout()
	for col := 0; col < l.columns; col++ {
		pos := m.view.start + col
		base := m.ref.Base(contig, pos)
		if base == 0 || chosenPalette.bases[base] == "" {
			continue
		}
		if got := m.frame[l.refRow][gutterWidth+col].s.fg; got != chosenPalette.bases[base] {
			t.Fatalf("the reference base at %d drew %s, want %s", pos, got, chosenPalette.bases[base])
		}
		break
	}
	for i, line := range coloringPanel {
		if line.name != "Color-blind safe" {
			continue
		}
		if bg := m.frame[l.listTop+i][m.width-3].s.bg; bg != colorLine {
			t.Fatalf("the chosen palette is drawn with background %s", bg)
		}
	}

	// The strip says which palette is on, and pressing another switches to it.
	if summary := m.coloringSummary(); !strings.Contains(summary, "Color-blind safe palette") {
		t.Fatalf("the coloring summary is %q", summary)
	}
	m.pickColoring(first + 4)
	if m.palette != 4 || chosenPalette.name != "Meadow" {
		t.Fatalf("clicking a palette gave %d (%s)", m.palette, chosenPalette.name)
	}
	click := func(name string) {
		for i, line := range coloringPanel {
			if line.name == name {
				mouse(m, m.width-4, l.listTop+i, tea.MouseActionPress, tea.MouseButtonLeft)
				return
			}
		}
		t.Fatalf("the panel has no row called %s", name)
	}
	click("Default")
	if m.palette != 0 || chosenPalette.name != "Default" {
		t.Fatalf("clicking Default gave %d (%s)", m.palette, chosenPalette.name)
	}
}

// TestColorPalettePerformance draws a frame in each palette, so that a palette
// cannot quietly make drawing slower.
func TestColorPalettePerformance(t *testing.T) {
	m := ready(t, "supertest")
	m.view = viewState{0, 1}
	for i, p := range palettes {
		m.choosePalette(i)
		start := time.Now()
		m.reload()
		m.View()
		elapsed := time.Since(start)
		t.Logf("%-16s %6d reads  %s", p.name, len(m.reads), elapsed.Round(time.Microsecond))
		if elapsed > 100*time.Millisecond {
			t.Fatalf("drawing a frame in the %s palette took %v", p.name, elapsed)
		}
	}
}

// TestColorPalettePanelScrollsInAShortWindow checks that the last palettes can
// still be reached when the window is too short to show the whole list.
func TestColorPalettePanelScrollsInAShortWindow(t *testing.T) {
	m := ready(t, "supertest")
	m.width, m.height = 100, 30
	run(m, "colors")
	m.focusColoring(len(highlights) + len(colors) + len(palettes) - 1) // Meadow, the last entry
	if m.palette != len(palettes)-1 {
		t.Fatalf("the last palette gave %d", m.palette)
	}
	m.View()
	l := m.layout()
	offset := m.coloringOffset(l)
	if offset == 0 {
		t.Fatal("the panel did not scroll to the cursor")
	}
	shown := coloringPanel[offset : offset+l.listRows]
	for i, line := range shown {
		if line.name == "Meadow" && !strings.Contains(string(frameRow(m, l.listTop+i)), "Meadow") {
			t.Fatal("the last palette is not drawn on screen")
		}
	}
}
