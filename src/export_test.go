package main

import (
	"bytes"
	"encoding/xml"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// svgText is one <text> element of an exported picture.
type svgText struct {
	X    float64 `xml:"x,attr"`
	Y    float64 `xml:"y,attr"`
	Fill string  `xml:"fill,attr"`
	Body string  `xml:",chardata"`
}

type svgRect struct {
	X    float64 `xml:"x,attr"`
	Y    float64 `xml:"y,attr"`
	W    float64 `xml:"width,attr"`
	H    float64 `xml:"height,attr"`
	Fill string  `xml:"fill,attr"`
}

type svgDoc struct {
	Texts []svgText `xml:"g>text"`
	Rects []svgRect `xml:"g>rect"`
}

// TestExportDrawsTheSameFrameAsTheScreen checks that every letter of the
// exported picture sits in the cell it filled on screen, in the same color,
// and that the block, line, arrow and marker characters became shapes.
func TestExportDrawsTheSameFrameAsTheScreen(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{250, 1}
	screen(t, m, 120, 30)
	path := filepath.Join(t.TempDir(), "view.svg")
	m.runExport(path)
	if !strings.Contains(m.message, "view.svg") || !strings.Contains(m.message, "light background") {
		t.Fatalf("the export said %q; light is the default", m.message)
	}
	if err := m.writeSVG(path, darkSVG); err != nil { // the dark paper keeps the screen's colors
		t.Fatal(err)
	}

	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if bytes.Contains(data, []byte{0x1b}) {
		t.Fatal("the picture holds a terminal escape")
	}
	doc := svgDoc{}
	if err := xml.Unmarshal(data, &doc); err != nil {
		t.Fatalf("the picture is not well formed: %v", err)
	}
	if len(doc.Texts) == 0 {
		t.Fatal("the picture has no text")
	}
	for _, r := range "▁▂▃▄▅▆▇█─━│·▪▼▶◀▸" {
		if strings.ContainsRune(string(data), r) {
			t.Fatalf("%q was left as a letter instead of a shape", r)
		}
	}
	if !strings.Contains(string(data), `id="background"`) || !strings.Contains(string(data), `<title>dnacli VIRUS1:251-`) {
		t.Fatal("the picture has no background group or title")
	}
	l := m.layout()
	top, bottom, right := m.contentRect(l)
	for _, r := range doc.Rects {
		if r.H > svgCellH+0.1 || r.Y < svgPad-0.1 || r.Y+r.H > svgPad+float64(bottom-top)*svgCellH+0.1 {
			t.Fatalf("a shape spills out of its cell: %+v", r)
		}
	}
	// The bars around the tracks, the tabs and any panel are left out.
	for _, gone := range []string{"viral_variants.bam", "type a command", "Visualizing", "←→ pan", "Contigs"} {
		if strings.Contains(string(data), gone) {
			t.Fatalf("the picture still holds %q", gone)
		}
	}

	// Put every letter back where it was drawn and compare with the screen.
	type cellText struct{ fill, body string }
	drawn := map[[2]int]cellText{}
	for _, run := range doc.Texts {
		row := int(math.Round((run.Y-svgPad-svgBaseline)/svgCellH)) + top
		col := int(math.Round((run.X - svgPad) / svgCellW))
		for i, r := range []rune(run.Body) {
			if row < top || row >= bottom || col+i >= right {
				t.Fatalf("%q was drawn outside the tracks at %d,%d", r, row, col+i)
			}
			drawn[[2]int{row, col + i}] = cellText{run.Fill, string(r)}
		}
	}
	for row := top; row < bottom; row++ {
		for col := 0; col < right; col++ {
			cell := m.frame[row][col]
			if !drawsText(cell.r) {
				continue
			}
			got := drawn[[2]int{row, col}]
			if got.body != string(cell.r) || got.fill != cell.s.fg {
				t.Fatalf("cell %d,%d drew %q in %q, want %q in %q",
					row, col, got.body, got.fill, cell.r, cell.s.fg)
			}
		}
	}
}

// TestExportLightThemePrintsOnWhite checks that the light theme has a white
// page and dark letters, with the same layout.
func TestExportLightThemePrintsOnWhite(t *testing.T) {
	m := ready(t, "supertest")
	m.view = viewState{0, 1}
	screen(t, m, 120, 30)
	path := filepath.Join(t.TempDir(), "print.svg")
	m.runExport("light " + path)
	if !strings.Contains(m.message, "light background") {
		t.Fatalf("the export said %q", m.message)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(data), `fill="#ffffff"`) {
		t.Fatal("the light theme has no white page")
	}
	if strings.Contains(string(data), chosenPalette.bases['A']) {
		t.Fatalf("the light theme kept the bright base color %s", chosenPalette.bases['A'])
	}
	doc := svgDoc{}
	if err := xml.Unmarshal(data, &doc); err != nil {
		t.Fatal(err)
	}
	for _, run := range doc.Texts {
		if c := contrastWhite(run.Fill); c < 4.4 {
			t.Fatalf("%q is too pale to read on white: %s is %.1f:1", run.Body, run.Fill, c)
		}
	}
	// The palette keeps its character: the green of A stays green, only darker.
	hue, saturation, _ := hsl(chosenPalette.bases['A'])
	shade := lightSVG.ink(chosenPalette.bases['A'])
	shadeHue, shadeSaturation, _ := hsl(shade)
	if math.Abs(hue-shadeHue) > 0.03 || shadeSaturation < 0.5*saturation {
		t.Fatalf("the base color %s turned into %s", chosenPalette.bases['A'], shade)
	}
}

// TestExportLeavesOutFilteredReads checks that the picture follows the filters
// in force, like the screen does.
func TestExportLeavesOutFilteredReads(t *testing.T) {
	m := ready(t, "supertest")
	m.panel = "" // the contigs panel would cover some of the read rows
	m.view = viewState{0, 1}
	screen(t, m, 120, 40)
	dir := t.TempDir()
	m.runExport(filepath.Join(dir, "all.svg"))
	all := pileupText(t, m, filepath.Join(dir, "all.svg"))

	m.filter.MinMeanQual = 100 // higher than any base quality in the file
	m.reload()
	screen(t, m, 120, 40)
	m.runExport(filepath.Join(dir, "none.svg"))
	none := pileupText(t, m, filepath.Join(dir, "none.svg"))
	if none >= all/2 || none > 40 { // only the "no reads" line is left
		t.Fatalf("filtered reads still fill the picture: %d of %d letters", none, all)
	}
}

// pileupText counts the letters an exported picture draws over the reads.
func pileupText(t *testing.T, m *model, path string) int {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	doc := svgDoc{}
	if err := xml.Unmarshal(data, &doc); err != nil {
		t.Fatal(err)
	}
	l := m.layout()
	top, _, _ := m.contentRect(l)
	count := 0
	for _, run := range doc.Texts {
		row := math.Round((run.Y-svgPad-svgBaseline)/svgCellH) + float64(top)
		if row >= float64(l.pileupTop) && row < float64(l.pileupTop+l.pileupRows) {
			count += len([]rune(run.Body))
		}
	}
	return count
}

// TestExportFileNameNamesTheView checks the file an export writes when no name
// is given.
func TestExportFileNameNamesTheView(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{299, 1}
	m.reload()
	want := "dnacli-VIRUS1-300-" + comma(m.view.start+m.span()) + ".svg"
	want = strings.ReplaceAll(want, ",", "")
	if got := m.exportName(); got != want {
		t.Fatalf("the file would be called %q, want %q", got, want)
	}
}

// TestExportIsOnline checks that the command is in the help and offers examples.
func TestExportIsOnline(t *testing.T) {
	m := ready(t, "viral_variants")
	for _, r := range "export " {
		press(m, string(r))
	}
	list := m.suggestions()
	if len(list) != 2 || list[0].label != "figure.svg light" || list[1].label != "figure.svg dark" {
		t.Fatalf("export offered %v", list)
	}
	if list[0].help != "save view as SVG (light)" || list[1].help != "save view as SVG (dark)" {
		t.Fatalf("export described its examples as %q and %q", list[0].help, list[1].help)
	}
	for _, line := range helpLines {
		if strings.HasPrefix(line[0], "export") {
			return
		}
	}
	t.Fatal("the help does not mention export")
}

// TestExportPerformance times a picture twice the usual window, at the largest
// size the interface allows.
func TestExportPerformance(t *testing.T) {
	m := ready(t, "supertest")
	screen(t, m, 200, 60)
	path := filepath.Join(t.TempDir(), "big.svg")
	start := time.Now()
	m.writeSVG(path, darkSVG)
	elapsed := time.Since(start)
	info, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("200 × 60 cells: %v, %d kB", elapsed.Round(time.Millisecond), info.Size()/1000)
	if elapsed > 250*time.Millisecond {
		t.Fatalf("exporting took %v", elapsed)
	}
}

// TestExportFileNameKeepsItsDots checks that a file name can be typed into the
// bar: the dot and slash of a path must not throw the command away.
func TestExportFileNameKeepsItsDots(t *testing.T) {
	m := ready(t, "viral_variants")
	dir := t.TempDir()
	t.Chdir(dir)
	for _, r := range "export example.svg" {
		press(m, string(r))
	}
	if m.command != "export example.svg" {
		t.Fatalf("the bar holds %q", m.command)
	}
	press(m, "enter")
	if _, err := os.Stat(filepath.Join(dir, "example.svg")); err != nil {
		t.Fatalf("the picture was not written: %v", err)
	}

	for _, path := range []string{"a/b_c-d.svg", "light ../x~1.svg"} {
		m.command = ""
		for _, r := range "export " + path {
			press(m, string(r))
		}
		if m.command != "export "+path {
			t.Fatalf("the bar holds %q, want %q", m.command, "export "+path)
		}
		m.command = ""
	}
}
