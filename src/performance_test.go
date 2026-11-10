package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	tea "github.com/charmbracelet/bubbletea"
)

// TestPerformance records how long the work behind one screen takes, so that a
// change that slows the program down shows up in the test output.
func TestPerformance(t *testing.T) {
	for _, stem := range []string{"viral_variants", "supertest", "multi_contig_uneven"} {
		m := testModel(t, stem)
		m.width, m.height = 150, 44
		for _, bpp := range []int{1, 20, 500} {
			m.view = viewState{0, bpp}
			load := time.Now()
			m.reload()
			loaded := time.Since(load)
			frame := time.Now()
			m.View()
			drawn := time.Since(frame)
			t.Log(fmt.Sprintf("%-20s %5d bases per column  %6d reads  load %7s  draw %7s",
				stem, bpp, len(m.reads), loaded.Round(time.Microsecond), drawn.Round(time.Microsecond)))
			if loaded+drawn > 300*time.Millisecond {
				t.Errorf("%s at %d bases per column took %v", stem, bpp, loaded+drawn)
			}
		}
	}
}

// TestQualityFilterPerformance reads the same regions with the average quality
// filter off and on, so that the cost of the read's average shows beside the
// cost of the reads themselves.
func TestQualityFilterPerformance(t *testing.T) {
	for _, stem := range []string{"viral_variants", "supertest", "multi_contig_uneven"} {
		for _, minimum := range []int{0, 30} {
			m := testModel(t, stem)
			m.width, m.height = 150, 44
			m.filter.MinMeanQual = minimum
			for _, bpp := range []int{1, 20, 500} {
				m.view = viewState{0, bpp}
				load := time.Now()
				m.reload()
				loaded := time.Since(load)
				frame := time.Now()
				m.View()
				drawn := time.Since(frame)
				t.Log(fmt.Sprintf("%-20s quality at least %2d  %5d bases per column  %6d reads  load %7s  draw %7s",
					stem, minimum, bpp, len(m.reads), loaded.Round(time.Microsecond), drawn.Round(time.Microsecond)))
				if loaded+drawn > 300*time.Millisecond {
					t.Errorf("%s at %d bases per column with quality at least %d took %v", stem, bpp, minimum, loaded+drawn)
				}
			}
		}
	}
}

// TestOpenFilesPerformance records the cost of opening each kind of file: an
// alignment file joins the tabs, a reference or a VCF is read and matched.
func TestOpenFilesPerformance(t *testing.T) {
	second := filepath.Join(t.TempDir(), "second.sam")
	sam := "@HD\tVN:1.6\n@SQ\tSN:CHR2\tLN:8000\nr1\t0\tCHR2\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
	if err := os.WriteFile(second, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}
	cases := []struct{ stem, path string }{
		{"multi_contig_uneven", second},
		{"viral_variants", "testdata/viral_variants.gb"},
		{"viral_variants", "testdata/viral_variants.vcf"},
	}
	for _, c := range cases {
		m := testModel(t, c.stem)
		m.width, m.height = 150, 44
		m.reload()
		start := time.Now()
		m.runOpen(c.path)
		took := time.Since(start)
		t.Log(fmt.Sprintf("%-32s opened in %7s, tabs %d", c.path, took.Round(time.Microsecond), len(m.tabs)))
		if took > 300*time.Millisecond {
			t.Errorf("opening %s took %v", c.path, took)
		}
	}
}

// TestMismatchModesPerformance draws the same regions with both mismatch
// colorings, so that the cost of finding the consensus shows up beside the cost
// of reading the reference.
func TestMismatchModesPerformance(t *testing.T) {
	for _, stem := range []string{"viral_variants", "supertest"} {
		for _, name := range []string{"Mismatch (from consensus)", "Mismatch (from reference)"} {
			m := testModel(t, stem)
			m.width, m.height = 150, 44
			m.highlight = namedIndex(highlights, name)
			for _, bpp := range []int{1, 20, 500} {
				m.view = viewState{0, bpp}
				load := time.Now()
				m.reload()
				loaded := time.Since(load)
				frame := time.Now()
				m.View()
				drawn := time.Since(frame)
				t.Log(fmt.Sprintf("%-14s %-26s %5d bases per column  %6d reads  load %7s  draw %7s",
					stem, name, bpp, len(m.reads), loaded.Round(time.Microsecond), drawn.Round(time.Microsecond)))
				if loaded+drawn > 300*time.Millisecond {
					t.Errorf("%s with %s at %d bases per column took %v", stem, name, bpp, loaded+drawn)
				}
			}
		}
	}
}

// TestDenseContigCheckpoints checks that a region deep inside a small, dense
// contig starts at a checkpoint instead of at the contig's first record once
// the coverage sketch has run, and that it still finds exactly the same reads.
func TestDenseContigCheckpoints(t *testing.T) {
	src, err := OpenAlignments("testdata/dense_mito.bam")
	if err != nil {
		t.Skip(err)
	}
	defer src.Close()
	if src.index == nil {
		t.Skip("dense_mito.bam has no index")
	}
	var contig Contig
	for _, c := range src.Contigs {
		if c.Reads > 0 {
			contig = c
			break
		}
	}
	if contig.Length == 0 {
		t.Skip("dense_mito.bam has no mapped reads")
	}
	start := contig.Length / 2
	want, wantSeen, _ := src.Reads(contig.Name, start, start+100, defaultFilter(), nil, true)

	begin := time.Now()
	m := newModel(src, nil, "", nil)
	m.width, m.height = 150, 44
	m.Update(sketchCoverage(src)())
	sketch := time.Since(begin)

	begin = time.Now()
	got, gotSeen, _ := src.Reads(contig.Name, start, start+100, defaultFilter(), nil, true)
	read := time.Since(begin)

	if len(src.checkpoints[contig.Name]) == 0 {
		t.Fatal("the coverage sketch made no checkpoints")
	}
	chunks, err := src.index.Chunks(src.refs[contig.Name], start, start+100)
	if err != nil || len(chunks) == 0 {
		t.Fatalf("no index chunk covers the region: %v", err)
	}
	at, ok := src.checkpointBefore(contig.Name, start)
	if !ok || !offsetBefore(chunks[0].Begin, at.Offset) {
		t.Fatal("the region cannot start at a checkpoint")
	}
	if len(got) != len(want) || gotSeen != wantSeen {
		t.Fatalf("the checkpointed read found %d of %d reads, want %d of %d", len(got), gotSeen, len(want), wantSeen)
	}
	for i := range want {
		if got[i].Name != want[i].Name {
			t.Fatalf("read %d is %s, want %s", i, got[i].Name, want[i].Name)
		}
	}
	t.Logf("dense contig: sketch %v, region of %d reads in %v",
		sketch.Round(time.Microsecond), len(got), read.Round(time.Microsecond))
}

// TestDenseContigUsesTheSketchForAThinSample checks that a screen holding far
// more reads than the cap shows the sketch instead of drawing an unrepresentative
// sample, that the filter line says how many were sampled, and that a pan out
// there costs a profile lookup rather than a BAM scan.
func TestDenseContigUsesTheSketchForAThinSample(t *testing.T) {
	var dense screenCase
	for _, c := range screenCases {
		if c.sketch {
			dense = c
		}
	}
	m := screenModel(t, dense)
	if !m.showReads() || !m.truncated {
		t.Fatalf("the base-level view shows %d reads, sampled %v, want a sample", len(m.reads), m.truncated)
	}
	baseReads := m.readsSeen
	summary := strings.Join(m.filterSummary(), " ")
	want := fmt.Sprintf("Visualizing %s sampled reads of %s total in this region",
		comma(len(m.reads)), comma(baseReads))
	if !strings.Contains(summary, want) {
		t.Fatalf("the filter line says %q, want %q", summary, want)
	}

	m.view.bpp = 10 // a screen now holds hundreds of thousands of reads
	m.held = window{}
	m.reload()
	if m.showReads() || m.covPeak == 0 {
		t.Fatalf("a screen of %s reads still draws a pileup instead of the sketch", comma(baseReads))
	}

	begin := time.Now()
	m.pan(1000)
	took := time.Since(begin)
	t.Logf("dense sketch: a base screen holds %s reads, pan at 10 bases per column %v",
		comma(baseReads), took.Round(time.Microsecond))
	if took > 50*time.Millisecond {
		t.Fatalf("a pan in the sketch took %v", took)
	}
}

// TestDenseSampledViewFillsTheScreen checks that a sample is spread through each
// slice rather than bunched at its left edge. The first rows of a wide pileup
// then cover the whole width, so it turns into the sketch without a seam.
func TestDenseSampledViewFillsTheScreen(t *testing.T) {
	var dense screenCase
	for _, c := range screenCases {
		if c.sketch {
			dense = c
		}
	}
	m := screenModel(t, dense)
	m.view = viewState{m.view.start, 5}
	m.held = window{}
	begin := time.Now()
	m.reload()
	loaded := time.Since(begin)
	if !m.showReads() || m.readsSeen <= len(m.reads) {
		t.Fatalf("the wide view has %d of %d reads and is not sampled", len(m.reads), m.readsSeen)
	}

	m.View()
	l := m.layout()
	covered := make([]bool, l.columns)
	for row := l.pileupTop; row < min(l.pileupTop+30, l.pileupTop+l.pileupRows); row++ {
		for col := 0; col < l.columns; col++ {
			covered[col] = covered[col] || m.frame[row][gutterWidth+col].r != ' '
		}
	}
	shown := 0
	for _, ok := range covered {
		if ok {
			shown++
		}
	}
	t.Logf("dense sample: %d of %d reads at 5 bases per column, load %v, %d of %d columns in the first rows",
		len(m.reads), m.readsSeen, loaded.Round(time.Microsecond), shown, l.columns)
	if shown < l.columns*9/10 {
		t.Fatalf("the first rows cover %d of %d columns, so the sample bunches at the slice edges", shown, l.columns)
	}
	if loaded > 500*time.Millisecond {
		t.Fatalf("loading a sampled screen took %v", loaded)
	}
}

func TestPanningStaysSmooth(t *testing.T) {
	m := ready(t, "multi_contig_uneven")
	start := time.Now()
	for i := 0; i < 60; i++ {
		press(m, "right")
		m.View()
	}
	each := time.Since(start) / 60
	t.Log("one pan step, load and draw:", each.Round(time.Microsecond))
	if each > 50*time.Millisecond {
		t.Errorf("panning takes %v a step", each)
	}
}

// TestHoveringStaysFast keeps one mouse move cheap: the terminal sends one for
// every cell the pointer crosses, so anything slow here shows up as highlights
// trailing behind the mouse.
func TestHoveringStaysFast(t *testing.T) {
	m := ready(t, "viral_variants")
	for _, bpp := range []int{1, 200} { // zoomed out there are far more reads on screen
		m.view.bpp = bpp
		m.reload()
		start := time.Now()
		for i := 0; i < 200; i++ {
			mouse(m, gutterWidth+i%80, m.layout().pileupTop+i%8, tea.MouseActionMotion, tea.MouseButtonNone)
			m.View()
		}
		each := time.Since(start) / 200
		t.Logf("one mouse move at %d bases per column, hover and draw: %v", bpp, each.Round(time.Microsecond))
		if each > 2*time.Millisecond {
			t.Errorf("a mouse move at %d bases per column takes %v", bpp, each)
		}
	}
}

// TestLongReadZoomingStaysFast walks every zoom rung on a long-read file,
// where one read can span tens of thousands of bases. It records the time one
// keypress takes and pins the rule that a notch from the whole-contig view
// moves to the neighbouring rung rather than skipping one.
func TestLongReadZoomingStaysFast(t *testing.T) {
	src, err := OpenAlignments("testdata/bf_fragilis.bam")
	if err != nil {
		t.Skip("missing long-read test data")
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	m.panel, m.width, m.height = "", 150, 44
	for ci := range src.Contigs {
		if src.Contigs[ci].Length == 0 {
			continue
		}
		m.contig, m.contigIndex = ci, ci
		m.view = viewState{max(0, src.Contigs[ci].Length/2-75), 1}
		m.reload()
		fit, worst := m.fitZoom(), time.Duration(0)
		sawReads := false
		for m.view.bpp < fit {
			before := m.view.bpp
			start := time.Now()
			m.zoom(1, m.layout().columns/2)
			m.View()
			took := time.Since(start)
			worst = max(worst, took)
			if len(m.reads) > 0 {
				sawReads = true
			}
			t.Logf("%s %6d -> %6d bases per column: %v, %d reads",
				src.Contigs[ci].Name, before, m.view.bpp, took.Round(time.Microsecond), len(m.reads))
		}
		if !sawReads {
			t.Errorf("%s never showed a read while zooming", src.Contigs[ci].Name)
		}
		if worst > 100*time.Millisecond {
			t.Errorf("%s took %v for one zoom step", src.Contigs[ci].Name, worst)
		}
		before := m.view.bpp
		m.zoom(-1, m.layout().columns/2)
		if want := nextZoom(before, -1, fit); m.view.bpp != want {
			t.Errorf("%s stepped in from %d to %d, want the neighbouring rung %d",
				src.Contigs[ci].Name, before, m.view.bpp, want)
		}
	}
}

// screenCase is one shape of work the screen meets: long reads, a dense
// pileup, or a big short-read genome.
type screenCase struct {
	name, align, ref string
	contig           string // "" opens the contig newModel picks
	bpp              int
	middle           bool // start at the middle of the contig, where reads are dense
	long             bool // whether this is one of the ultra-long-read screens
	sketch           bool // build the coverage sketch and its checkpoints first
}

var screenCases = []screenCase{
	{"long-read plasmid", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "NC_006873.1", 200, false, true, false},
	{"long-read chromosome", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.fa", "NC_003228.3", 150, true, true, false},
	{"short-read control", "testdata/bf_fragilis_short_20x.bam", "testdata/bf_fragilis.fa", "NC_003228.3", 150, true, false, false},
	{"ONT 203x", "testdata/hiv1.bam", "testdata/hiv1.gb", "", 50, false, false, false},
	{"dense 300x", "testdata/influenza_h1n1.bam", "", "", 20, false, false, false},
	{"genome 40x", "testdata/ecoli_k12.bam", "testdata/ecoli_k12.gb", "", 50, false, false, false},
	{"dense mitochondrial 36,000x", "testdata/dense_mito.bam", "", "NC_013700.1", 1, true, false, true},
}

// screenModel prepares one case the way opening the file does: the panel is
// closed, the named contig is chosen, and a base-level screen is loaded so the
// model learns the file's read length before it is asked for a wide view.
func screenModel(tb testing.TB, c screenCase) *model {
	tb.Helper()
	src, err := OpenAlignments(c.align)
	if err != nil {
		tb.Skip(err)
	}
	tb.Cleanup(src.Close)
	var ref *Reference
	if c.ref != "" {
		if ref, err = LoadReference(c.ref); err != nil {
			tb.Fatal(err)
		}
		ref.MatchContigs(src.Contigs)
	}
	m := newModel(src, ref, c.ref, nil)
	m.panel, m.panelFocus = "", false
	m.width, m.height = 150, 44
	if c.contig != "" {
		found := false
		for i, contig := range src.Contigs {
			if contig.Name == c.contig {
				m.contig, m.contigIndex, found = i, i, true
				break
			}
		}
		if !found {
			tb.Fatalf("%s has no contig %s", c.align, c.contig)
		}
	}
	contig := m.contigs()[m.contig]
	if contig.Length == 0 {
		tb.Fatalf("%s chose the empty contig %q", c.name, contig.Name)
	}
	if c.sketch {
		m.Update(sketchCoverage(src)())
	}
	start := 0
	if c.middle {
		start = max(0, contig.Length/2-75*c.bpp)
	}
	m.view = viewState{start, 1}
	m.reload() // learn the read length before asking for a view that may hide reads
	m.held = window{}
	m.view = viewState{start, c.bpp}
	m.reload()
	return m
}

// BenchmarkLoad is a cached reload: the region is already held and the view
// moves only a few bases.
func BenchmarkLoad(b *testing.B) {
	for _, c := range screenCases {
		b.Run(c.name, func(b *testing.B) {
			m := screenModel(b, c)
			start := m.view.start
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				m.view = viewState{start + i*7, c.bpp}
				m.reload()
			}
		})
	}
}

// BenchmarkColdLoad asks for a region the model has not held, the way a zoom
// or a filter change does.
func BenchmarkColdLoad(b *testing.B) {
	for _, c := range screenCases {
		b.Run(c.name, func(b *testing.B) {
			m := screenModel(b, c)
			start := m.view.start
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				m.view = viewState{start, c.bpp}
				m.held = window{}
				m.reload()
			}
		})
	}
}

// TestScreenCasesShowReads guards the benchmark setup: every case must load
// reads at its measured zoom, and the pileup must not be blank. It also records
// one cold load and frame per case.
func TestScreenCasesShowReads(t *testing.T) {
	for _, c := range screenCases {
		t.Run(c.name, func(t *testing.T) {
			m := screenModel(t, c)
			if !m.showReads() {
				t.Fatalf("the view at %d bases per column does not show reads", c.bpp)
			}
			if len(m.reads) == 0 || m.rowCount == 0 {
				t.Fatalf("the view holds %d reads in %d rows", len(m.reads), m.rowCount)
			}
			lines := screen(t, m, 150, 44)
			l := m.layout()
			seen := false
			for row := l.pileupTop; row < l.pileupTop+l.pileupRows && row < len(lines); row++ {
				if strings.ContainsAny(lines[row], "ACGT▶◀━") {
					seen = true
					break
				}
			}
			if !seen {
				t.Fatal("the pileup area is blank")
			}
			m.view = viewState{m.view.start, c.bpp}
			m.held = window{}
			begin := time.Now()
			m.reload()
			loaded := time.Since(begin)
			begin = time.Now()
			m.View()
			drawn := time.Since(begin)
			t.Logf("%d reads, %d rows, learned read length %d, cold load %v, draw %v",
				len(m.reads), m.rowCount, m.readLength, loaded.Round(time.Microsecond), drawn.Round(time.Microsecond))
			if loaded+drawn > 300*time.Millisecond {
				t.Errorf("one cold screen took %v", loaded+drawn)
			}
		})
	}
}

func BenchmarkFrame(b *testing.B) {
	for _, c := range screenCases {
		b.Run(c.name, func(b *testing.B) {
			m := screenModel(b, c)
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				m.View()
			}
		})
	}
}

// The dense mitochondrial case is the one that stresses pan and zoom: a pan
// reuses rows, while the sketch takes over once a screen holds far more reads
// than the sample cap.
func BenchmarkDensePan(b *testing.B) {
	m := screenModel(b, screenCases[len(screenCases)-1])
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		m.pan(10)
		m.View()
	}
}

func BenchmarkDenseSketchPan(b *testing.B) {
	m := screenModel(b, screenCases[len(screenCases)-1])
	m.view.bpp = 10
	m.held = window{}
	m.reload()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		m.pan(100)
		m.View()
	}
}

// The two long-read screens use the real key step, with the frame that follows.
func BenchmarkLongReadPan(b *testing.B)  { benchLongRead(b, "right", "left") }
func BenchmarkLongReadZoom(b *testing.B) { benchLongRead(b, "+", "-") }

func benchLongRead(b *testing.B, forward, back string) {
	for _, c := range screenCases {
		if !c.long {
			continue
		}
		b.Run(c.name, func(b *testing.B) {
			m := screenModel(b, c)
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				if i%2 == 0 {
					press(m, forward)
				} else {
					press(m, back)
				}
				m.View()
			}
		})
	}
}
