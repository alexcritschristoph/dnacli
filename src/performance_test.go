package main

import (
	"fmt"
	"os"
	"path/filepath"
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
	m.width, m.height = 150, 44
	for ci := range src.Contigs {
		m.contig, m.contigIndex = ci, ci
		m.view = viewState{0, 1}
		m.reload()
		fit, worst := m.fitZoom(), time.Duration(0)
		for m.view.bpp < fit {
			before := m.view.bpp
			start := time.Now()
			m.zoom(1, m.layout().columns/2)
			m.View()
			took := time.Since(start)
			worst = max(worst, took)
			t.Logf("%s %6d -> %6d bases per column: %v, %d reads",
				src.Contigs[ci].Name, before, m.view.bpp, took.Round(time.Microsecond), len(m.reads))
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

// The benchmarks below time the two things a keystroke costs: loading the
// region behind the screen, and drawing it. Each case is a shape of work the
// screen meets: long reads, a dense pileup, and a big short-read genome.
var screenCases = []struct {
	name, align, ref string
	plasmid, bpp     int
}{
	{"long-read plasmid", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", 1, 200},
	{"ONT 203x", "testdata/hiv1.bam", "testdata/hiv1.gb", 0, 50},
	{"dense 300x", "testdata/influenza_h1n1.bam", "", 0, 20},
	{"genome 40x", "testdata/ecoli_k12.bam", "testdata/ecoli_k12.gb", 0, 50},
}

func screenModel(b *testing.B, c struct {
	name, align, ref string
	plasmid, bpp     int
}) *model {
	src, err := OpenAlignments(c.align)
	if err != nil {
		b.Skip(err)
	}
	b.Cleanup(src.Close)
	var ref *Reference
	if c.ref != "" {
		if ref, err = LoadReference(c.ref); err != nil {
			b.Fatal(err)
		}
		ref.MatchContigs(src.Contigs)
	}
	m := newModel(src, ref, c.ref, nil)
	m.width, m.height = 150, 44
	if c.plasmid == 1 {
		for i, cg := range src.Contigs {
			if cg.Length < 100000 {
				m.contig, m.contigIndex = i, i
			}
		}
	}
	m.view = viewState{0, c.bpp}
	m.reload()
	return m
}

func BenchmarkLoad(b *testing.B) {
	for _, c := range screenCases {
		b.Run(c.name, func(b *testing.B) {
			m := screenModel(b, c)
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				m.view = viewState{i * 7, c.bpp}
				m.reload()
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
