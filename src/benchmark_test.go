package main

import (
	"fmt"
	"testing"

	tea "github.com/charmbracelet/bubbletea"
)

// The benchmarks below measure the work behind one user action: a key press, a
// mouse move, a command, a resize, or a whole startup. They always use the same
// window size, datasets and zoom levels, so two runs can be compared directly.
var (
	benchDatasets  = []string{"viral_variants", "supertest", "multi_contig_uneven"}
	benchZooms     = []int{1, 20, 500} // 500 leaves the reads for the depth sketch
	benchReadZooms = []int{1, 20}      // the pileup is only drawn this close in
)

// loadTestModel opens one test dataset the way main does at startup.
func loadTestModel(stem string) (*model, error) {
	path := "testdata/" + stem + ".bam"
	if stem == "paired_reads" {
		path = "testdata/" + stem + ".sam"
	}
	src, err := OpenAlignments(path)
	if err != nil {
		return nil, err
	}
	refPath := "testdata/" + stem + ".gb"
	ref, err := LoadReference(refPath)
	if err != nil {
		refPath = "testdata/" + stem + ".fa"
		ref, _ = LoadReference(refPath)
	}
	ref.MatchContigs(src.Contigs)
	variants, _ := LoadVariants("testdata/" + stem + ".vcf")
	return newModel(src, ref, refPath, variants), nil
}

// benchReady loads a dataset and frames the first screen, outside the timer.
func benchReady(b *testing.B, stem string, bpp int) *model {
	b.Helper()
	m, err := loadTestModel(stem)
	if err != nil {
		b.Fatal(err)
	}
	b.Cleanup(func() { m.src.Close() })
	m.width, m.height = 150, 44
	m.view = viewState{0, bpp}
	m.reload()
	m.View() // fill the shared color and escape caches before timing
	b.ReportAllocs()
	b.ResetTimer()
	return m
}

// benchScreens runs one action across datasets and zoom levels.
func benchScreens(b *testing.B, datasets []string, zooms []int, action func(*model, *testing.B)) {
	for _, stem := range datasets {
		for _, bpp := range zooms {
			b.Run(fmt.Sprintf("%s/%d", stem, bpp), func(b *testing.B) {
				m := benchReady(b, stem, bpp)
				action(m, b)
			})
		}
	}
}

// BenchmarkStartup is opening the files, working out the variants and building
// the first screen.
func BenchmarkStartup(b *testing.B) {
	for _, stem := range benchDatasets {
		b.Run(stem, func(b *testing.B) {
			if warm, err := loadTestModel(stem); err == nil {
				warm.src.Close()
			}
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				m, err := loadTestModel(stem)
				if err != nil {
					b.Fatal(err)
				}
				m.src.Close()
			}
		})
	}
}

// BenchmarkMeasureImpacts is the one-off cost of annotating every variant,
// with the CDS position cache warm.
func BenchmarkMeasureImpacts(b *testing.B) {
	ref, err := LoadReference("testdata/bf_fragilis.gb")
	if err != nil {
		b.Fatal(err)
	}
	variants, err := LoadVariants("testdata/bf_fragilis_multisample.vcf")
	if err != nil {
		b.Fatal(err)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		m := &model{ref: ref, variants: variants}
		m.measureImpacts()
	}
}

// BenchmarkReload is the work a pan, a zoom or a filter change does before the
// next frame is drawn.
func BenchmarkReload(b *testing.B) {
	benchScreens(b, benchDatasets, benchZooms, func(m *model, b *testing.B) {
		for i := 0; i < b.N; i++ {
			m.reload()
		}
	})
}

// BenchmarkDraw is one frame, with the reads already loaded.
func BenchmarkDraw(b *testing.B) {
	benchScreens(b, benchDatasets, benchZooms, func(m *model, b *testing.B) {
		for i := 0; i < b.N; i++ {
			m.View()
		}
	})
}

func BenchmarkPan(b *testing.B) {
	benchScreens(b, benchDatasets, benchReadZooms, func(m *model, b *testing.B) {
		for i := 0; i < b.N; i++ {
			if i%2 == 0 {
				press(m, "right")
			} else {
				press(m, "left")
			}
			m.View()
		}
	})
}

func BenchmarkZoom(b *testing.B) {
	benchScreens(b, benchDatasets, benchReadZooms, func(m *model, b *testing.B) {
		for i := 0; i < b.N; i++ {
			if i%2 == 0 {
				press(m, "+")
			} else {
				press(m, "-")
			}
			m.View()
		}
	})
}

// BenchmarkHover is one mouse move: the terminal sends one for every cell the
// pointer crosses, so the highlight has to keep up.
func BenchmarkHover(b *testing.B) {
	benchScreens(b, benchDatasets, benchZooms, func(m *model, b *testing.B) {
		l := m.layout()
		for i := 0; i < b.N; i++ {
			x := gutterWidth + (i*13)%l.columns
			y := l.pileupTop + (i*3)%max(1, l.pileupRows)
			mouse(m, x, y, tea.MouseActionMotion, tea.MouseButtonNone)
			m.View()
		}
	})
}

func BenchmarkScrollReads(b *testing.B) {
	m := benchReady(b, "multi_contig_uneven", 20)
	for i := 0; i < b.N; i++ {
		if i%2 == 0 {
			press(m, "down")
		} else {
			press(m, "up")
		}
		m.View()
	}
}

// BenchmarkCommandBar is one keystroke while a filter command is being typed,
// with its list of suggestions redrawn.
func BenchmarkCommandBar(b *testing.B) {
	m := benchReady(b, "viral_variants", 20)
	text := "filters mapq 30"
	for i := 0; i < b.N; i++ {
		if i%len(text) == 0 {
			m.command, m.choice = "", 0
		}
		press(m, string(text[i%len(text)]))
		m.View()
	}
}

// BenchmarkFilterCommand is a whole typed filter change, reload included.
func BenchmarkFilterCommand(b *testing.B) {
	m := benchReady(b, "supertest", 20)
	for i := 0; i < b.N; i++ {
		m.filter = defaultFilter()
		m.command, m.choice = "", 0
		run(m, "filters mapq 30")
		m.View()
	}
}

// BenchmarkFilterArrow is changing the chosen filter with the arrow keys.
func BenchmarkFilterArrow(b *testing.B) {
	m := benchReady(b, "supertest", 20)
	b.StopTimer()
	m.openPanel("filters")
	m.panelFocus, m.filterIndex = true, 0
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		if i%2 == 0 {
			press(m, "right")
		} else {
			press(m, "left")
		}
		m.View()
	}
}

// BenchmarkPanel opens and closes a panel, redrawing the narrower or wider
// pileup each time.
func BenchmarkPanel(b *testing.B) {
	for _, name := range []string{"variants", "contigs", "colors", "filters"} {
		b.Run(name, func(b *testing.B) {
			m := benchReady(b, "multi_contig_uneven", 20)
			for i := 0; i < b.N; i++ {
				m.togglePanel(name)
				m.View()
			}
		})
	}
}

// BenchmarkColorSchemes changes the coloring, which some schemes have to
// measure the reads for before the next frame.
func BenchmarkColorSchemes(b *testing.B) {
	for i, c := range colors {
		b.Run(c, func(b *testing.B) {
			m := benchReady(b, "viral_variants", 20)
			b.StopTimer()
			m.color = i
			m.reload()
			m.View()
			b.StartTimer()
			for n := 0; n < b.N; n++ {
				m.color = i
				m.reload()
				m.View()
			}
		})
	}
}

func BenchmarkContigSwitch(b *testing.B) {
	m := benchReady(b, "multi_contig_uneven", 1)
	for i := 0; i < b.N; i++ {
		m.selectContig(i % len(m.contigs()))
		m.View()
	}
}

func BenchmarkVariantJump(b *testing.B) {
	m := benchReady(b, "viral_variants", 20)
	for i := 0; i < b.N; i++ {
		m.showVariant(i % len(m.variants))
		m.View()
	}
}

// BenchmarkResize is a change of window size, which reflows every row.
func BenchmarkResize(b *testing.B) {
	m := benchReady(b, "viral_variants", 20)
	for i := 0; i < b.N; i++ {
		if i%2 == 0 {
			m.width, m.height = 150, 44
		} else {
			m.width, m.height = 100, 30
		}
		m.reload()
		m.View()
	}
}
