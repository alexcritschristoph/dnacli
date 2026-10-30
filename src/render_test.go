package main

import (
	"os"
	"regexp"
	"strings"
	"testing"
	"time"
)

var ansi = regexp.MustCompile("\x1b\\[[0-9;]*m")

// screen renders the whole interface at a fixed size, the way the terminal
// would draw it, and returns the plain text.
func screen(t *testing.T, m *model, width, height int) []string {
	t.Helper()
	m.width, m.height = width, height
	m.reload()
	start := time.Now()
	out := m.View()
	if elapsed := time.Since(start); elapsed > 250*time.Millisecond {
		t.Errorf("drawing one frame took %v", elapsed)
	}
	var lines []string
	for _, line := range strings.Split(out, "\n") {
		lines = append(lines, strings.TrimRight(ansi.ReplaceAllString(line, ""), " "))
	}
	if len(lines) != height {
		t.Fatalf("drew %d lines, want %d", len(lines), height)
	}
	return lines
}

func testModel(t *testing.T, stem string) *model {
	t.Helper()
	path := "testdata/" + stem + ".bam"
	if stem == "paired_reads" {
		path = "testdata/" + stem + ".sam"
	}
	src, err := OpenAlignments(path)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(src.Close)
	refPath := "testdata/" + stem + ".gb"
	ref, err := LoadReference(refPath)
	if err != nil {
		refPath = "testdata/" + stem + ".fa"
		ref, _ = LoadReference(refPath)
	}
	ref.MatchContigs(src.Contigs)
	variants, _ := LoadVariants("testdata/" + stem + ".vcf")
	return newModel(src, ref, refPath, variants)
}

func TestScreenshots(t *testing.T) {
	dir := os.Getenv("SHOTS")
	if dir == "" {
		t.Skip("set SHOTS to a directory to save screenshots")
	}
	shots := map[string]func() *model{
		"variants": func() *model {
			m := testModel(t, "viral_variants")
			m.showVariant(0)
			return m
		},
		"zoomed-out": func() *model {
			m := testModel(t, "viral_variants")
			m.view = viewState{0, 20}
			return m
		},
		"contigs": func() *model {
			m := testModel(t, "multi_contig_uneven")
			m.panel = "contigs"
			return m
		},
		"supertest": func() *model {
			m := testModel(t, "supertest")
			m.panel = "filters"
			return m
		},
		"pairs": func() *model {
			m := testModel(t, "paired_reads")
			m.view = viewState{0, 5}
			return m
		},
		"help": func() *model {
			m := testModel(t, "plasmid_even")
			m.panel = "help"
			return m
		},
	}
	for name, build := range shots {
		m := build()
		m.width, m.height = 150, 44
		m.reload()
		if err := os.WriteFile(dir+"/"+name+".txt", []byte(m.View()+"\n"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
}

func TestIndelsAreDrawn(t *testing.T) {
	m := ready(t, "viral_variants")
	for _, at := range []int{2, 3} { // the deletion at 3,100 and the insertion at 4,300
		m.showVariant(at)
		lines := screen(t, m, 150, 44)
		found := false
		for _, line := range lines[10:] {
			if strings.Contains(line, "---") {
				found = true
			}
		}
		insertion := false
		for _, a := range m.reads {
			insertion = insertion || len(a.Insertions) > 0
		}
		t.Logf("%s deletion gaps drawn: %v, reads carrying an insertion: %v", m.variants[at].Label(), found, insertion)
		if at == 2 && !found {
			t.Error("no deleted bases are drawn at the deletion")
		}
		if at == 3 && !insertion {
			t.Error("no read carries the inserted bases")
		}
	}
}

func TestReverseStrandProtein(t *testing.T) {
	m := ready(t, "viral_variants")
	m.showVariant(7) // the deletion inside the reverse-strand gene nucleo
	lines := screen(t, m, 150, 44)
	t.Log("reference", lines[5])
	t.Log("protein  ", lines[6])
	impact := m.ref.CodingImpact("VIRUS1", 6402, "C", "A")
	t.Log("a substitution inside nucleo:", impact)
	if impact.Gene != "nucleo" || impact.ToAmino == 0 {
		t.Fatalf("the reverse-strand gene gave %v", impact)
	}
}

// TestRulerSitsAboveTheGenes pins down the order of the bars at the top:
// the ruler, then every gene row, then the variants.
func TestRulerSitsAboveTheGenes(t *testing.T) {
	m := ready(t, "viral_variants")
	m.view = viewState{250, 1} // the first gene starts at 300
	screen(t, m, 150, 44)      // draws one frame and times it
	l := m.layout()
	ruler, genes, variants := rowWith(m, "position"), rowWith(m, "genes"), rowWith(m, "variants")
	if l.geneRows < 1 {
		t.Fatal("the example drew no gene rows")
	}
	if ruler != trackTop {
		t.Fatalf("the ruler is on row %d, want %d", ruler, trackTop)
	}
	if genes != ruler+1 {
		t.Fatalf("the genes are on row %d, want %d", genes, ruler+1)
	}
	last := genes + l.geneRows - 1
	if variants != last+1 || variants != l.refRow-1 {
		t.Fatalf("the variants are on row %d, want %d, just below the last gene row %d", variants, last+1, last)
	}
}
