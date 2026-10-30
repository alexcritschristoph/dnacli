package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// These tests pin the fixes made after the audit, so a later change that
// brings one back shows up in the ordinary test run.

func TestFastaHeaderAndInnerWhitespace(t *testing.T) {
	dir := t.TempDir()
	bad := filepath.Join(dir, "bad.fa")
	os.WriteFile(bad, []byte(">\nACGT\n"), 0o644)
	if _, err := LoadReference(bad); err == nil {
		t.Fatal("a FASTA record with no name did not give an error")
	}

	spaced := filepath.Join(dir, "spaced.fa")
	os.WriteFile(spaced, []byte(">c\nAC GT\nA\tC\n"), 0o644)
	ref, err := LoadReference(spaced)
	if err != nil {
		t.Fatal(err)
	}
	if got := string(ref.Seq["c"]); got != "ACGTAC" {
		t.Fatalf("spaced FASTA loaded %q, want ACGTAC", got)
	}
}

func TestCommandBarTakesPositionPunctuation(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "")
	m.width, m.height = 150, 44
	for _, text := range []string{"go NC_003228.3:5000-6000", "gene BF9343_RS00005"} {
		m.command = ""
		typeText(m, text)
		if m.command != text {
			t.Fatalf("typing %q left the bar holding %q", text, m.command)
		}
		press(m, "enter")
		if m.command != "" || m.view.bpp == 1 || strings.HasPrefix(m.message, "No ") {
			t.Fatalf("%q did not run: command %q, bpp %d, message %q", text, m.command, m.view.bpp, m.message)
		}
	}
}

func TestHeaderlessSamIsExplained(t *testing.T) {
	_, err := OpenAlignments("testdata/bf_fragilis_noheader.sam")
	if err == nil || !strings.Contains(err.Error(), "@SQ") {
		t.Fatalf("a SAM without @SQ gave %v", err)
	}
}

func TestWrappedGenBankQualifierKeepsAllItsText(t *testing.T) {
	p := filepath.Join(t.TempDir(), "wrapped.gb")
	body := "LOCUS       c 12 bp DNA linear BCT 01-JAN-2020\n" +
		"FEATURES             Location/Qualifiers\n" +
		"     CDS             1..9\n" +
		"                     /product=\"a long product name that\n" +
		"                     wraps onto the next line\"\n" +
		"ORIGIN\n        1 acgtacgtacgt\n//\n"
	os.WriteFile(p, []byte(body), 0o644)
	ref, err := LoadReference(p)
	if err != nil {
		t.Fatal(err)
	}
	if got := ref.Features["c"][0].Product; got != "a long product name that wraps onto the next line" {
		t.Fatalf("wrapped product loaded %q", got)
	}
}

func TestExportLeavesOverlaysOut(t *testing.T) {
	m := ready(t, "viral_variants")
	m.width, m.height = 150, 44
	m.panel = "help"
	m.reload()
	m.View()
	p := filepath.Join(t.TempDir(), "help.svg")
	if err := m.writeSVG(p, darkSVG); err != nil {
		t.Fatal(err)
	}
	svg, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(svg), "Navigation") {
		t.Fatal("the exported picture contains the help overlay")
	}
}

func TestMissingBaseQualitiesAreNeutral(t *testing.T) {
	a := &Alignment{Start: 0, End: 1, Bases: []byte{'A'}, Quals: []byte{0xff}}
	scheme := namedIndex(colors, "Base quality")
	if got := readColor(scheme, a, 0, 150); got != quietGrey {
		t.Fatalf("a base with no quality is coloured %s, want the neutral %s", got, quietGrey)
	}
	a.Quals[0] = 40
	if got := readColor(scheme, a, 0, 150); got == quietGrey {
		t.Fatal("a Q40 base is coloured as neutral")
	}
}

func TestVariantImpactsAreMeasured(t *testing.T) {
	start := time.Now()
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis_multisample.vcf")
	took := time.Since(start)
	if m.variants[0].Impact.Kind == "" {
		t.Fatal("the variants were not measured against the reference")
	}
	t.Logf("21,690 variants loaded and measured in %v", took)
}

func TestZoomedOutReadsKeepOnlyTheirDifferences(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "")
	m.width, m.height = 150, 44
	for i, c := range m.contigs() {
		if c.Length < 100000 {
			m.contig, m.contigIndex = i, i
		}
	}
	m.view = viewState{0, m.fitZoom()} // the whole small contig, still a read view
	m.reload()
	if !m.showReads() || len(m.reads) == 0 {
		t.Skip("the small contig is not in a read view")
	}
	differences := 0
	for _, a := range m.reads {
		if a.Bases != nil {
			t.Fatalf("read %s at %d bases per column still keeps a base for every position", a.Name, m.view.bpp)
		}
		for i, at := range a.MismatchAt {
			if got := a.BaseAt(at); got != a.MismatchBase[i] {
				t.Fatalf("read %s shows %q at %s, want %q", a.Name, got, comma(at+1), a.MismatchBase[i])
			}
			differences++
		}
	}
	if differences == 0 {
		t.Fatal("no read had a difference to check")
	}

	m.view = viewState{0, 1}
	m.reload()
	if len(m.reads) == 0 {
		t.Skip("no reads at one base per column")
	}
	for _, a := range m.reads {
		if a.Bases == nil {
			t.Fatalf("read %s at one base per column lost its bases", a.Name)
		}
	}
}

func TestReadsGiveWayToTheDepthTrack(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "")
	m.width, m.height = 150, 44
	small := 0
	for i, c := range m.contigs() {
		if c.Length < 100000 {
			small = i
		}
	}
	m.contig, m.contigIndex = small, small
	m.view = viewState{0, m.fitZoom()}
	m.reload()
	if !m.showReads() {
		t.Fatalf("the %s contig lost its reads when shown whole", m.contigs()[small].Name)
	}

	big := 1 - small
	m.contig, m.contigIndex = big, big
	m.view = viewState{0, 500}
	m.reload()
	if m.showReads() || len(m.reads) != 0 {
		t.Fatalf("the %s contig still loads reads at 500 bases per column", m.contigs()[big].Name)
	}
}

func TestCoarseAndDetailedReadsSpanTheSameBases(t *testing.T) {
	dir := t.TempDir()
	fa := filepath.Join(dir, "ref.fa")
	os.WriteFile(fa, []byte(">c\n"+strings.Repeat("A", 100)+"\n"), 0o644)
	ref, err := LoadReference(fa)
	if err != nil {
		t.Fatal(err)
	}
	sam := filepath.Join(dir, "reads.sam")
	body := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\n" +
		"leading\t0\tc\t1\t60\t3D10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n" +
		"middle\t0\tc\t31\t60\t10M\t*\t0\t0\tAAAAACAAAA\tIIIIIIIIII\n" +
		"trailing\t0\tc\t61\t60\t10M3D\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n"
	os.WriteFile(sam, []byte(body), 0o644)
	src, err := OpenAlignments(sam)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	full, _, _ := src.Reads("c", 0, 100, defaultFilter(), ref, true)
	coarse, _, _ := src.Reads("c", 0, 100, defaultFilter(), ref, false)
	if len(full) != len(coarse) {
		t.Fatalf("%d detailed reads and %d coarse ones", len(full), len(coarse))
	}
	for i := range full {
		if full[i].Start != coarse[i].Start || full[i].End != coarse[i].End {
			t.Fatalf("%s spans %d–%d in detail and %d–%d coarse",
				full[i].Name, full[i].Start, full[i].End, coarse[i].Start, coarse[i].End)
		}
		for j, at := range coarse[i].MismatchAt {
			if got, want := coarse[i].BaseAt(at), full[i].BaseAt(at); got != want {
				t.Fatalf("%s coarse shows %q at %s, detailed %q", full[i].Name, got, comma(at+1), want)
			}
			if coarse[i].MismatchBase[j] != full[i].BaseAt(at) {
				t.Fatalf("%s kept %q as its difference", full[i].Name, coarse[i].MismatchBase[j])
			}
		}
	}
}

func TestWheelAtTheZoomLimitDoesNotPan(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "")
	m.width, m.height = 150, 44
	m.contig, m.contigIndex = 0, 0
	m.view = viewState{100000, 1}
	m.reload()
	if len(m.inserts) == 0 {
		t.Skip("this region has no insertions")
	}
	anchor, start := 90, m.view.start
	for i := 0; i < 20; i++ {
		m.zoom(-1, anchor) // a trackpad's momentum keeps sending these
	}
	if m.view.start != start {
		t.Fatalf("twenty wheel-in events at one base per column panned from %d to %d", start, m.view.start)
	}

	m.view = viewState{0, m.fitZoom()}
	m.reload()
	start = m.view.start
	for i := 0; i < 20; i++ {
		m.zoom(1, anchor)
	}
	if m.view.start != start {
		t.Fatalf("twenty wheel-out events at the whole contig panned from %d to %d", start, m.view.start)
	}
}

func TestOnlyTheReadsOnScreenLookForConsensusDifferences(t *testing.T) {
	m := auditOpen(t, "testdata/influenza_h1n1.bam", "", "") // no reference, so the consensus is the guide
	m.width, m.height = 150, 44
	m.view = viewState{0, 20}
	m.reload()
	m.View()

	visible := 0
	for _, a := range m.reads {
		if a.Row >= m.rowOffset && a.Row < m.rowOffset+m.layout().pileupRows {
			visible++
		}
	}
	computed := 0
	for _, a := range m.reads {
		if a.consensusDone {
			computed++
		}
	}
	if computed != visible {
		t.Fatalf("%d of %d reads were compared against the consensus, but %d are on screen", computed, len(m.reads), visible)
	}
	if computed == 0 || computed == len(m.reads) {
		t.Fatalf("%d of %d reads were compared, so the work is not limited to those drawn", computed, len(m.reads))
	}

	m.scroll(5) // the rows that come into view are compared as they are drawn
	m.View()
	scrolled := 0
	for _, a := range m.reads {
		if a.consensusDone {
			scrolled++
		}
	}
	if scrolled <= computed {
		t.Fatalf("scrolling showed %d more rows but no more reads were compared (%d before, %d after)", 5, computed, scrolled)
	}
}
