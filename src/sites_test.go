package main

import (
	"fmt"
	"strings"
	"testing"
	"time"
)

// TestVariantMarkers pins down when a position is marked from the reads and
// which colour each level of support gets.
func TestVariantMarkers(t *testing.T) {
	ref := testReference(t, genBankRecord("MARK", strings.Repeat("A", 40), ""))

	cases := []struct {
		depth, support int
		marked         bool
		ratio          float64
	}{
		{40, 2, false, 0},    // two reads are not enough, however deep the site
		{40, 3, true, 0.075}, // three reads of a deep site
		{10, 3, true, 0.3},   // the 5% floor is passed
		{4, 2, false, 0},     // half the reads, but only two of them
		{3, 3, true, 1},      // every read differs
	}
	for _, c := range cases {
		counts, depth := oneSite(c.depth, 'C', c.support)
		sites := callSites(counts, depth, ref, "MARK", 0)
		s, found := sites[0]
		if found != c.marked {
			t.Fatalf("%d of %d reads marked %v, want %v", c.support, c.depth, found, c.marked)
		}
		if found && (s.Alt != 'C' || s.Depth != c.depth || fmt.Sprintf("%.3f", s.Ratio) != fmt.Sprintf("%.3f", c.ratio)) {
			t.Fatalf("%d of %d reads gave %+v", c.support, c.depth, s)
		}
	}

	// A base that matches the reference is never a variant, however common.
	counts, depth := oneSite(10, 'A', 6)
	if sites := callSites(counts, depth, ref, "MARK", 0); len(sites) != 0 {
		t.Fatalf("the reference base was called as a variant: %+v", sites)
	}

	// The marker row colors by support: grey to 30%, orange to 50%, red above.
	m := ready(t, "viral_variants")
	m.view = viewState{0, 1}
	m.reload()
	m.sites = map[int]site{
		10: {Depth: 40, Ratio: 0.1},
		20: {Depth: 40, Ratio: 0.3},
		30: {Depth: 40, Ratio: 0.6},
	}
	m.View()
	row := rowWith(m, "variants")
	if row < 0 {
		t.Fatal("no variants row was drawn")
	}
	for pos, want := range map[int]string{10: colorMuted, 20: colorWarn, 30: colorAlert} {
		cell := m.frame[row][gutterWidth+m.column(pos)]
		if cell.r != '▪' || cell.s.fg != want {
			t.Fatalf("a site with %.0f%% support drew %q in %s, want a square in %s",
				100*m.sites[pos].Ratio, cell.r, cell.s.fg, want)
		}
	}
	if cell := m.frame[row][gutterWidth+m.column(15)]; cell.r != ' ' {
		t.Fatalf("a column without a site drew %q", cell.r)
	}
}

// oneSite is a pileup of one reference position, with a second base seen a
// given number of times.
func oneSite(depth int, alt byte, support int) ([][4]int, []int) {
	counts := make([][4]int, 1)
	counts[0][0] = depth - support
	counts[0][baseIndex(alt)] = support
	return counts, []int{depth}
}

// TestVariantMarkerPerformance times the call over a whole 8 kb contig at a
// depth of 30, which is the most a single screen ever asks it to do.
func TestVariantMarkerPerformance(t *testing.T) {
	ref := testReference(t, genBankRecord("MARK", strings.Repeat("ACGT", 2000), ""))
	counts := make([][4]int, 8000)
	depth := make([]int, 8000)
	for i := range counts {
		counts[i][0], counts[i][1], counts[i][3] = 28, 1, 1
		depth[i] = 30
	}
	start := time.Now()
	sites := callSites(counts, depth, ref, "MARK", 0)
	elapsed := time.Since(start)
	t.Logf("8,000 positions at depth 30: %d sites in %s", len(sites), elapsed.Round(time.Microsecond))
	if elapsed > 5*time.Millisecond {
		t.Fatalf("calling sites over a contig took %v", elapsed)
	}
}
