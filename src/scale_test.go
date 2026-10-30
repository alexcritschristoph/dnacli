package main

// scale_test.go stress-tests dnacli with large synthetic fixtures made by
// /tmp/opencode/scale/make_fixtures.py: a 250 Mb contig, 50,000 long-ish
// reads (the sampling path), 1,000 contigs, 10,000 CDS and 1,000,000 variants.
//
//	go test . -run TestScale -v

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
	"time"
)

const scaleDir = "/tmp/opencode/scale"

func scalePath(name string) string { return filepath.Join(scaleDir, name) }

func needScale(t *testing.T, names ...string) {
	t.Helper()
	for _, n := range names {
		if _, err := os.Stat(scalePath(n)); err != nil {
			t.Skip("scale fixtures not built")
		}
	}
}

func heapMB() float64 {
	runtime.GC()
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	return float64(m.HeapAlloc) / (1 << 20)
}

func frameTime(m *model, w, h int) time.Duration {
	m.width, m.height = w, h
	m.reload()
	start := time.Now()
	m.View()
	return time.Since(start)
}

func TestScaleLargeContig(t *testing.T) {
	needScale(t, "large.fa", "large.bam")
	base = heapMB()

	start := time.Now()
	ref, err := LoadReference(scalePath("large.fa"))
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("250 Mb FASTA: %d contigs, %d bp, load %v, heap +%.0f MB", len(ref.Seq),
		len(ref.Seq["chrBig"]), time.Since(start), heapMB()-base)

	start = time.Now()
	src, err := OpenAlignments(scalePath("large.bam"))
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	t.Logf("large BAM: load %v, reads %d", time.Since(start), src.Contigs[0].Reads)
	ref.MatchContigs(src.Contigs)

	start = time.Now()
	m := newModel(src, ref, scalePath("large.fa"), nil)
	m.width, m.height = 150, 44
	m.reload()
	t.Logf("model: %v, fitZoom=%d bpp, heap +%.0f MB", time.Since(start), m.fitZoom(), heapMB()-base)

	start = time.Now()
	profile := src.CoverageProfile("chrBig", profileBins)
	t.Logf("coverage sketch over 250 Mb: %v (%d bins)", time.Since(start), len(profile))
	src.Contigs[0].Profile = profile
	m.reload()

	for _, step := range []string{"z", "0"} {
		press(m, step)
		t.Logf("frame at %d bp/col after %q: %v (reads=%d)", m.view.bpp, step, frameTime(m, 150, 44), len(m.reads))
	}

	press(m, "end")
	out := plain(m.View())
	t.Logf("end of a 250 Mb contig: start=%d span=%d bpp=%d", m.view.start, m.span(), m.view.bpp)
	if m.view.start+m.span() < m.contigs()[m.contig].Length {
		t.Logf("FINDING: the end key does not reach the end of a 250 Mb contig")
	}
	for _, line := range strings.Split(out, "\n") {
		if strings.Contains(line, "250,000,000") {
			t.Logf("header at the end: %s", strings.TrimSpace(line))
			break
		}
	}

	press(m, "home")
	press(m, "+")
	press(m, "+")
	t.Logf("zoom from fit with +: %d bp/col (ladder max is %d)", m.view.bpp, zoomLadder[len(zoomLadder)-1])
	if m.view.bpp == zoomLadder[len(zoomLadder)-1] && m.fitZoom() > zoomLadder[len(zoomLadder)-1] {
		t.Logf("note: pressing + right after 0 drops from %d to %d bp/col in one step", m.fitZoom(), m.view.bpp)
	}

}

func TestScaleSampling(t *testing.T) {
	needScale(t, "dense.fa", "dense.bam")
	ref, err := LoadReference(scalePath("dense.fa"))
	if err != nil {
		t.Fatal(err)
	}
	src, err := OpenAlignments(scalePath("dense.bam"))
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	ref.MatchContigs(src.Contigs)

	m := newModel(src, ref, scalePath("dense.fa"), nil)
	m.width, m.height = 150, 44
	m.reload()
	t.Logf("read length measured: %d; reads at base level: %d", m.readLength, len(m.reads))

	for _, bpp := range []int{1000, 500, 100, 10} {
		m.view = viewState{0, bpp}
		m.reload()
		dup := 0
		seenNames := map[string]int{}
		for _, a := range m.reads {
			seenNames[a.Name]++
			if seenNames[a.Name] > 1 {
				dup++
			}
		}
		d := frameTime(m, 150, 44)
		t.Logf("bpp=%4d: visible=%d reads, seen=%d, truncated=%v, duplicates=%d, frame=%v",
			bpp, len(m.reads), m.readsSeen, m.truncated, dup, d)
		if m.truncated {
			pct := 0
			if m.readsSeen > 0 {
				pct = 100 * len(m.reads) / m.readsSeen
			}
			t.Logf("  filter line percentage would read %d%% of %d; true region holds ~%d",
				pct, m.readsSeen, m.readsSeen)
			if pct > 95 {
				t.Logf("FINDING: when sampling, readsSeen counts only until the per-slice limit, so the 'x of y (p%%)' line reports ~100%% even though most reads are not drawn")
			}
			out := plain(m.View())
			for _, line := range strings.Split(out, "\n") {
				if strings.Contains(line, "Visualizing") {
					t.Logf("  filter strip says: %s", strings.TrimSpace(line))
					break
				}
			}
		}
		if dup > 0 {
			t.Logf("FINDING: sampling returned %d duplicate reads", dup)
		}
		if m.truncated && bpp == 1000 {
			m.filter.MinMapQ = 1 // a non-default filter makes the strip show x of y
			m.reload()
			out := plain(m.View())
			for _, line := range strings.Split(out, "\n") {
				if strings.Contains(line, "Visualizing") {
					t.Logf("  with a filter set the strip says: %s", strings.TrimSpace(line))
					break
				}
			}
			m.filter = defaultFilter()
		}
	}
}

func TestScaleManyContigs(t *testing.T) {
	needScale(t, "many.fa", "many.bam", "many.gb", "lastcontig.vcf")
	base = heapMB()
	start := time.Now()
	ref, err := LoadReference(scalePath("many.gb"))
	if err != nil {
		t.Fatal(err)
	}
	nfeat := 0
	for _, fs := range ref.Features {
		nfeat += len(fs)
	}
	t.Logf("GenBank with %d records / %d features: %v, heap +%.0f MB", len(ref.Seq), nfeat, time.Since(start), heapMB()-base)

	src, err := OpenAlignments(scalePath("many.bam"))
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	t.Logf("BAM with %d contigs: loaded", len(src.Contigs))
	ref.MatchContigs(src.Contigs)
	m := newModel(src, ref, scalePath("many.gb"), nil)
	m.width, m.height = 150, 44
	m.reload()
	t.Logf("start contig %q; panel=%q", m.contigs()[m.contig].Name, m.panel)

	start = time.Now()
	m.Update(sketchCoverage(src)())
	t.Logf("coverage sketch of 1000 contigs (one open): %v", time.Since(start))
	if d := time.Since(start); d > 2*time.Second {
		t.Fatalf("sketching 1000 contigs took %v", d)
	}

	start = time.Now()
	m.panel, m.panelFocus = "contigs", true
	m.selectContig(999)
	t.Logf("select last of 1000 contigs: %v, now %q", time.Since(start), m.contigs()[m.contig].Name)
	t.Logf("frame with 1000-contig panel: %v", frameTime(m, 150, 44))

	start = time.Now()
	genes := m.genes("")
	t.Logf("genes(\"\") across 1000 contigs: %d genes, %v", len(genes), time.Since(start))
	start = time.Now()
	m.command = "gene "
	suggest := m.suggestions()
	t.Logf("gene suggestions: %d shown, %v", len(suggest), time.Since(start))
	m.command = ""

	variants, err := LoadVariants(scalePath("lastcontig.vcf"))
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("lastcontig.vcf: %d variants", len(variants))
	start = time.Now()
	err = variantsMatch(m.contigs(), variants)
	t.Logf("variantsMatch(1000 contigs, 100k variants): %v, err=%v", time.Since(start), err)
	if d := time.Since(start); d > 500*time.Millisecond {
		t.Logf("FINDING: openVariants validates each variant by scanning every contig (O(variants×contigs)); %v for 100k×1000", d)
	}
}

func TestScaleVariants(t *testing.T) {
	needScale(t, "million.vcf", "variants.fa", "variants.gb")
	base = heapMB()

	start := time.Now()
	vs, err := LoadVariants(scalePath("million.vcf"))
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("1,000,000-record VCF: load %v, heap +%.0f MB (%.1f bytes/variant)",
		time.Since(start), heapMB()-base, (heapMB()-base)*(1<<20)/float64(len(vs)))

	ref, err := LoadReference(scalePath("variants.gb"))
	if err != nil {
		t.Fatal(err)
	}
	nfeat := 0
	for _, fs := range ref.Features {
		nfeat += len(fs)
	}
	t.Logf("reference: %d features", nfeat)

	// an alignment on the same contig
	dir := t.TempDir()
	sam := filepath.Join(dir, "v.sam")
	read := strings.Repeat("ACGT", 13)[:50]
	qual := strings.Repeat("I", 50)
	body := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrVar\tLN:1000000\n" +
		"r1\t0\tchrVar\t100\t60\t50M\t*\t0\t0\t" + read + "\t" + qual + "\n"
	os.WriteFile(sam, []byte(body), 0o644)
	src, err := OpenAlignments(sam)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	ref.MatchContigs(src.Contigs)

	start = time.Now()
	m := newModel(src, ref, scalePath("variants.gb"), vs)
	m.width, m.height = 150, 44
	m.reload()
	t.Logf("newModel with 1,000,000 variants (cached CDS positions): %v, heap +%.0f MB",
		time.Since(start), heapMB()-base)
	if d := time.Since(start); d > 2*time.Second {
		t.Logf("FINDING: opening a 1M-variant VCF costs %v before the first frame", d)
	}

	t.Logf("frame at 1 bp/col with 1M variants: %v", frameTime(m, 150, 44))
	press(m, "0")
	t.Logf("frame at fit with 1M variants: %v", frameTime(m, 150, 44))

	start = time.Now()
	m.showVariant(len(m.variants) - 1)
	t.Logf("select variant 1,000,000: %v, start=%d bpp=%d", time.Since(start), m.view.start, m.view.bpp)
	t.Logf("frame at the last variant: %v", frameTime(m, 150, 44))

	press(m, "-")
	press(m, "-")
	t.Logf("frame after zooming out: %v", frameTime(m, 150, 44))

	// marker drawing walks every variant every frame; compare with no variants
	m.variants = nil
	m.reload()
	t.Logf("frame with no variants: %v", frameTime(m, 150, 44))
}

var base float64

// TestScaleSamplingIsOnePass checks the sampler on 50,000 reads over 100 kb:
// it must keep a share from the whole region, draw no read twice, and say how
// many the region really held. Reading each slice of the region separately
// used to decompress the same BAM blocks over and over: 734 ms for this view.
func TestScaleSamplingIsOnePass(t *testing.T) {
	needScale(t, "dense.bam")
	src, err := OpenAlignments(scalePath("dense.bam"))
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	m.width, m.height = 150, 44
	m.view = viewState{0, 1} // measures the read length, which decides how far out reads still show
	m.reload()
	m.view = viewState{0, 400}
	start := time.Now()
	m.reload()
	took := time.Since(start)
	t.Logf("sampled view of 50,000 reads: %d drawn of %d seen in %v", len(m.reads), m.readsSeen, took)
	if !m.truncated {
		t.Fatalf("this view holds %d reads and was not sampled", len(m.reads))
	}
	if len(m.reads) > maxAlignments {
		t.Fatalf("the sample kept %d reads, more than the %d it may hold", len(m.reads), maxAlignments)
	}
	if m.readsSeen <= len(m.reads) {
		t.Fatalf("the region reported %d reads but drew %d, so the count stops at the sample", m.readsSeen, len(m.reads))
	}
	names, halves := map[string]bool{}, [2]int{}
	length := src.Contigs[0].Length
	for _, a := range m.reads {
		if names[a.Name] {
			t.Fatalf("read %s was drawn twice", a.Name)
		}
		names[a.Name] = true
		halves[max(0, a.Start*2/length)]++
	}
	if halves[0] == 0 || halves[1] == 0 {
		t.Fatalf("the sample covers %d reads in the first half and %d in the second", halves[0], halves[1])
	}
	if took > 500*time.Millisecond {
		t.Fatalf("sampling this view took %v", took)
	}
}
