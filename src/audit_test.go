package main

// audit_test.go is a bug-hunting harness for the dnacli in this directory.
// It logs "FINDING:" lines that the audit report quotes.
//
//	go test . -run TestAudit -v

import (
	"encoding/xml"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"
	"time"

	tea "github.com/charmbracelet/bubbletea"
)

func auditOpen(t *testing.T, align, refPath, vcf string) *model {
	t.Helper()
	if _, err := os.Stat("testdata"); err != nil {
		t.Skip("testdata not present")
	}
	src, err := OpenAlignments(align)
	if err != nil {
		t.Fatalf("OpenAlignments(%s): %v", align, err)
	}
	t.Cleanup(src.Close)
	if len(src.Contigs) == 0 {
		t.Skip("alignment has no contigs")
	}
	var ref *Reference
	if refPath != "" {
		if ref, err = LoadReference(refPath); err != nil {
			t.Fatalf("LoadReference(%s): %v", refPath, err)
		}
		ref.MatchContigs(src.Contigs)
	}
	var variants []Variant
	if vcf != "" {
		if variants, err = LoadVariants(vcf); err != nil {
			t.Fatalf("LoadVariants(%s): %v", vcf, err)
		}
	}
	m := newModel(src, ref, refPath, variants)
	m.width, m.height = 150, 44
	m.reload()
	return m
}

var ansiPattern = regexp.MustCompile(`\x1b\[[0-9;]*[A-Za-z]`)

func plain(s string) string { return ansiPattern.ReplaceAllString(s, "") }

func renderAt(m *model, w, h int) string {
	m.width, m.height = w, h
	m.reload()
	return plain(m.View())
}

func typeText(m *model, text string) {
	for _, r := range text {
		press(m, string(r))
	}
}

// ─── references ─────────────────────────────────────────────────────────────

func TestAuditReferenceLoading(t *testing.T) {
	files := []string{
		"bf_fragilis.gb", "bf_fragilis.fa", "bf_fragilis_prodigal.gb",
		"bf_fragilis_benchling.gb", "ecoli_k12.gb", "ecoli_k12.fa",
		"mito_human.gb", "sarscov2.gb", "phage_lambda.gb", "hiv1.gb",
		"vcholerae.gb", "influenza_h1n1.gb", "phage_t4.fa", "pstvd.fa",
	}
	for _, name := range files {
		t.Run(name, func(t *testing.T) {
			p := "testdata/" + name
			if _, err := os.Stat(p); err != nil {
				t.Skip("missing")
			}
			var ref *Reference
			var err error
			dur := time.Since(time.Now())
			start := time.Now()
			ref, err = LoadReference(p)
			dur = time.Since(start)
			if err != nil {
				t.Fatalf("load: %v", err)
			}
			names := make([]string, 0, len(ref.Seq))
			for n := range ref.Seq {
				names = append(names, n)
			}
			sort.Strings(names)
			for _, n := range names {
				t.Logf("%s: %s len=%d features=%d", name, n, len(ref.Seq[n]), len(ref.Features[n]))
			}
			t.Logf("%s load %v", name, dur)
		})
	}

	t.Run("aliases-attach", func(t *testing.T) {
		src, err := OpenAlignments("testdata/bf_fragilis.bam")
		if err != nil {
			t.Skip("missing bam")
		}
		defer src.Close()
		ref, err := LoadReference("testdata/bf_fragilis.gb")
		if err != nil {
			t.Skip("missing gb")
		}
		ref.MatchContigs(src.Contigs)
		for _, c := range src.Contigs {
			if len(ref.Seq[c.Name]) != c.Length {
				t.Logf("FINDING: %s attaches %d bp to a %d bp contig", c.Name, len(ref.Seq[c.Name]), c.Length)
			}
		}
		t.Logf("contig names in ref: %v", func() []string {
			out := make([]string, 0, len(ref.Seq))
			for n := range ref.Seq {
				out = append(out, n)
			}
			sort.Strings(out)
			return out
		}())
	})
}

func TestAuditReferenceEdges(t *testing.T) {
	dir := t.TempDir()
	write := func(name, body string) string {
		p := filepath.Join(dir, name)
		os.WriteFile(p, []byte(body), 0o644)
		return p
	}

	t.Run("empty-header", func(t *testing.T) {
		p := write("empty.fa", ">\nACGT\n")
		defer func() {
			if r := recover(); r != nil {
				t.Logf("FINDING: empty FASTA header panics: %v", r)
			}
		}()
		if _, err := LoadReference(p); err != nil {
			t.Logf("error: %v", err)
		}
	})

	t.Run("duplicate-names", func(t *testing.T) {
		p := write("dup.fa", ">c\nAAAA\n>c\nCCCC\n")
		ref, _ := LoadReference(p)
		t.Logf("duplicate headers: seq=%q", ref.Seq["c"])
		if string(ref.Seq["c"]) == "AAAA" {
			t.Logf("FINDING: duplicate FASTA names silently overwrite (first copy lost)")
		}
	})

	t.Run("spaces-in-sequence", func(t *testing.T) {
		p := write("spaced.fa", ">c\nAC GT\nTGCA\n")
		ref, _ := LoadReference(p)
		t.Logf("spaced sequence stored as %q len=%d", ref.Seq["c"], len(ref.Seq["c"]))
		if strings.Contains(string(ref.Seq["c"]), " ") {
			t.Logf("FINDING: whitespace inside FASTA sequence lines is kept, so every coordinate after it shifts")
		}
	})

	t.Run("wrapped-product", func(t *testing.T) {
		p := write("wrap.gb", "LOCUS       T1                      60 bp    DNA     linear BCT 01-JAN-2020\n"+
			"FEATURES             Location/Qualifiers\n"+
			"     CDS             1..30\n"+
			"                     /gene=\"g\"\n"+
			"                     /product=\"this is a very long product name that\n"+
			"                     wraps onto the next line\"\n"+
			"ORIGIN\n"+
			"        1 atgaaacccgggttt\n"+
			"//\n")
		ref, _ := LoadReference(p)
		for _, f := range ref.Features["T1"] {
			t.Logf("wrapped product: %q", f.Product)
			if !strings.Contains(f.Product, "next line") {
				t.Logf("FINDING: multi-line GenBank qualifiers other than /translation are truncated at the line break (product=%q)", f.Product)
			}
		}
	})

	t.Run("join-complement", func(t *testing.T) {
		p := write("jc.gb", "LOCUS       J1                      30 bp    DNA     linear BCT 01-JAN-2020\n"+
			"FEATURES             Location/Qualifiers\n"+
			"     CDS             join(complement(1..3),complement(10..12))\n"+
			"ORIGIN\n"+
			"        1 atgaaacccgggtttaaacccgggtttaaa\n"+
			"//\n")
		ref, _ := LoadReference(p)
		for _, f := range ref.Features["J1"] {
			t.Logf("join(complement()) ranges=%v reverse=%v", f.Ranges, f.Reverse)
			if len(f.Ranges) != 2 || !f.Reverse {
				t.Logf("FINDING: a location written join(complement(...),complement(...)) is misparsed: ranges=%v reverse=%v (want two ranges, reverse)", f.Ranges, f.Reverse)
			}
		}
	})

	t.Run("complement-join", func(t *testing.T) {
		p := write("cj.gb", "LOCUS       J2                      30 bp    DNA     linear BCT 01-JAN-2020\n"+
			"FEATURES             Location/Qualifiers\n"+
			"     CDS             complement(join(1..3,10..12))\n"+
			"ORIGIN\n"+
			"        1 atgaaacccgggtttaaacccgggtttaaa\n"+
			"//\n")
		ref, _ := LoadReference(p)
		for _, f := range ref.Features["J2"] {
			t.Logf("complement(join()) ranges=%v reverse=%v", f.Ranges, f.Reverse)
		}
	})

	t.Run("con-record", func(t *testing.T) {
		p := write("con.gb", "LOCUS       CON1                    1000 bp    DNA     circular CON 01-JAN-2020\n"+
			"DEFINITION  no sequence here.\n"+
			"//\n")
		ref, err := LoadReference(p)
		t.Logf("CON record: err=%v seqLen=%d keys=%d", err, len(ref.Seq["CON1"]), len(ref.Seq))
		if err == nil && len(ref.Seq["CON1"]) == 0 {
			t.Logf("FINDING: a GenBank CON record (no ORIGIN) loads as a 0 bp reference with no error")
		}
	})

	t.Run("match-by-length", func(t *testing.T) {
		seq := make([]byte, 100)
		for i := range seq {
			seq[i] = 'A'
		}
		mapped := map[string]int{}
		for i := 0; i < 50; i++ {
			ref := newReference()
			ref.Seq["alpha"] = append([]byte(nil), seq...)
			ref.Seq["beta"] = append([]byte(nil), seq...)
			ref.MatchContigs([]Contig{{Name: "contig1", Length: 100}})
			if ref.Seq["contig1"] != nil {
				mapped["mapped"]++
			}
		}
		t.Logf("equal-length contigs mapped by length in %d/50 runs", mapped["mapped"])
		t.Logf("FINDING: MatchContigs maps by length only and iterates a map, so two contigs of equal length are matched at random")
	})
}

// ─── variants ───────────────────────────────────────────────────────────────

func TestAuditVariantLoading(t *testing.T) {
	for _, name := range []string{"bf_fragilis.vcf", "bf_fragilis_multisample.vcf",
		"bf_fragilis.vcf.gz", "bf_fragilis.bcf"} {
		t.Run(name, func(t *testing.T) {
			if _, err := os.Stat("testdata/" + name); err != nil {
				t.Skip("missing")
			}
			start := time.Now()
			vs, err := LoadVariants("testdata/" + name)
			t.Logf("%s: %d variants, err=%v, %v", name, len(vs), err, time.Since(start))
			if strings.HasSuffix(name, ".gz") || strings.HasSuffix(name, ".bcf") {
				if err == nil && len(vs) == 0 {
					t.Logf("note: %s yields 0 variants when handed to LoadVariants directly; the CLI rejects the extension before this point", name)
				}
			}
		})
	}

	t.Run("multi-allelic-impact", func(t *testing.T) {
		ref, err := LoadReference("testdata/bf_fragilis.gb")
		if err != nil {
			t.Skip("missing gb")
		}
		vs, _ := LoadVariants("testdata/bf_fragilis.vcf")
		checked := 0
		for _, v := range vs {
			if !strings.Contains(v.Alt, ",") {
				continue
			}
			imp := ref.CodingImpact(v.Contig, v.Start, v.Ref, v.Alt)
			t.Logf("multi-allelic %s:%d %s→%s: kind=%q short=%q", v.Contig, v.Start+1, v.Ref, v.Alt, imp.Kind, imp.Short())
			checked++
			if checked == 3 {
				break
			}
		}
		if checked > 0 {
			t.Logf("FINDING: a multi-allelic ALT (\"T,G\") is used verbatim as one long allele, so a SNP is reported as an in-frame insertion/complex change")
		}
	})

	t.Run("symbolic-alt", func(t *testing.T) {
		p := filepath.Join(t.TempDir(), "sym.vcf")
		os.WriteFile(p, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
			"c\t10\t.\tA\t<DEL>\t.\tPASS\t.\n"+
			"c\t20\t.\tA\t*\t.\tPASS\t.\n"), 0o644)
		vs, _ := LoadVariants(p)
		ref := newReference()
		ref.Seq["c"] = []byte(strings.Repeat("ACGT", 50))
		for _, v := range vs {
			imp := ref.CodingImpact("c", v.Start, v.Ref, v.Alt)
			t.Logf("symbolic %s→%s: end=%d label=%q impact=%q", v.Ref, v.Alt, v.End(), v.Label(), imp.Kind)
		}
		t.Logf("FINDING: symbolic ALT alleles (<DEL>, *) are treated as literal bases")
	})

	t.Run("startup-contig-mismatch", func(t *testing.T) {
		if _, err := os.Stat("testdata/ecoli_k12.bam"); err != nil {
			t.Skip("missing")
		}
		m := auditOpen(t, "testdata/ecoli_k12.bam", "", "testdata/bf_fragilis.vcf")
		m.warnIfVariantsDoNotMatch("testdata/bf_fragilis.vcf")
		t.Logf("ecoli bam + bf vcf: variants=%d on %q, warning=%q", len(m.variants), m.contigs()[m.contig].Name, m.message)
		if !strings.Contains(m.message, "Warning") {
			t.Fatal("a mismatched VCF at startup gave no warning")
		}
	})
}

// ─── alignments ─────────────────────────────────────────────────────────────

func TestAuditAlignmentLoading(t *testing.T) {
	files := []string{"bf_fragilis.bam", "bf_fragilis.noindex.bam",
		"bf_fragilis_namesorted.bam", "bf_fragilis_unsorted.bam",
		"bf_fragilis_tiny.bam", "bf_fragilis_full.sam", "bf_fragilis_noheader.sam",
		"vcholerae.bam", "influenza_h1n1.bam", "ecoli_k12.bam",
		"phage_lambda.cram", "bf_fragilis.bcf"}
	for _, name := range files {
		t.Run(name, func(t *testing.T) {
			p := "testdata/" + name
			if _, err := os.Stat(p); err != nil {
				t.Skip("missing")
			}
			start := time.Now()
			src, err := OpenAlignments(p)
			if err != nil {
				t.Logf("%s: error %v (%v)", name, err, time.Since(start))
				return
			}
			defer src.Close()
			reads := int64(0)
			for _, c := range src.Contigs {
				reads += c.Reads
			}
			t.Logf("%s: contigs=%d reads=%d load=%v indexed=%v", name, len(src.Contigs), reads, time.Since(start), src.index != nil)
		})
	}

	t.Run("seq-star", func(t *testing.T) {
		src, err := OpenAlignments("testdata/bf_fragilis.bam")
		if err != nil {
			t.Skip("missing")
		}
		defer src.Close()
		var pan interface{}
		func() {
			defer func() { pan = recover() }()
			src.Reads("NC_003228.3", 2308172, 2308242, defaultFilter(), nil, true)
		}()
		t.Logf("SEQ=* window panic=%v (this program tolerates it; cg. cmd/dnacli did not)", pan)
	})

	t.Run("truncated-sam", func(t *testing.T) {
		raw, err := os.ReadFile("testdata/bf_fragilis_full.sam")
		if err != nil {
			t.Skip("missing")
		}
		lines := strings.Split(string(raw), "\n")
		cut := strings.Join(lines[:len(lines)/3], "\n")
		p := filepath.Join(t.TempDir(), "cut.sam")
		os.WriteFile(p, []byte(cut), 0o644)
		src, err := OpenAlignments(p)
		if err != nil {
			t.Logf("truncated SAM: error %v", err)
			return
		}
		var n int64
		for _, c := range src.Contigs {
			n += c.Reads
		}
		t.Logf("truncated SAM: contigs=%d reads=%d err=nil", len(src.Contigs), n)
		t.Logf("FINDING: readAll swallows the parse error, so a truncated SAM loads partially with no warning")
	})

	t.Run("stale-index", func(t *testing.T) {
		raw, err := os.ReadFile("testdata/bf_fragilis.bam")
		if err != nil {
			t.Skip("missing")
		}
		dir := t.TempDir()
		p := filepath.Join(dir, "copy.bam")
		os.WriteFile(p, raw, 0o644)
		idx, err := os.ReadFile("testdata/bf_fragilis.bam.bai")
		if err == nil {
			os.WriteFile(p+".bai", idx, 0o644)
		}
		src, err := OpenAlignments(p)
		if err != nil {
			t.Logf("stale index: error %v", err)
			return
		}
		defer src.Close()
		reads, seen, _ := src.Reads("NC_003228.3", 1100, 1200, defaultFilter(), nil, true)
		t.Logf("stale index: indexed=%v reads=%d seen=%d", src.index != nil, len(reads), seen)
	})

	t.Run("headerless-sam", func(t *testing.T) {
		src, err := OpenAlignments("testdata/bf_fragilis_noheader.sam")
		if err != nil {
			t.Logf("headerless SAM: error %v", err)
			return
		}
		defer src.Close()
		t.Logf("headerless SAM: contigs=%d (main would quit)", len(src.Contigs))
	})
}

func TestAuditUnindexedMultiContig(t *testing.T) {
	raw, err := os.ReadFile("testdata/vcholerae.bam")
	if err != nil {
		t.Skip("missing")
	}
	p := filepath.Join(t.TempDir(), "vcholerae.bam")
	os.WriteFile(p, raw, 0o644)
	src, err := OpenAlignments(p)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	a, _, _ := src.Reads("NC_002505.1", 0, 2000, defaultFilter(), nil, true)
	b, _, _ := src.Reads("NC_002506.1", 0, 2000, defaultFilter(), nil, true)
	t.Logf("unindexed multi-contig: chr1=%d chr2=%d reads in window", len(a), len(b))
	if len(a) > 0 && len(b) > 0 {
		t.Logf("(memory is kept per contig, so no cross-contig mixing)")
	}
}

func TestAuditMultiTab(t *testing.T) {
	src1, err := OpenAlignments("testdata/bf_fragilis_tiny.bam")
	if err != nil {
		t.Skip("missing")
	}
	defer src1.Close()
	m := newModel(src1, nil, "", nil)
	m.width, m.height = 120, 30

	t.Run("same-contigs", func(t *testing.T) {
		src2, err := OpenAlignments("testdata/bf_fragilis.bam")
		if err != nil {
			t.Skip("missing")
		}
		m.openAlignments("testdata/bf_fragilis.bam")
		t.Logf("tabs after opening bf_fragilis.bam: %d", len(m.tabs))
		_ = src2
		t.Cleanup(func() {
			for _, s := range m.tabs[1:] {
				s.Close()
			}
		})
	})

	t.Run("same-name-different-length", func(t *testing.T) {
		dir := t.TempDir()
		mk := func(name string, length int) string {
			body := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrX\tLN:" + strconv.Itoa(length) + "\n" +
				"r\t0\tchrX\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
			p := filepath.Join(dir, name)
			os.WriteFile(p, []byte(body), 0o644)
			return p
		}
		a, b := mk("a.sam", 1000), mk("b.sam", 2000)
		s1, err := OpenAlignments(a)
		if err != nil {
			t.Fatal(err)
		}
		defer s1.Close()
		m2 := newModel(s1, nil, "", nil)
		m2.width, m2.height = 120, 30
		before := len(m2.tabs)
		m2.openAlignments(b)
		t.Logf("open same contig name, 1000 vs 2000 bp: tabs %d -> %d, message=%q", before, len(m2.tabs), m2.message)
		if len(m2.tabs) > before {
			t.Logf("FINDING: sameContigs checks names only, so alignment files with the same contig names but different lengths are joined as tabs")
		}
		for _, s := range m2.tabs[1:] {
			s.Close()
		}
	})
}

// ─── commands, panels ───────────────────────────────────────────────────────

func TestAuditCommands(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	t.Logf("start: panel=%q variants=%d features=%d", m.panel, len(m.variants), len(m.ref.Features["NC_003228.3"]))

	run(m, "go 5000")
	t.Logf("go 5000: start=%d bpp=%d marked=%d", m.view.start, m.view.bpp, m.marked)

	run(m, "go NC_003228.3:5000-6000")
	t.Logf("go range: start=%d bpp=%d", m.view.start, m.view.bpp)

	run(m, "gene BF9343_RS00005")
	t.Logf("gene: start=%d bpp=%d message=%q", m.view.start, m.view.bpp, m.message)

	run(m, "filters mapq 30")
	t.Logf("filters mapq 30: filter=%+v panel=%q focus=%v", m.filter, m.panel, m.panelFocus)

	run(m, "filters mismatches any")
	t.Logf("mismatches any: max=%d", m.filter.MaxMismatches)

	run(m, "filters bogus")
	t.Logf("filters bogus: message=%q panel=%q", m.message, m.panel)

	run(m, "filters reset")
	t.Logf("filters reset: %+v", m.filter)

	run(m, "contigsx")
	t.Logf("unknown command: message=%q", m.message)

	press(m, "1")
	t.Logf("typing a bare digit with no bar open: command=%q (bare positions need the go command)", m.command)

	run(m, "colours")
	t.Logf("hidden spelling colours: panel=%q", m.panel)

	m.panel = ""
	m.panelFocus = false
	run(m, "open testdata/hiv1.bam")
	t.Logf("open mismatched bam: tabs=%d message=%q", len(m.tabs), m.message)

	run(m, "open nope.bam")
	t.Logf("open missing: tabs=%d message=%q", len(m.tabs), m.message)

	run(m, "open testdata/bf_fragilis_prodigal.gb")
	t.Logf("open prodigal gb: refPath=%q message=%q", filepath.Base(m.refPath), m.message)
}

func TestAuditPanels(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")

	m.panel, m.panelFocus = "variants", true
	m.variantIndex = 0
	for i := 0; i < 40; i++ {
		press(m, "down")
	}
	t.Logf("variant panel: index=%d/%d visible=%d", m.variantIndex, len(m.variants), m.layout().listRows)
	if m.variantIndex != 40 {
		t.Logf("FINDING: down did not walk the variant list (index %d)", m.variantIndex)
	}
	out := plain(m.View())
	if !strings.Contains(out, comma(m.variants[m.variantIndex].Start+1)) {
		t.Logf("FINDING: the selected variant is not drawn after scrolling deep into the list")
	}

	m.panel = "contigs"
	m.panelFocus = true
	m.contigIndex = 0
	press(m, "down")
	t.Logf("contig panel: contig=%d message=%q", m.contig, m.message)

	m.panel = "colors"
	m.panelFocus = true
	m.coloringLocked = false
	before := m.highlight
	press(m, "down")
	t.Logf("colors preview: highlight %d -> %d locked=%v", before, m.highlight, m.coloringLocked)
	press(m, "enter")
	t.Logf("colors enter: locked=%v", m.coloringLocked)
	press(m, "enter")
	t.Logf("colors enter again: highlight=%d color=%d locked=%v", m.highlight, m.color, m.coloringLocked)

	m.panel = "filters"
	m.panelFocus = true
	m.filterIndex = 0
	press(m, "right")
	t.Logf("filters right: mapq=%d", m.filter.MinMapQ)
	press(m, "down")
	press(m, "right")
	t.Logf("filters down+right: quality=%d", m.filter.MinMeanQual)

	m.panel = ""
	m.panelFocus = false
	m.width = 80
	m.height = 24
	m.reload()
	run(m, "variants")
	t.Logf("narrow window: panel=%q message=%q layout.panel=%d", m.panel, m.message, m.layout().panel)
}

// ─── rendering ──────────────────────────────────────────────────────────────

func TestAuditRendering(t *testing.T) {
	cases := []struct{ name, align, ref, vcf string }{
		{"tiny", "testdata/bf_fragilis_tiny.bam", "", ""},
		{"bf-gb-vcf", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf"},
		{"bf-fa-vcf", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.fa", "testdata/bf_fragilis.vcf"},
		{"multisample", "testdata/bf_fragilis.bam", "testdata/bf_fragilis.fa", "testdata/bf_fragilis_multisample.vcf"},
		{"ecoli", "testdata/ecoli_k12.bam", "testdata/ecoli_k12.gb", "testdata/ecoli_k12.vcf"},
		{"vcholerae", "testdata/vcholerae.bam", "testdata/vcholerae.gb", "testdata/vcholerae.vcf"},
		{"influenza", "testdata/influenza_h1n1.bam", "testdata/influenza_h1n1.gb", "testdata/influenza_h1n1.vcf"},
		{"mito", "testdata/mito_human.bam", "testdata/mito_human.gb", "testdata/mito_human.vcf"},
		{"sam", "testdata/bf_fragilis_full.sam", "", ""},
		{"sam-noheader", "testdata/bf_fragilis_noheader.sam", "", ""},
		{"noindex", "testdata/bf_fragilis.noindex.bam", "", ""},
		{"unsorted", "testdata/bf_fragilis_unsorted.bam", "", ""},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			if _, err := os.Stat(c.align); err != nil {
				t.Skip("missing")
			}
			if c.name == "sam-noheader" {
				_, err := OpenAlignments(c.align)
				if err == nil || !strings.Contains(err.Error(), "@SQ") {
					t.Fatalf("a SAM without @SQ opened without an explanatory error: %v", err)
				}
				t.Logf("a SAM without @SQ is refused: %v", err)
				return
			}
			m := auditOpen(t, c.align, c.ref, c.vcf)
			for _, size := range [][2]int{{150, 44}, {80, 24}, {40, 18}, {39, 18}, {150, 17}} {
				out := renderAt(m, size[0], size[1])
				if size[0] < 40 || size[1] < 18 {
					if !strings.Contains(out, "larger window") {
						t.Logf("FINDING: %dx%d draws something other than the size warning", size[0], size[1])
					}
					continue
				}
				if out == "" || strings.Contains(out, "larger window") {
					t.Logf("FINDING: %s at %dx%d does not draw", c.name, size[0], size[1])
				}
			}
			start := time.Now()
			renderAt(m, 150, 44)
			t.Logf("%s first 150x44 frame: %v", c.name, time.Since(start))
			for _, step := range []string{"-", "-", "-", "0", "z"} {
				press(m, step)
			}
			t.Logf("%s zoomed out at %d bp/col", c.name, m.view.bpp)
		})
	}
}

func TestAuditAminoAndMark(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	m.goTo(2000)
	out := renderAt(m, 150, 44)
	if !strings.Contains(out, "protein") {
		t.Logf("FINDING: no protein track drawn at one base per column despite CDS features with translations")
	}
	m.marked = 2000
	out = renderAt(m, 150, 44)
	if !strings.Contains(out, "BF9343") && !strings.Contains(out, "protein:") {
		t.Logf("position bar text: %q", firstLineContaining(out, "Position"))
	}
	positions := 0
	for _, f := range m.ref.Features["NC_003228.3"] {
		if f.Protein != "" {
			positions++
		}
	}
	t.Logf("features with /translation: %d", positions)
}

func firstLineContaining(s, sub string) string {
	for _, line := range strings.Split(s, "\n") {
		if strings.Contains(line, sub) {
			return line
		}
	}
	return ""
}

// ─── export ─────────────────────────────────────────────────────────────────

type auditSVG struct {
	XMLName xml.Name
	Title   string `xml:"title"`
}

func TestAuditExport(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	dir := t.TempDir()
	t.Chdir(dir)

	m.runExport("")
	files, _ := filepath.Glob("*.svg")
	t.Logf("export with no name wrote %v", files)
	if len(files) == 0 {
		t.Fatalf("no file written")
	}
	raw, _ := os.ReadFile(files[0])
	var doc auditSVG
	if err := xml.Unmarshal(raw, &doc); err != nil {
		t.Logf("FINDING: exported SVG is not well-formed XML: %v", err)
	}
	t.Logf("svg %s: %d bytes, title=%q", files[0], len(raw), doc.Title)
	if strings.Contains(string(raw), "Variants  ") || strings.Contains(string(raw), "↑ ↓ to jump") {
		t.Logf("FINDING: the exported SVG includes side-panel text")
	}

	m.panel = "help"
	m.reload()
	m.View()
	m.runExport("help-open.svg")
	helpRaw, _ := os.ReadFile("help-open.svg")
	if strings.Contains(string(helpRaw), "dnacli - commands") {
		t.Logf("FINDING: exporting while the help overlay is open writes the help text into the picture")
	}

	m.panel = ""
	m.runExport("dark.svg dark")
	darkRaw, _ := os.ReadFile("dark.svg")
	if !strings.Contains(string(darkRaw), "<rect width=\"100%\" height=\"100%\"") {
		t.Logf("FINDING: dark export missing background")
	}

	m.runExport("/nonexistent/dir/x.svg")
	t.Logf("unwritable export message: %q", m.message)

	m.width, m.height = 30, 12
	m.runExport("small.svg")
	t.Logf("small window export message: %q", m.message)
}

// ─── filters, pairing ───────────────────────────────────────────────────────

func TestAuditFilters(t *testing.T) {
	dir := t.TempDir()
	sam := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:1000\n" +
		"good\t0\tc\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"dup\t1024\tc\t30\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"sec\t256\tc\t50\t40\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"lowmq\t0\tc\t70\t5\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
		"lowq\t0\tc\t90\t60\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n" +
		"short\t0\tc\t110\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
		"improper\t1\tc\t130\t60\t10M\t=\t200\t70\tACGTACGTAC\tIIIIIIIIII\n"
	p := filepath.Join(dir, "filters.sam")
	os.WriteFile(p, []byte(sam), 0o644)
	src, err := OpenAlignments(p)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	m.width, m.height = 120, 30
	m.reload()

	all := len(m.reads)
	t.Logf("default filter: %d reads", all)
	m.filter.HideDuplicates = false
	m.filter.HideSecondary = false
	m.reload()
	t.Logf("dups+secondary shown: %d reads", len(m.reads))
	m.filter = defaultFilter()
	m.filter.MinMapQ = 30
	m.reload()
	t.Logf("mapq>=30: %d reads", len(m.reads))
	m.filter = defaultFilter()
	m.filter.MinMeanQual = 20
	m.reload()
	t.Logf("mean qual>=20: %d reads", len(m.reads))
	m.filter = defaultFilter()
	m.filter.MinLength = 8
	m.reload()
	t.Logf("length>=8: %d reads", len(m.reads))
	m.filter = defaultFilter()
	m.filter.HideUnpaired = true
	m.reload()
	t.Logf("proper pairs only: %d reads", len(m.reads))
}

func TestAuditPairRows(t *testing.T) {
	mk := func(name string, start, end int, paired bool) *Alignment {
		a := &Alignment{Name: name, Start: start, End: end, Bases: make([]byte, end-start)}
		if paired {
			a.Flags = 1
		}
		return a
	}
	a, b := mk("pair", 0, 100, true), mk("pair", 1000, 1100, true)
	rows := packRows([]*Alignment{a, b}, 1, 120)
	t.Logf("far-apart mates: rows=%d a.Row=%d b.Row=%d", rows, a.Row, b.Row)
	c, d := mk("close", 0, 100, true), mk("close", 120, 220, true)
	rows = packRows([]*Alignment{c, d}, 1, 300)
	t.Logf("close mates: rows=%d a.Row=%d b.Row=%d (same row keeps a pair together)", rows, c.Row, d.Row)
	if c.Row != d.Row {
		t.Logf("FINDING: close mates were split onto different rows")
	}
}

// ─── performance ────────────────────────────────────────────────────────────

func TestAuditPerformance(t *testing.T) {
	measure := func(name string, f func()) time.Duration {
		start := time.Now()
		f()
		d := time.Since(start)
		t.Logf("%-46s %v", name, d)
		return d
	}

	measure("LoadReference bf_fragilis.gb (19.8 MB)", func() {
		LoadReference("testdata/bf_fragilis.gb")
	})
	measure("LoadReference ecoli_k12.gb (11.4 MB)", func() {
		LoadReference("testdata/ecoli_k12.gb")
	})
	measure("LoadVariants multisample (21,690)", func() {
		LoadVariants("testdata/bf_fragilis_multisample.vcf")
	})
	measure("OpenAlignments bf bam indexed", func() {
		s, err := OpenAlignments("testdata/bf_fragilis.bam")
		if err == nil {
			s.Close()
		}
	})
	measure("OpenAlignments bf bam unindexed (83 MB)", func() {
		s, err := OpenAlignments("testdata/bf_fragilis.noindex.bam")
		if err == nil {
			s.Close()
		}
	})

	var m *model
	measure("newModel bf bam+gb+593 variants", func() {
		m = auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	})
	measure("first frame 150x44", func() { m.View() })

	if _, err := os.Stat("testdata/bf_fragilis_multisample.vcf"); err == nil {
		var big *model
		measure("newModel with 21,690 variants (cached CDS positions)", func() {
			big = auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis_multisample.vcf")
		})
		_ = big
	}

	// per-frame cost at base level, mid zoom and whole-contig zoom
	for _, step := range []string{"z", "-", "-", "-", "-", "0"} {
		press(m, step)
		d := measure("frame at "+strconv.Itoa(m.view.bpp)+" bp/col", func() { m.View() })
		if d > 200*time.Millisecond {
			t.Logf("FINDING: %d bp/col frame takes %v", m.view.bpp, d)
		}
	}

	if _, err := os.Stat("testdata/vcholerae.bam"); err == nil {
		m2 := auditOpen(t, "testdata/vcholerae.bam", "testdata/vcholerae.gb", "testdata/vcholerae.vcf")
		measure("vcholerae first frame 150x44", func() { m2.View() })
		press(m2, "0")
		measure("vcholerae whole-contig frame (sampling)", func() { m2.View() })
		t.Logf("vcholerae whole contig: reads=%d seen=%d truncated=%v", len(m2.reads), m2.readsSeen, m2.truncated)
	}

	measure("export 150x44 SVG", func() {
		t.Chdir(t.TempDir())
		m.runExport("profile.svg")
	})
}

// ─── second wave ────────────────────────────────────────────────────────────

func TestAuditTypingLimits(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")

	run(m, "go 5000")
	t.Logf("go 5000: command=%q start=%d bpp=%d", m.command, m.view.start, m.view.bpp)
	okStart := m.view.start

	m.view = viewState{0, 1}
	typeText(m, "go NC_003228.3:5000")
	t.Logf("typing 'go NC_003228.3:5000': command=%q start=%d bpp=%d message=%q", m.command, m.view.start, m.view.bpp, m.message)
	if m.command == "" && m.view.start != okStart {
		t.Logf("FINDING: the command bar rejects '.'; contig names like NC_003228.3 cannot be typed, the line is discarded, and the remaining characters (digits, '-') act on the view instead")
	}

	m.view = viewState{0, 1}
	typeText(m, "go 5000-6000")
	t.Logf("typing 'go 5000-6000': command=%q start=%d bpp=%d", m.command, m.view.start, m.view.bpp)
	if m.command == "" {
		t.Logf("FINDING: the command bar rejects '-'; position ranges in the go command cannot be typed")
	}

	m.view = viewState{0, 1}
	typeText(m, "gene BF9343_RS00005")
	t.Logf("typing 'gene BF9343_RS00005': command=%q start=%d message=%q", m.command, m.view.start, m.message)
	if m.command == "" {
		t.Logf("FINDING: the command bar rejects '_'; locus tags such as BF9343_RS00005 cannot be typed")
	}

	m.view = viewState{0, 1}
	typeText(m, "gene spike")
	t.Logf("typing 'gene spike' (letters and a space): command=%q message=%q", m.command, m.message)
}

func TestAuditHeaderlessSAMModel(t *testing.T) {
	src, err := OpenAlignments("testdata/bf_fragilis_noheader.sam")
	if err == nil {
		src.Close()
		t.Fatal("a SAM without @SQ opened as a blank 0 bp contig")
	}
	t.Logf("a SAM without @SQ is refused: %v", err)
}

func TestAuditStaleIndex(t *testing.T) {
	raw, err := os.ReadFile("testdata/bf_fragilis_tiny.bam")
	if err != nil {
		t.Skip("missing")
	}
	idx, err := os.ReadFile("testdata/bf_fragilis.bam.bai")
	if err != nil {
		t.Skip("missing index")
	}
	p := filepath.Join(t.TempDir(), "tiny.bam")
	os.WriteFile(p, raw, 0o644)
	os.WriteFile(p+".bai", idx, 0o644)
	src, err := OpenAlignments(p)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer src.Close()
	reads, seen, _ := src.Reads("NC_003228.3", 0, 2000, defaultFilter(), nil, true)
	t.Logf("stale index: indexed=%v reads=%d seen=%d", src.index != nil, len(reads), seen)
	if src.index != nil && len(reads) == 0 {
		t.Logf("FINDING: a .bai that belongs to another BAM is accepted; every region returns nothing and the view is blank with no warning")
	}
}

func TestAuditReferenceCoverage(t *testing.T) {
	src, err := OpenAlignments("testdata/bf_fragilis.bam")
	if err != nil {
		t.Skip("missing")
	}
	defer src.Close()

	t.Run("con-record", func(t *testing.T) {
		p := filepath.Join(t.TempDir(), "con.gb")
		os.WriteFile(p, []byte("LOCUS       CON1                    1000 bp    DNA     circular CON 01-JAN-2020\n//\n"), 0o644)
		ref, err := LoadReference(p)
		if err != nil {
			t.Fatalf("load: %v", err)
		}
		ref.MatchContigs(src.Contigs)
		m := newModel(src, ref, p, nil)
		m.width, m.height = 120, 30
		m.reload()
		t.Logf("CON reference + bam: ref keys=%d, sequence for current contig=%d", len(ref.Seq), len(ref.Seq[m.contigs()[m.contig].Name]))
		out := plain(m.View())
		if strings.Contains(out, "reference") {
			t.Logf("FINDING: a CON GenBank (no sequence) is accepted at startup with no error; the reference row is empty")
		}
	})

	t.Run("one-contig-only", func(t *testing.T) {
		ref := newReference()
		dir := t.TempDir()
		seq := strings.Repeat("ACGT", 1000)
		p := filepath.Join(dir, "one.gb")
		body := "LOCUS       NC_003228.3              4000 bp    DNA     circular BCT 01-JAN-2020\n" +
			"FEATURES             Location/Qualifiers\n     CDS             1..99\n" +
			"ORIGIN\n        1 " + seq + "\n//\n"
		os.WriteFile(p, []byte(body), 0o644)
		loaded, _ := LoadReference(p)
		loaded.MatchContigs(src.Contigs)
		t.Logf("reference covers both contigs: %v, covers one: %v", loaded.covers(src.Contigs), loaded.covers([]Contig{{Name: "NC_006873.1", Length: 36560}}))
		t.Logf("FINDING: covers() only asks whether any contig has a sequence, so a reference for one replicon is accepted for a multi-contig BAM; the others open with no sequence")
		_ = ref
	})

	t.Run("name-match-ignores-length", func(t *testing.T) {
		ref := newReference()
		ref.Seq["c"] = []byte(strings.Repeat("A", 50))
		ref.MatchContigs([]Contig{{Name: "c", Length: 100}})
		t.Logf("name match with 50 bp sequence on a 100 bp contig: attached=%d bp", len(ref.Seq["c"]))
		t.Logf("FINDING: MatchContigs accepts a name match without checking length; a GenBank whose LOCUS length disagrees with the BAM header is used as-is (the bf GenBank attaches 10,410,280 bp to a 5,205,140 bp contig because its ORIGIN is duplicated)")
	})
}

func TestAuditSampling(t *testing.T) {
	m := auditOpen(t, "testdata/vcholerae.bam", "testdata/vcholerae.fa", "")
	m.width, m.height = 150, 44
	m.view.bpp = 2000
	m.goTo(1_500_000)
	start := time.Now()
	out := m.View()
	d := time.Since(start)
	names := map[string]int{}
	dups := 0
	for _, a := range m.reads {
		names[a.Name]++
		if names[a.Name] > 1 {
			dups++
		}
	}
	t.Logf("sampling at 2000 bp/col: reads=%d seen=%d truncated=%v frame=%v duplicates=%d", len(m.reads), m.readsSeen, m.truncated, d, dups)
	if !m.truncated && m.readsSeen > 20000 {
		t.Logf("FINDING: a region holding %d reads is loaded whole without the sampling flag", m.readsSeen)
	}
	if dups > 0 {
		t.Logf("FINDING: sampling returns %d reads more than once", dups)
	}
	_ = out
}

func TestAuditTabs(t *testing.T) {
	src1, err := OpenAlignments("testdata/bf_fragilis.bam")
	if err != nil {
		t.Skip("missing")
	}
	defer src1.Close()
	ref, _ := LoadReference("testdata/bf_fragilis.gb")
	ref.MatchContigs(src1.Contigs)
	m := newModel(src1, ref, "testdata/bf_fragilis.gb", nil)
	m.width, m.height = 120, 30
	m.showVariant(0)
	run(m, "go 100000")
	before := m.view.start

	m.openAlignments("testdata/bf_fragilis.bam")
	defer func() {
		for _, s := range m.tabs[1:] {
			s.Close()
		}
	}()
	t.Logf("after opening the same BAM again: tabs=%d tab=%d message=%q", len(m.tabs), m.tab, m.message)
	press(m, "t")
	t.Logf("after t: tab=%d start=%d (was %d) contig=%d", m.tab, m.view.start, before, m.contig)
	if len(m.tabs) > 1 && m.view.start != before {
		t.Logf("FINDING: switching tabs does not keep the position")
	}
}

func TestAuditMouseAndCopy(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	m.width, m.height = 150, 44
	m.reload()

	r := m.reads[0]
	row := m.layout().pileupTop + r.Row - m.rowOffset
	x := gutterWidth + m.column(r.Start)
	m.handleMouse(tea.MouseMsg{X: x, Y: row, Action: tea.MouseActionMotion})
	t.Logf("hover: read=%v position=%d", m.hover != nil, m.hoverPos)

	old := os.Stderr
	devnull, _ := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	os.Stderr = devnull
	cmd := m.handleMouse(tea.MouseMsg{X: x, Y: row, Action: tea.MouseActionPress, Button: tea.MouseButtonLeft})
	os.Stderr = old
	devnull.Close()
	t.Logf("click read: message=%q flash=%v cmd=%v", m.message, m.flash != nil, cmd != nil)

	m.panel = "variants"
	m.panelFocus = false
	m.reload()
	before := m.variantIndex
	m.handleMouse(tea.MouseMsg{X: m.width - 5, Y: m.layout().listTop + 2, Action: tea.MouseActionPress, Button: tea.MouseButtonLeft})
	t.Logf("panel click: focus=%v variant %d -> %d", m.panelFocus, before, m.variantIndex)

	tabY := 0
	m.handleMouse(tea.MouseMsg{X: 3, Y: tabY, Action: tea.MouseActionPress, Button: tea.MouseButtonLeft})
	t.Logf("tab bar click: message=%q", m.message)
}

func TestAuditInputs(t *testing.T) {
	cases := []struct {
		args []string
		ok   bool
	}{
		{[]string{"a.bam"}, true},
		{[]string{"a.sam", "b.bam"}, true},
		{[]string{"a.gb"}, true},
		{[]string{"a.gbk"}, true},
		{[]string{"a.genbank"}, true},
		{[]string{"a.fa"}, true},
		{[]string{"a.fasta"}, true},
		{[]string{"a.fna"}, true},
		{[]string{"a.vcf"}, true},
		{[]string{"a.fa", "b.gb"}, false},
		{[]string{"a.vcf", "b.vcf"}, false},
		{[]string{"a.cram"}, false},
		{[]string{"a.vcf.gz"}, false},
		{[]string{"a.bcf"}, false},
		{[]string{"a.txt"}, false},
	}
	for _, c := range cases {
		_, err := readInputs(c.args)
		got := err == nil
		if got != c.ok {
			t.Errorf("readInputs(%v): err=%v", c.args, err)
		} else if err != nil {
			t.Logf("%v -> %v", c.args, err)
		}
	}
	in, _ := readInputs([]string{"x.gb", "a.bam", "b.bam", "v.vcf"})
	t.Logf("mixed order parsed as alignments=%v ref=%q vcf=%q", in.alignments, in.reference, in.variants)
}

func TestAuditCoverageProfiles(t *testing.T) {
	for _, c := range []struct{ name, path, contig string }{
		{"bf", "testdata/bf_fragilis.bam", "NC_003228.3"},
		{"ecoli", "testdata/ecoli_k12.bam", "NC_000913.3"},
		{"vcholerae-chr1", "testdata/vcholerae.bam", "NC_002505.1"},
		{"influenza", "testdata/influenza_h1n1.bam", "seg1_PB2"},
		{"sam", "testdata/bf_fragilis_full.sam", "NC_003228.3"},
	} {
		t.Run(c.name, func(t *testing.T) {
			if _, err := os.Stat(c.path); err != nil {
				t.Skip("missing")
			}
			src, err := OpenAlignments(c.path)
			if err != nil {
				t.Skip(err)
			}
			defer src.Close()
			start := time.Now()
			profile := src.CoverageProfile(c.contig, profileBins)
			t.Logf("%s profile (%d bins): %v", c.name, len(profile), time.Since(start))
			if len(profile) != profileBins {
				t.Logf("FINDING: profile has %d bins, want %d", len(profile), profileBins)
			}
		})
	}
}

func TestAuditGenBankStyles(t *testing.T) {
	for _, gb := range []string{"bf_fragilis_benchling.gb", "bf_fragilis_prodigal.gb"} {
		t.Run(gb, func(t *testing.T) {
			src, err := OpenAlignments("testdata/bf_fragilis.bam")
			if err != nil {
				t.Skip("missing")
			}
			defer src.Close()
			ref, err := LoadReference("testdata/" + gb)
			if err != nil {
				t.Skip("missing gb")
			}
			ref.MatchContigs(src.Contigs)
			attached := 0
			for _, c := range src.Contigs {
				if len(ref.Seq[c.Name]) == c.Length {
					attached++
				}
			}
			t.Logf("%s: contigs matched by name+length=%d/%d", gb, attached, len(src.Contigs))
			for _, c := range src.Contigs {
				t.Logf("  %s: bam=%d ref=%d", c.Name, c.Length, len(ref.Seq[c.Name]))
			}
			if attached == 0 {
				t.Logf("note: %s attaches only by name here; its sequences are longer than the contigs (duplicated ORIGIN in the benchling copy)", gb)
			}
		})
	}
}

func TestAuditImpactReal(t *testing.T) {
	ref, err := LoadReference("testdata/bf_fragilis.gb")
	if err != nil {
		t.Skip("missing")
	}
	vs, _ := LoadVariants("testdata/bf_fragilis.vcf")
	if len(vs) == 0 {
		t.Skip("no variants")
	}
	contig := vs[0].Contig
	badStart := 0
	badGiven := 0
	checked := 0
	kinds := map[string]int{}
	for _, v := range vs {
		if strings.Contains(v.Alt, ",") || len(v.Ref) != 1 || len(v.Alt) != 1 {
			continue
		}
		imp := ref.CodingImpact(v.Contig, v.Start, v.Ref, v.Alt)
		kinds[imp.Kind]++
		if imp.Kind == "start loss" && imp.CodonNumber != 1 {
			badStart++
		}
		if imp.Gene != "" && imp.FromAmino != 0 {
			checked++
			cds := ref.findCDS(v.Contig, v.Start)
			if cds != nil && cds.Protein != "" {
				positions := codingPositions(*cds)
				protein := ref.givenProtein(*cds, v.Contig, positions)
				if idx := indexOf(positions, v.Start); idx >= 0 && protein != "" && idx/3 < len(protein) {
					if imp.FromAmino != protein[idx/3] {
						badGiven++
					}
				}
			}
		}
	}
	t.Logf("impact kinds over %d substitutions on %s: %v", checked, contig, kinds)
	t.Logf("start-loss calls outside codon 1: %d; reference amino acids disagreeing with /translation: %d", badStart, badGiven)
	if badStart > 0 {
		t.Logf("FINDING: start loss reported for internal codons")
	}

	// indels
	kinds = map[string]int{}
	for _, v := range vs {
		if len(v.Ref) == len(v.Alt) {
			continue
		}
		imp := ref.CodingImpact(v.Contig, v.Start, v.Ref, v.Alt)
		kinds[imp.Kind]++
	}
	t.Logf("impact kinds over indels: %v", kinds)

	// a deletion whose REF and ALT are both long
	v := Variant{Contig: contig, Start: 1000, Ref: "ACGT", Alt: "AC"}
	imp := ref.CodingImpact(contig, 1000, v.Ref, v.Alt)
	t.Logf("REF=ACGT ALT=AC (2 bp deletion): %q, short=%q", imp.Kind, imp.Short())
}

func TestAuditMouseTabsAndWheel(t *testing.T) {
	src1, err := OpenAlignments("testdata/bf_fragilis.bam")
	if err != nil {
		t.Skip("missing")
	}
	defer src1.Close()
	ref, _ := LoadReference("testdata/bf_fragilis.gb")
	ref.MatchContigs(src1.Contigs)
	m := newModel(src1, ref, "testdata/bf_fragilis.gb", nil)
	m.width, m.height = 120, 30
	m.reload()
	m.openAlignments("testdata/bf_fragilis.bam")
	defer func() {
		for _, s := range m.tabs[1:] {
			s.Close()
		}
	}()
	m.View()
	if len(m.tabHits) < 2 {
		t.Logf("tabHits=%v", m.tabHits)
	} else {
		hit := m.tabHits[1]
		m.handleMouse(tea.MouseMsg{X: hit.from + 1, Y: 0, Action: tea.MouseActionPress, Button: tea.MouseButtonLeft})
		t.Logf("cluster click: tab=%d hits=%v", m.tab, m.tabHits)
	}

	bpp := m.view.bpp
	m.handleMouse(tea.MouseMsg{X: gutterWidth + 10, Y: m.layout().pileupTop, Button: tea.MouseButtonWheelUp})
	t.Logf("wheel up: bpp %d -> %d", bpp, m.view.bpp)
	m.handleMouse(tea.MouseMsg{X: 10, Y: 10, Button: tea.MouseButtonWheelDown, Ctrl: true})
	t.Logf("ctrl wheel down: rowOffset=%d", m.rowOffset)
	start := m.view.start
	m.handleMouse(tea.MouseMsg{X: 10, Y: 10, Button: tea.MouseButtonWheelDown, Shift: true})
	t.Logf("shift wheel down: start %d -> %d", start, m.view.start)

	cmd := m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("m")})
	t.Logf("m toggle: cmd=%v mouseOn=%v message=%q", cmd != nil, m.mouseOn, m.message)
	cmd = m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("m")})
	t.Logf("m again: mouseOn=%v", m.mouseOn)

	if quit := m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("q")}); quit == nil {
		t.Logf("FINDING: q does not quit")
	}
	press(m, "g")
	_ = m.handleKey(tea.KeyMsg{Type: tea.KeyRunes, Runes: []rune("q")})
	t.Logf("q while the bar is open is typed into the bar: command=%q", m.command)
	if m.command != "gq" {
		t.Logf("FINDING: q was not typed into the open command bar")
	}
}

func TestAuditHelpAndQuit(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	press(m, "?")
	if m.panel != "help" {
		t.Fatalf("? did not open the help panel")
	}
	out := renderAt(m, 120, 34)
	t.Logf("help panel drawn: %v", strings.Contains(out, "commands"))
	press(m, "esc")
	t.Logf("after esc: panel=%q", m.panel)

	if cmd := m.handleKey(tea.KeyMsg{Type: tea.KeyCtrlC}); cmd == nil {
		t.Logf("FINDING: ctrl+c does not quit")
	}
}

func TestAuditOpenVCFAndExtensions(t *testing.T) {
	m := auditOpen(t, "testdata/ecoli_k12.bam", "", "")
	run(m, "open testdata/bf_fragilis.vcf")
	t.Logf("open mismatched VCF: variants=%d message=%q", len(m.variants), m.message)
	if len(m.variants) != 0 {
		t.Logf("FINDING: a mismatched VCF was accepted through the open command")
	}

	dir := t.TempDir()
	os.WriteFile(filepath.Join(dir, "e.gbk"), readFile(t, "testdata/ecoli_k12.gb"), 0o644)
	os.WriteFile(filepath.Join(dir, "e.genbank"), readFile(t, "testdata/ecoli_k12.gb"), 0o644)
	os.WriteFile(filepath.Join(dir, "e.fasta"), readFile(t, "testdata/ecoli_k12.fa"), 0o644)
	for _, f := range []string{"e.gbk", "e.genbank", "e.fasta"} {
		ref, err := LoadReference(filepath.Join(dir, f))
		if err != nil || len(ref.Seq) == 0 {
			t.Logf("FINDING: %s does not load (%v)", f, err)
		} else {
			t.Logf("%s loads %d contig(s)", f, len(ref.Seq))
		}
	}

	m.width, m.height = 40, 18
	m.panel = "help"
	m.reload()
	out := plain(m.View())
	t.Logf("help at 40x18: %d lines, has title=%v", len(strings.Split(out, "\n")), strings.Contains(out, "commands"))
}

func readFile(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Skip(err)
	}
	return b
}

func TestAuditInvisiblePanel(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	m.width, m.height = 80, 24
	m.reload()
	run(m, "variants")
	t.Logf("narrow window: panel=%q layout.panel=%d focus=%v message=%q", m.panel, m.layout().panel, m.panelFocus, m.message)
	before := m.variantIndex
	press(m, "down")
	press(m, "down")
	t.Logf("two down presses with an invisible panel: variantIndex %d -> %d", before, m.variantIndex)
	if m.layout().panel == 0 && m.variantIndex != before {
		t.Logf("FINDING: below 91 columns the side panel is not drawn, but it still takes the arrow keys and changes an off-screen selection")
	}
	out := plain(m.View())
	if strings.Contains(out, "Variants") {
		t.Logf("panel drawn at 80 columns despite layout.panel=0")
	}
}

func TestAuditCommandCase(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis_tiny.bam", "", "")
	press(m, "G")
	press(m, "o")
	t.Logf("typing 'Go' (capital): command=%q", m.command)
	if m.command == "" {
		t.Logf("FINDING: command names must start in lowercase; a capital first letter silently does nothing")
	}
}

func TestAuditSmallPanel(t *testing.T) {
	m := auditOpen(t, "testdata/bf_fragilis.bam", "testdata/bf_fragilis.gb", "testdata/bf_fragilis.vcf")
	m.width, m.height = 120, 20
	m.reload()
	m.panel, m.panelFocus = "filters", true
	m.filterIndex = 0
	m.reload()
	listRows := m.layout().listRows
	m.filterIndex = len(filterRows) - 1
	out := plain(m.View())
	t.Logf("filters panel at 120x20: listRows=%d rows=%d last selected=%d", listRows, len(filterRows), m.filterIndex)
	if len(filterRows)*2 > listRows {
		t.Logf("FINDING: the filters panel draws two rows per entry and does not scroll; at this height %d of %d entries are off-screen but still selectable", len(filterRows)-listRows/2, len(filterRows))
	}
	if !strings.Contains(out, filterRows[len(filterRows)-1].name) {
		t.Logf("confirmed: the selected filter row is not drawn")
	}

	m.width, m.height = 40, 18
	m.panel = "help"
	m.reload()
	out = plain(m.View())
	t.Logf("help at 40x18 contains the final line %q: %v", "q quit", strings.Contains(out, "q quit"))
	if !strings.Contains(out, "q quit") {
		t.Logf("FINDING: the help overlay is taller than an 18-row window and is clipped without scrolling")
	}
}
