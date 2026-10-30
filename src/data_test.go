package main

import (
	"os"
	"path/filepath"
	"testing"
)

func TestGenBankKeepsTranslationQualifiers(t *testing.T) {
	const gb = `LOCUS       TEST                      18 bp    DNA     linear   BCT 01-JAN-2026
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="a"
                     /transl_table=11
                     /translation="M
                     K"
     misc_feature    10..18
                     /locus_tag="b"
                     /translation="WW"
     CDS             10..18
                     /locus_tag="c"
                     /transl_table=2
                     /translation="WW"
ORIGIN
        1 gtgaaataatggtggtaa
//`
	path := filepath.Join(t.TempDir(), "test.gb")
	if err := os.WriteFile(path, []byte(gb), 0o644); err != nil {
		t.Fatal(err)
	}
	ref, err := LoadReference(path)
	if err != nil {
		t.Fatal(err)
	}
	features := map[string]Feature{}
	for _, f := range ref.Features["TEST"] {
		features[f.Name] = f
	}
	if f := features["a"]; f.Protein != "MK" || f.Table != 11 {
		t.Fatalf("CDS a came out as %q, table %d", f.Protein, f.Table)
	}
	if f := features["b"]; f.Kind != "misc_feature" || f.Protein != "WW" {
		t.Fatalf("the misc_feature lost its translation: %q", f.Protein)
	}
	if f := features["c"]; f.Table != 2 {
		t.Fatalf("the unsupported table was not recorded: %d", f.Table)
	}
	if ref.UnsupportedTable != 2 {
		t.Fatalf("the reference reports table %d as unsupported, want 2", ref.UnsupportedTable)
	}
	if got := ref.givenProtein(features["a"], "TEST", codingPositions(features["a"])); got != "MK" {
		t.Fatalf("the CDS translation came back as %q", got)
	}
	if got := ref.givenProtein(features["b"], "TEST", codingPositions(features["b"])); got != "WW" {
		t.Fatalf("the misc_feature translation came back as %q", got)
	}
}

func TestBacterialStartCodonImpacts(t *testing.T) {
	ref := newReference()
	ref.Seq["X"] = []byte("GTGAAATAA")
	gene := Feature{Kind: "CDS", Name: "gtg", Table: 11, Protein: "MK", Ranges: []Range{{0, 9}}}
	ref.Features["X"] = []Feature{gene}
	ref.longest["X"] = 9

	cases := []struct {
		pos      int
		ref, alt string
		kind     string
	}{
		{0, "G", "T", "synonymous"}, // GTG -> TTG, both first codons
		{0, "G", "A", "synonymous"}, // GTG -> ATG, both first codons
		{2, "G", "A", "start loss"}, // GTG -> GTA, no longer a start
		{3, "A", "G", "missense"},   // AAA -> GAA, K -> E
		{5, "A", "G", "synonymous"}, // AAA -> AAG
		{3, "A", "T", "nonsense"},   // AAA -> TAA
	}
	for _, c := range cases {
		impact := ref.CodingImpact("X", c.pos, c.ref, c.alt)
		if impact.Kind != c.kind {
			t.Fatalf("%s>%s at %d: %s, want %s", c.ref, c.alt, c.pos, impact.Kind, c.kind)
		}
		if c.pos < 3 && impact.FromAmino != 'M' {
			t.Fatalf("the first residue came out %c, want M", impact.FromAmino)
		}
	}
	if got := ref.givenProtein(gene, "X", codingPositions(gene)); got != "MK" {
		t.Fatalf("the file's translation came back as %q", got)
	}

	bare := newReference()
	bare.Seq["X"] = ref.Seq["X"]
	bare.Features["X"] = []Feature{{Kind: "CDS", Name: "bare", Table: 11, Ranges: []Range{{0, 9}}}}
	bare.longest["X"] = 9
	if got := bare.givenProtein(bare.Features["X"][0], "X", codingPositions(gene)); got != "" {
		t.Fatalf("a CDS without a translation still gave %q", got)
	}
	if impact := bare.CodingImpact("X", 0, "G", "T"); impact.Kind != "synonymous" || impact.FromAmino != 'M' {
		t.Fatalf("without a translation, the first codon came out %c %s", impact.FromAmino, impact.Kind)
	}
}

func TestAlternativeStartImpliesBacterialTable(t *testing.T) {
	const gb = `LOCUS       TEST                       9 bp    DNA     linear   BCT 01-JAN-2026
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="a"
ORIGIN
        1 gtgaaataa
//`
	path := filepath.Join(t.TempDir(), "test.gb")
	if err := os.WriteFile(path, []byte(gb), 0o644); err != nil {
		t.Fatal(err)
	}
	ref, err := LoadReference(path)
	if err != nil {
		t.Fatal(err)
	}
	if f := ref.Features["TEST"][0]; f.Table != 11 {
		t.Fatalf("a GTG-start CDS with no table came out as table %d", f.Table)
	}
}

func TestBacterialExampleCarriesItsOwnTranslations(t *testing.T) {
	ref, err := LoadReference("testdata/bf_fragilis.gb")
	if err != nil {
		t.Fatal(err)
	}
	for _, f := range ref.Features["NC_003228"] {
		if f.Kind != "CDS" || f.Table != 11 || f.Protein == "" {
			continue
		}
		positions := codingPositions(f)
		if string(codonAt(ref, f, "NC_003228", positions)) != "GTG" {
			continue
		}
		if f.Protein[0] != 'M' {
			t.Fatalf("%s begins %c, want the initiator M", f.Name, f.Protein[0])
		}
		return
	}
	t.Fatal("no GTG-start CDS in the bacterial example")
}

func TestGenBankFeaturesAndSequence(t *testing.T) {
	ref, err := LoadReference("testdata/viral_variants.gb")
	if err != nil {
		t.Fatal(err)
	}
	if got := len(ref.Seq["VIRUS1"]); got != 8000 {
		t.Fatalf("sequence length = %d, want 8000", got)
	}
	var spike Feature
	for _, f := range ref.Features["VIRUS1"] {
		if f.Kind == "CDS" && f.Name == "spike" {
			spike = f
		}
	}
	if spike.Start() != 3199 || spike.End() != 5000 {
		t.Fatalf("spike CDS = %d..%d, want 3199..5000", spike.Start(), spike.End())
	}
	for _, f := range ref.Features["VIRUS1"] {
		if f.Name == "nucleo" && f.Kind == "CDS" && !f.Reverse {
			t.Fatal("nucleo CDS should be on the reverse strand")
		}
	}
}

func TestVariantImpacts(t *testing.T) {
	ref, _ := LoadReference("testdata/viral_variants.gb")
	variants, err := LoadVariants("testdata/viral_variants.vcf")
	if err != nil {
		t.Fatal(err)
	}
	if len(variants) != 9 {
		t.Fatalf("loaded %d variants, want 9", len(variants))
	}
	for _, v := range variants {
		impact := ref.CodingImpact(v.Contig, v.Start, v.Ref, v.Alt)
		if impact.Kind == "" {
			t.Fatalf("no impact worked out for %s:%d", v.Contig, v.Start+1)
		}
		t.Logf("%d %s %s", v.Start+1, v.Label(), impact)
	}
}

func TestReadsAreLaidOutOnTheReference(t *testing.T) {
	ref, _ := LoadReference("testdata/viral_variants.gb")
	src, err := OpenAlignments("testdata/viral_variants.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	if len(src.Contigs) != 1 || src.Contigs[0].Name != "VIRUS1" {
		t.Fatalf("contigs = %v", src.Contigs)
	}
	reads, _, _ := src.Reads("VIRUS1", 1150, 1250, defaultFilter(), ref, true)
	if len(reads) == 0 {
		t.Fatal("no reads near the planted variant at 1201")
	}
	alt := 0
	for _, r := range reads {
		if r.BaseAt(1200) == 'C' {
			alt++
		}
		if r.End-r.Start != len(r.Bases) {
			t.Fatalf("read %s spans %d but has %d bases", r.Name, r.End-r.Start, len(r.Bases))
		}
	}
	if alt == 0 {
		t.Fatal("expected reads carrying the A>C variant at position 1201")
	}
}

func TestFilters(t *testing.T) {
	src, err := OpenAlignments("testdata/supertest.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	contig := src.Contigs[0].Name
	all, _, _ := src.Reads(contig, 0, src.Contigs[0].Length, defaultFilter(), nil, true)
	strict := defaultFilter()
	strict.MinMapQ = 30
	kept, _, _ := src.Reads(contig, 0, src.Contigs[0].Length, strict, nil, true)
	if len(kept) >= len(all) || len(kept) == 0 {
		t.Fatalf("mapping quality filter kept %d of %d reads", len(kept), len(all))
	}
	for _, r := range kept {
		if r.MapQ < 30 {
			t.Fatalf("read %s has mapping quality %d", r.Name, r.MapQ)
		}
	}
}

func TestCoverageProfile(t *testing.T) {
	src, err := OpenAlignments("testdata/multi_contig_uneven.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	empty := 0
	for _, c := range src.Contigs {
		profile := src.CoverageProfile(c.Name, 24)
		total := 0.0
		for _, v := range profile {
			total += v
		}
		if total == 0 {
			empty++
		}
	}
	if empty != 1 {
		t.Fatalf("%d contigs had no coverage, want exactly 1", empty)
	}
}

func TestFeatureSearchFindsEveryOverlap(t *testing.T) {
	ref, err := LoadReference("testdata/viral_variants.gb")
	if err != nil {
		t.Fatal(err)
	}
	for start := 0; start < 8000; start += 137 {
		for _, width := range []int{1, 100, 3000} {
			found := map[string]bool{}
			for _, f := range ref.featuresIn("VIRUS1", start, start+width) {
				found[f.Kind+f.Name] = true
			}
			for _, f := range ref.Features["VIRUS1"] {
				if f.End() > start && f.Start() < start+width && !found[f.Kind+f.Name] {
					t.Fatalf("the search of %d..%d missed %s %s", start, start+width, f.Kind, f.Name)
				}
			}
		}
	}
}
