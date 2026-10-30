package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// revComp is the reverse complement of a short test sequence.
func revComp(s string) string {
	b := []byte(s)
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = complement(b[j]), complement(b[i])
	}
	if len(b)%2 == 1 {
		b[len(b)/2] = complement(b[len(b)/2])
	}
	return string(b)
}

// testReference writes a GenBank file, loads it, and returns the reference.
func testReference(t *testing.T, record string) *Reference {
	t.Helper()
	path := filepath.Join(t.TempDir(), "test.gb")
	if err := os.WriteFile(path, []byte(record), 0o644); err != nil {
		t.Fatal(err)
	}
	ref, err := LoadReference(path)
	if err != nil {
		t.Fatal(err)
	}
	return ref
}

// genBankRecord wraps a sequence and its features in a one-contig record.
func genBankRecord(name, sequence, features string) string {
	return fmt.Sprintf(`LOCUS       %-16s %6d bp    DNA     linear   SYN 01-JAN-2026
FEATURES             Location/Qualifiers
%sORIGIN
        1 %s
//`, name, len(sequence), features, sequence)
}

// codingProtein translates a feature from its reference sequence, in coding order.
func codingProtein(ref *Reference, f Feature, contig string) string {
	var out []byte
	positions := codingPositions(f)
	for i := 0; i+2 < len(positions); i += 3 {
		out = append(out, aminoAt(f.Table, codonAt(ref, f, contig, positions[i:i+3]), i/3+1))
	}
	return string(out)
}

func TestReverseLocationsAndTranslations(t *testing.T) {
	coding := "ATGAAACCCGGGTTTGAAGATTCCAAGTAA"
	sequence := coding + "AAAAAAAAAA" + revComp(coding) + "AAAAAAAAAA"
	ref := testReference(t, genBankRecord("REV", sequence, `     CDS             1..30
                     /locus_tag="plus"
     CDS             complement(41..70)
                     /locus_tag="minus"
`))
	if len(ref.Features["REV"]) != 2 {
		t.Fatalf("loaded %d features", len(ref.Features["REV"]))
	}
	plus, minus := ref.Features["REV"][0], ref.Features["REV"][1]
	if plus.Reverse || !minus.Reverse {
		t.Fatalf("strands came out plus %v, minus %v", plus.Reverse, minus.Reverse)
	}
	if minus.Start() != 40 || minus.End() != 70 || len(minus.Ranges) != 1 {
		t.Fatalf("the reverse location came out %v (%d..%d)", minus.Ranges, minus.Start(), minus.End())
	}
	if got := codingProtein(ref, plus, "REV"); got != "MKPGFEDSK*" {
		t.Fatalf("the forward gene translates as %q", got)
	}
	if got := codingProtein(ref, minus, "REV"); got != "MKPGFEDSK*" {
		t.Fatalf("the reverse gene translates as %q", got)
	}
	positions := codingPositions(minus)
	if positions[0] != 69 || positions[len(positions)-1] != 40 {
		t.Fatalf("the reverse coding positions run from %d to %d", positions[0], positions[len(positions)-1])
	}
	for i := 0; i+2 < len(positions); i += 3 {
		forward := codonAt(ref, plus, "REV", []int{i, i + 1, i + 2})
		reverse := codonAt(ref, minus, "REV", positions[i:i+3])
		if string(forward) != string(reverse) {
			t.Fatalf("codon %d is %s forward and %s reverse", i/3+1, forward, reverse)
		}
	}
}

func TestReverseImpactsMirrorForward(t *testing.T) {
	coding := "ATGAAACCCGGGTTTGAAGATTCCAAGTAA"
	sequence := coding + "AAAAAAAAAA" + revComp(coding) + "AAAAAAAAAA"
	ref := testReference(t, genBankRecord("REV", sequence, `     CDS             1..30
                     /locus_tag="plus"
     CDS             complement(41..70)
                     /locus_tag="minus"
`))
	minus := ref.Features["REV"][1]
	positions := codingPositions(minus)
	for i := 0; i < len(coding); i++ {
		for _, alt := range "ACGT" {
			if byte(alt) == coding[i] {
				continue
			}
			forward := ref.CodingImpact("REV", i, string(coding[i]), string(alt))
			reverse := ref.CodingImpact("REV", positions[i], string(complement(coding[i])), string(complement(byte(alt))))
			if forward.Kind != reverse.Kind || forward.GenePos != reverse.GenePos ||
				forward.CodonNumber != reverse.CodonNumber || forward.FromAmino != reverse.FromAmino ||
				forward.ToAmino != reverse.ToAmino {
				t.Fatalf("coding base %d %c>%c: forward %+v, reverse %+v",
					i, coding[i], alt, forward, reverse)
			}
		}
	}
}

func TestReverseMultiExonCodingAndImpacts(t *testing.T) {
	coding := "ATGAAACCCGGGTTTGAAGATTCCAAGTAA"
	split := 13
	first, second := revComp(coding[split:]), revComp(coding[:split])
	sequence := coding + "AAAAAAAAAA" + first + "AAAAAAAAAA" + second + "AAAAAAAAAA"
	ref := testReference(t, genBankRecord("REVJ", sequence, `     CDS             1..30
                     /locus_tag="whole"
     CDS             complement(join(41..57,68..80))
                     /locus_tag="spliced"
`))
	spliced := ref.Features["REVJ"][1]
	if !spliced.Reverse || len(spliced.Ranges) != 2 {
		t.Fatalf("the spliced location came out %v, reverse %v", spliced.Ranges, spliced.Reverse)
	}
	if got := codingProtein(ref, spliced, "REVJ"); got != "MKPGFEDSK*" {
		t.Fatalf("the spliced reverse gene translates as %q", got)
	}
	positions := codingPositions(spliced)
	if positions[split-1]-positions[split] < 2 {
		t.Fatalf("the exon junction does not jump: %d then %d", positions[split-1], positions[split])
	}
	// The codon straddling the junction reads one exon and then the other.
	codon := codonAt(ref, spliced, "REVJ", positions[12:15])
	if string(codon) != coding[12:15] {
		t.Fatalf("the junction codon is %s, want %s", codon, coding[12:15])
	}

	for i := 0; i < len(coding); i++ {
		for _, alt := range "ACGT" {
			if byte(alt) == coding[i] {
				continue
			}
			forward := ref.CodingImpact("REVJ", i, string(coding[i]), string(alt))
			reverse := ref.CodingImpact("REVJ", positions[i], string(complement(coding[i])), string(complement(byte(alt))))
			if forward.Kind != reverse.Kind || forward.GenePos != reverse.GenePos ||
				forward.CodonNumber != reverse.CodonNumber || forward.FromAmino != reverse.FromAmino ||
				forward.ToAmino != reverse.ToAmino {
				t.Fatalf("coding base %d: whole %+v, spliced %+v", i, forward, reverse)
			}
		}
	}
}

func TestReverseBacterialStartsAndTableInference(t *testing.T) {
	ref := testReference(t, genBankRecord("BACT", revComp("GTGAAATAA"), `     CDS             complement(1..9)
                     /locus_tag="b"
`))
	gene := ref.Features["BACT"][0]
	if gene.Table != 11 {
		t.Fatalf("a reverse GTG start gave table %d", gene.Table)
	}
	if got := codingProtein(ref, gene, "BACT"); got != "MK*" {
		t.Fatalf("the reverse GTG gene translates as %q", got)
	}
	positions := codingPositions(gene)
	start := ref.CodingImpact("BACT", positions[2], string(complement('G')), string(complement('A')))
	if start.Kind != "start loss" || start.FromAmino != 'M' || start.ToAmino != 'V' {
		t.Fatalf("GTG>GTA on the reverse strand gave %+v", start)
	}
	other := ref.CodingImpact("BACT", positions[0], string(complement('G')), string(complement('T')))
	if other.Kind != "synonymous" || other.FromAmino != 'M' {
		t.Fatalf("GTG>TTG on the reverse strand gave %+v", other)
	}

	standard := testReference(t, genBankRecord("STD", revComp("ATGAAATAA"), `     CDS             complement(1..9)
                     /locus_tag="s"
`))
	gene = standard.Features["STD"][0]
	positions = codingPositions(gene)
	loss := standard.CodingImpact("STD", positions[2], string(complement('G')), string(complement('A')))
	if loss.Kind != "start loss" || loss.FromAmino != 'M' {
		t.Fatalf("ATG>ATA on the reverse strand gave %+v", loss)
	}
}

func TestReverseTranslationQualifierIsUsed(t *testing.T) {
	coding := "ATGAAACCCGGGTTTGAAGATTCCAAGTAA"
	sequence := "AAAAAAAAAA" + revComp(coding)
	ref := testReference(t, genBankRecord("TRANS", sequence, `     CDS             complement(11..40)
                     /locus_tag="minus"
                     /translation="WWWWWWWWWW"
`))
	gene := ref.Features["TRANS"][0]
	if got := ref.givenProtein(gene, "TRANS", codingPositions(gene)); got != "WWWWWWWWWW" {
		t.Fatalf("the reverse translation came back as %q", got)
	}
	model := modelFor(t, ref, "TRANS")
	letters := aminoLetters(model, gene)
	if len(letters) != 10 || strings.Count(letters, "W") != 10 {
		t.Fatalf("the reverse protein row drew %q", letters)
	}

	// A translation that fits neither the codons nor all but the stop is stale
	// and the sequence is translated instead.
	gene.Protein = "WWWWW"
	if got := ref.givenProtein(gene, "TRANS", codingPositions(gene)); got != "" {
		t.Fatalf("a stale reverse translation was taken: %q", got)
	}
	// A CDS translation that omits the closing stop fits one residue fewer.
	gene.Protein = "WWWWWWWWW"
	if got := ref.givenProtein(gene, "TRANS", codingPositions(gene)); got != "WWWWWWWWW" {
		t.Fatalf("a reverse translation without its stop came back as %q", got)
	}
}

func TestReverseGenesInTheView(t *testing.T) {
	coding := "ATGAAACCCGGGTTTGAAGATTCCAAGTAA"
	sequence := coding + "AAAAAAAAAA" + revComp(coding) + "AAAAAAAAAA"
	ref := testReference(t, genBankRecord("REV", sequence, `     gene            1..30
                     /gene="plus"
     CDS             1..30
                     /gene="plus"
     gene            complement(41..70)
                     /gene="minus"
     CDS             complement(41..70)
                     /gene="minus"
                     /translation="MKPGFEDSK"
`))
	model := modelFor(t, ref, "REV")
	lines := screen(t, model, 120, 30)
	geneRow := rowWith(model, "genes")
	if geneRow < 0 {
		t.Fatal("no gene row was drawn")
	}
	text := string(frameRow(model, geneRow))
	if !strings.Contains(text, "plus") || !strings.Contains(text, "minus") {
		t.Fatalf("the gene row does not name both genes: %q", text)
	}
	if !strings.Contains(text, "▶") || !strings.Contains(text, "◀") {
		t.Fatalf("the gene row does not show both strands: %q", text)
	}
	if !strings.Contains(strings.Join(lines, "\n"), "minus") {
		t.Fatal("the reverse gene is missing from the screen")
	}

	minus := ref.Features["REV"][1]
	want := "MKPGFEDSK*"
	positions := codingPositions(minus)
	proteinRow := rowWith(model, "protein")
	if proteinRow < 0 {
		t.Fatal("no protein row was drawn")
	}
	for i := 0; i < len(want); i++ {
		col := model.column(positions[3*i+1])
		if col < 0 {
			continue
		}
		if got := model.frame[proteinRow][gutterWidth+col].r; got != rune(want[i]) {
			t.Fatalf("reverse codon %d drew %c at column %d, want %c", i+1, got, col, want[i])
		}
	}

	// The protein row still lines up when the view opens inside the gene.
	model.view = viewState{45, 1}
	model.reload()
	model.View()
	proteinRow = rowWith(model, "protein")
	if proteinRow < 0 {
		t.Fatal("no protein row while inside the reverse gene")
	}
	for i := 0; i < len(want); i++ {
		col := model.column(positions[3*i+1])
		if col < 0 {
			continue
		}
		if got := model.frame[proteinRow][gutterWidth+col].r; got != rune(want[i]) {
			t.Fatalf("inside the gene, codon %d drew %c, want %c", i+1, got, want[i])
		}
	}

	run(model, "gene minus")
	if model.view.start > 40 || model.view.start+model.span() < 70 {
		t.Fatalf("gene minus opened at %v", model.view)
	}
}

func TestReverseLocationParsing(t *testing.T) {
	cases := []struct {
		location string
		ranges   int
		reverse  bool
	}{
		{"complement(41..70)", 1, true},
		{"complement(<41..>70)", 1, true},
		{"complement(join(41..50,60..70))", 2, true},
		{"join(41..50,60..70)", 2, false},
		{"complement(complement(41..70))", 1, false},
	}
	for _, c := range cases {
		ranges, reverse := parseLocation(c.location)
		if len(ranges) != c.ranges || reverse != c.reverse {
			t.Fatalf("%s gave %v, reverse %v", c.location, ranges, reverse)
		}
	}
	ref := testReference(t, genBankRecord("EDGE", revComp("ATGAAATAA")+"AAAAAAAAAA", `     CDS             complement(1..9)
`))
	gene := ref.Features["EDGE"][0]
	if gene.Start() != 0 || gene.End() != 9 {
		t.Fatalf("an edge gene came out %d..%d", gene.Start(), gene.End())
	}
}

func TestReverseAnnotationsMatchTheirProteins(t *testing.T) {
	ref, err := LoadReference("testdata/bf_fragilis.gb")
	if err != nil {
		t.Fatal(err)
	}
	start := time.Now()
	checked, matched := 0, 0
	for contig, features := range ref.Features {
		for _, f := range features {
			if f.Kind != "CDS" || !f.Reverse || f.Protein == "" {
				continue
			}
			checked++
			got := codingProtein(ref, f, contig)
			if len(got) < len(f.Protein) {
				t.Fatalf("%s translates to %d residues, the file gives %d", f.Name, len(got), len(f.Protein))
			}
			if got[:len(f.Protein)] != f.Protein {
				t.Fatalf("%s:%s translates as %.40s, the file says %.40s", contig, f.Name, got, f.Protein)
			}
			matched++
		}
	}
	t.Logf("%d reverse CDS checked against their translations in %s", matched, time.Since(start).Round(time.Microsecond))
	if checked == 0 || matched != checked {
		t.Fatalf("%d of %d reverse CDS match", matched, checked)
	}
}

// modelFor puts a reference on an empty SAM with the same contig, so that the
// coding rows can be drawn without any reads.
func modelFor(t *testing.T, ref *Reference, contig string) *model {
	t.Helper()
	path := filepath.Join(t.TempDir(), "empty.sam")
	header := fmt.Sprintf("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:%s\tLN:%d\n", contig, len(ref.Seq[contig]))
	if err := os.WriteFile(path, []byte(header), 0o644); err != nil {
		t.Fatal(err)
	}
	src, err := OpenAlignments(path)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(src.Close)
	ref.MatchContigs(src.Contigs)
	m := newModel(src, ref, "", nil)
	m.width, m.height = 120, 30
	m.view = viewState{0, 1}
	m.reload()
	return m
}

// rowWith finds a track row by the label in its gutter.
func rowWith(m *model, label string) int {
	col := gutterWidth - 2 - len(label)
	for row := range m.frame {
		if col >= 0 && m.frame[row][col].r == rune(label[0]) {
			return row
		}
	}
	return -1
}

func frameRow(m *model, row int) []rune {
	var out []rune
	for _, c := range m.frame[row] {
		out = append(out, c.r)
	}
	return out
}

// aminoLetters draws one feature's protein row and returns the letters, left to right.
func aminoLetters(m *model, f Feature) string {
	g := newGrid(m.width, 2)
	m.drawAminoAcids(g, 0, m.layout(), f)
	var out []rune
	for _, c := range g[0] {
		if c.r != ' ' {
			out = append(out, c.r)
		}
	}
	return string(out)
}
