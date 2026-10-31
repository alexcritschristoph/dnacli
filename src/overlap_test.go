package main

import (
	"strings"
	"testing"
)

// TestOverlappingCodingSequences pins down what overlapping CDSs look like:
// their bars go on separate rows, their proteins share one row where the
// longest wins each column, and an impact goes to the same longest CDS.
func TestOverlappingCodingSequences(t *testing.T) {
	model := ready(t, "overlapping")
	model.view = viewState{0, 1}
	model.reload()
	screen(t, model, 800, 44)

	contig := "OVERLAP"
	genes := map[string]Feature{}
	for _, f := range model.ref.Features[contig] {
		if f.Kind == "CDS" {
			genes[f.Name] = f
		}
	}
	if len(genes) != 3 {
		t.Fatalf("the example has %d CDS, want 3", len(genes))
	}
	if genes["alpha"].Reverse || genes["beta"].Reverse || !genes["gamma"].Reverse {
		t.Fatal("the strands of the example are wrong")
	}
	for name, f := range genes {
		if got := model.ref.givenProtein(f, contig, codingPositions(f)); got != f.Protein {
			t.Fatalf("%s does not use its own translation: %.30s vs %.30s", name, got, f.Protein)
		}
	}

	// alpha and beta overlap, as do beta and gamma; non-overlapping bars share a row.
	rows := map[string]int{}
	for _, name := range []string{"alpha", "beta", "gamma"} {
		rows[name] = geneRowNamed(model, name)
		if rows[name] < 0 {
			t.Fatalf("%s has no gene bar", name)
		}
	}
	if rows["alpha"] == rows["beta"] || rows["beta"] == rows["gamma"] {
		t.Fatalf("overlapping bars were not separated: %v", rows)
	}
	if rows["alpha"] != rows["gamma"] {
		t.Fatalf("non-overlapping bars did not share a row: %v", rows)
	}
	if !strings.Contains(string(frameRow(model, rows["gamma"])), "◀") {
		t.Fatal("the reverse gene has no left-pointing arrow")
	}

	// One protein row: features are painted shortest first, so the longest CDS
	// covering a column clears and overwrites it and wins there.
	columns := model.layout().columns
	proteinRow := rowWith(model, "protein")
	if proteinRow < 0 {
		t.Fatal("no protein row was drawn")
	}
	track := model.codingTrack(columns)
	want := make(map[int]rune)
	for _, f := range track {
		positions := codingPositions(f)
		protein := model.ref.givenProtein(f, contig, positions)
		for i := 0; i+2 < len(positions); i += 3 {
			col := model.column(positions[i+1])
			if col < 0 || col >= columns {
				continue
			}
			amino := byte('X')
			if i/3 < len(protein) {
				amino = protein[i/3]
			} else if f.Kind == "CDS" {
				amino = aminoAt(f.Table, codonAt(model.ref, f, contig, positions[i:i+3]), i/3+1)
			}
			for k := 0; k < 3; k++ {
				if c := model.column(positions[i+k]); c >= 0 && c < columns {
					want[c] = ' '
				}
			}
			want[col] = rune(amino)
		}
	}
	for col := 0; col < columns; col++ {
		expect := want[col]
		if expect == 0 {
			expect = ' '
		}
		if got := model.frame[proteinRow][gutterWidth+col].r; got != expect {
			t.Fatalf("column %d drew %q, want %q", col, got, expect)
		}
	}

	// The position belongs to two CDSs; the longer one owns it, and a tie keeps
	// the one the file lists first.
	cases := []struct {
		pos  int
		gene string
	}{
		{200, "alpha"}, // alpha alone
		{350, "beta"},  // alpha (300 bp) and beta (399 bp)
		{550, "beta"},  // beta and gamma, both 399 bp
		{800, "gamma"}, // gamma alone
	}
	for _, c := range cases {
		cds := model.ref.findCDS(contig, c.pos)
		if cds == nil || cds.Name != c.gene {
			t.Fatalf("position %d went to %v, want %s", c.pos, cds, c.gene)
		}
		base := string(model.ref.Base(contig, c.pos))
		alt := "A"
		if base == alt {
			alt = "C"
		}
		impact := model.ref.CodingImpact(contig, c.pos, base, alt)
		if impact.Gene != c.gene {
			t.Fatalf("the impact at %d came out %+v, want %s", c.pos, impact, c.gene)
		}
	}
}

// geneRowNamed finds the gene row that carries a feature's name.
func geneRowNamed(m *model, name string) int {
	label := gutterWidth - 2 - len("genes")
	for row := range m.frame {
		if m.frame[row][label].r != 'g' {
			continue
		}
		if strings.Contains(string(frameRow(m, row)), name) {
			return row
		}
	}
	return -1
}
