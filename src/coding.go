package main

import (
	"fmt"
	"sort"
	"strings"
)

// Impact describes what a variant does to a coding sequence.
type Impact struct {
	Gene        string
	Kind        string // synonymous, missense, nonsense, frameshift, noncoding, ...
	GenePos     int    // position within the coding sequence, counting from 1
	CodonNumber int
	FromAmino   byte
	ToAmino     byte
}

// Short is the one-letter summary asked for at a glance: S for synonymous,
// N for a change in the protein.
func (i Impact) Short() string {
	switch {
	case i.Kind == "synonymous":
		return "S"
	case strings.Contains(i.Kind, "insertion"), strings.Contains(i.Kind, "deletion"):
		return "I"
	case i.Kind != "" && i.Kind != "noncoding":
		return "N"
	}
	return "·"
}

func (i Impact) String() string {
	switch {
	case i.Gene == "":
		return kindOrNoncoding(i.Kind)
	case i.ToAmino == 0 || i.FromAmino == i.ToAmino: // no protein change worth spelling out
		return fmt.Sprintf("%s:%d %s", i.Gene, i.GenePos, kindOrNoncoding(i.Kind))
	}
	return fmt.Sprintf("%s:%d protein: %c%d%c %s", i.Gene, i.GenePos, i.FromAmino, i.CodonNumber, i.ToAmino, i.Kind)
}

func kindOrNoncoding(kind string) string {
	if kind == "" {
		return "noncoding"
	}
	return kind
}

// CodingImpact works out the effect of replacing refAllele with alt at a
// zero-based reference position.
func (r *Reference) CodingImpact(contig string, pos int, refAllele, alt string) Impact {
	if r == nil {
		return Impact{}
	}
	shift := len(alt) - len(refAllele)
	cds := r.findCDS(contig, pos)
	if cds == nil {
		return Impact{Kind: outsideGene(shift)}
	}
	coding := r.positionsOf(cds)
	index := indexOf(coding, pos)
	if index < 0 {
		return Impact{Gene: geneName(*cds), Kind: outsideGene(shift)}
	}
	impact := Impact{Gene: geneName(*cds), GenePos: index + 1, CodonNumber: index/3 + 1}

	switch {
	case shift > 0 && shift%3 != 0:
		impact.Kind = "frameshift insertion"
	case shift < 0 && shift%3 != 0:
		impact.Kind = "frameshift deletion"
	case shift > 0:
		impact.Kind = "in-frame insertion"
	case shift < 0:
		impact.Kind = "in-frame deletion"
	}
	if shift != 0 {
		return impact
	}
	if len(refAllele) != 1 || len(alt) != 1 {
		impact.Kind = "complex substitution"
		return impact
	}

	first := index - index%3
	if first+2 >= len(coding) {
		impact.Kind = "noncoding"
		return impact
	}
	from := codonAt(r, *cds, contig, coding[first:first+3])
	to := append([]byte(nil), from...)
	newBase := alt[0]
	if cds.Reverse {
		newBase = complement(newBase)
	}
	to[index%3] = newBase

	fromAA := aminoAt(cds.Table, from, impact.CodonNumber)
	if protein := r.givenProtein(*cds, contig, coding); protein != "" && index/3 < len(protein) {
		fromAA = protein[index/3]
	}
	toAA := aminoAt(cds.Table, to, impact.CodonNumber)
	impact.FromAmino, impact.ToAmino = fromAA, toAA
	switch {
	case impact.CodonNumber == 1 && isStart(cds.Table, string(from)) && !isStart(cds.Table, string(to)):
		impact.Kind = "start loss"
	case fromAA == toAA:
		impact.Kind = "synonymous"
	case toAA == '*':
		impact.Kind = "nonsense"
	case fromAA == '*':
		impact.Kind = "stop loss"
	default:
		impact.Kind = "missense"
	}
	return impact
}

func (r *Reference) findCDS(contig string, pos int) *Feature {
	feats := r.featuresIn(contig, pos, pos+1)
	for i, f := range feats {
		if f.Kind != "CDS" {
			continue
		}
		for _, rng := range f.Ranges {
			if pos >= rng.Start && pos < rng.End {
				return &feats[i]
			}
		}
	}
	return nil
}

// outsideGene names what an indel does where there is no coding sequence.
func outsideGene(shift int) string {
	switch {
	case shift > 0:
		return "insertion"
	case shift < 0:
		return "deletion"
	}
	return "noncoding"
}

// positionsOf lists a CDS's reference positions in translation order and
// remembers them, because a VCF asks about the same gene thousands of times.
func (r *Reference) positionsOf(f *Feature) []int {
	if r.coding == nil {
		r.coding = map[*Feature][]int{}
	}
	if positions, ok := r.coding[f]; ok {
		return positions
	}
	positions := codingPositions(*f)
	r.coding[f] = positions
	return positions
}

// codingPositions lists every reference position of a CDS in translation order.
func codingPositions(f Feature) []int {
	var out []int
	for _, rng := range f.Ranges {
		for p := rng.Start; p < rng.End; p++ {
			out = append(out, p)
		}
	}
	if f.Reverse {
		for i, j := 0, len(out)-1; i < j; i, j = i+1, j-1 {
			out[i], out[j] = out[j], out[i]
		}
	}
	return out
}

func indexOf(positions []int, pos int) int {
	n := len(positions)
	if n == 0 {
		return -1
	}
	ascending := positions[0] <= positions[n-1]
	i := sort.Search(n, func(i int) bool {
		if ascending {
			return positions[i] >= pos
		}
		return positions[i] <= pos
	})
	if i < n && positions[i] == pos {
		return i
	}
	return -1
}

func geneName(f Feature) string {
	if f.Name != "" {
		return f.Name
	}
	if f.Product != "" {
		return f.Product
	}
	return f.Kind
}

func complement(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	}
	return b
}

// The standard genetic code, indexed with bases ordered T, C, A, G.
const codonTable = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

// builtInTable reports whether a genetic code is built in. Tables 1 and 11
// differ only in which codons may start a protein, so one codon table serves
// both.
func builtInTable(n int) bool { return n == 0 || n == 1 || n == 11 }

// isStart reports whether a codon can begin a protein under a genetic code.
// Any unbuilt table is read as the standard code.
func isStart(table int, codon string) bool {
	switch codon {
	case "ATG":
		return true
	case "GTG", "TTG", "CTG", "ATT", "ATC", "ATA":
		return table == 11
	}
	return false
}

// aminoAt is the amino acid a codon contributes at one place in a protein: a
// start codon in the first place is the methionine that begins it.
func aminoAt(table int, codon []byte, number int) byte {
	if number == 1 && isStart(table, string(codon)) {
		return 'M'
	}
	return translate(codon)
}

// codonAt reads the codon beginning at a list of reference positions.
func codonAt(r *Reference, f Feature, contig string, positions []int) []byte {
	codon := make([]byte, 3)
	for i := range codon {
		codon[i] = r.Base(contig, positions[i])
		if f.Reverse {
			codon[i] = complement(codon[i])
		}
	}
	return codon
}

// givenProtein is the translation from the file, when it covers the feature:
// a residue for every codon, or one for every codon but a closing stop, which
// GenBank leaves out of a CDS. A length that fits neither means the qualifier
// is stale, and a CDS is translated from its sequence instead. Features that
// are not CDS keep whatever translation they give.
func (r *Reference) givenProtein(f Feature, contig string, positions []int) string {
	codons := len(positions) / 3
	n := len(f.Protein)
	if n == 0 || codons == 0 || n > codons {
		return ""
	}
	if n != codons && f.Kind == "CDS" {
		if n != codons-1 || translate(codonAt(r, f, contig, positions[(codons-1)*3:])) != '*' {
			return ""
		}
	}
	return f.Protein
}

// assignTable gives a record that names no genetic code the bacterial one when
// one of its CDS begins with a bacterial alternative start, which the standard
// code cannot explain. Sense codons are the same either way, so only the first
// codon is read differently.
func assignTable(feats []Feature, seq []byte) {
	for _, f := range feats {
		if f.Kind != "CDS" || f.Table != 0 {
			continue
		}
		positions := codingPositions(f)
		if len(positions) < 3 {
			continue
		}
		codon := make([]byte, 3)
		for i := range codon {
			if positions[i] >= len(seq) {
				return
			}
			codon[i] = upper(seq[positions[i]])
			if f.Reverse {
				codon[i] = complement(codon[i])
			}
		}
		switch string(codon) {
		case "GTG", "TTG", "CTG", "ATT", "ATC", "ATA":
			for i := range feats {
				if feats[i].Kind == "CDS" && feats[i].Table == 0 {
					feats[i].Table = 11
				}
			}
			return
		}
	}
}

func translate(codon []byte) byte {
	index := 0
	for _, b := range codon {
		i := map[byte]int{'T': 0, 'C': 1, 'A': 2, 'G': 3}[b]
		if b != 'T' && b != 'C' && b != 'A' && b != 'G' {
			return 'X'
		}
		index = index*4 + i
	}
	return codonTable[index]
}
