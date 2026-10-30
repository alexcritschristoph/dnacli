package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Range is a zero-based, half-open interval on a contig.
type Range struct{ Start, End int }

// Feature is one annotation from a GenBank file, such as a gene or a CDS.
type Feature struct {
	Kind    string
	Name    string
	Product string
	Protein string // the /translation qualifier, when the file carries one
	Table   int    // the /transl_table qualifier, 0 when the file names none
	Ranges  []Range
	Reverse bool
}

func (f Feature) Start() int { return f.Ranges[0].Start }
func (f Feature) End() int   { return f.Ranges[len(f.Ranges)-1].End }

// Reference holds the sequence and annotations of every contig we know about.
type Reference struct {
	Seq              map[string][]byte
	Features         map[string][]Feature // sorted by start position
	longest          map[string]int       // longest feature on a contig, so a search knows how far back to look
	coding           map[*Feature][]int   // each CDS's coding positions, worked out once
	UnsupportedTable int                  // the first genetic code named that is not built in
}

func newReference() *Reference {
	return &Reference{Seq: map[string][]byte{}, Features: map[string][]Feature{}, longest: map[string]int{}}
}

// featuresIn returns the features that may overlap a range. Looking only at
// those keeps a reference with thousands of annotations from being read end to
// end for every frame.
func (r *Reference) featuresIn(contig string, start, end int) []Feature {
	feats := r.Features[contig]
	first := sort.Search(len(feats), func(i int) bool { return feats[i].Start() > start-r.longest[contig] })
	last := sort.Search(len(feats), func(i int) bool { return feats[i].Start() >= end })
	return feats[first:last]
}

// Base returns the reference base at a zero-based position, or 0 if unknown.
func (r *Reference) Base(contig string, pos int) byte {
	if r == nil {
		return 0
	}
	s := r.Seq[contig]
	if pos < 0 || pos >= len(s) {
		return 0
	}
	return upper(s[pos])
}

// Span returns the reference bases over a range.
func (r *Reference) Span(contig string, from, to int) string {
	var out []byte
	for pos := from; pos < to; pos++ {
		if base := r.Base(contig, pos); base != 0 {
			out = append(out, base)
		}
	}
	return string(out)
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// LoadReference reads a GenBank or FASTA file, chosen by file extension.
func LoadReference(path string) (*Reference, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	lower := strings.ToLower(path)
	if strings.HasSuffix(lower, ".gb") || strings.HasSuffix(lower, ".gbk") || strings.HasSuffix(lower, ".genbank") {
		return readGenBank(f)
	}
	return readFasta(f)
}

func readFasta(f *os.File) (*Reference, error) {
	ref := newReference()
	name := ""
	var seq []byte
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	line := 0
	for sc.Scan() {
		line++
		text := sc.Text()
		if strings.HasPrefix(text, ">") {
			if name != "" {
				ref.Seq[name] = seq
			}
			fields := strings.Fields(text[1:])
			if len(fields) == 0 {
				return nil, fmt.Errorf("line %d: a FASTA record has no name after >", line)
			}
			name, seq = fields[0], nil
			continue
		}
		text = strings.TrimSpace(text)
		if strings.ContainsAny(text, " \t") { // aligned FASTA pads the bases
			text = strings.Join(strings.Fields(text), "")
		}
		seq = append(seq, text...)
	}
	if name != "" {
		ref.Seq[name] = seq
	}
	return ref, sc.Err()
}

// readGenBank handles multi-record files; each record contributes one contig.
func readGenBank(f *os.File) (*Reference, error) {
	ref := newReference()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)

	section := ""
	var names []string
	var seq []byte
	var feats []Feature
	var location string  // location text of the feature being read, possibly wrapped
	var qualifier string // the qualifier whose value may continue on the next lines

	finishFeature := func() {
		if location != "" {
			feats[len(feats)-1].Ranges, feats[len(feats)-1].Reverse = parseLocation(location)
			location = ""
		}
		if len(feats) > 0 {
			f := &feats[len(feats)-1]
			f.Name = strings.Trim(f.Name, "\" \t")
			f.Product = strings.Trim(f.Product, "\" \t")
			f.Protein = strings.Trim(f.Protein, "\" \t")
		}
		qualifier = ""
	}
	// A record is stored under its locus, accession and version names, so that
	// an alignment file naming the contig any of those ways still lines up.
	finishRecord := func() {
		finishFeature()
		kept, longest := feats[:0], 0
		for _, f := range feats {
			if len(f.Ranges) > 0 && f.Kind != "source" { // source spans the whole contig and says nothing
				kept, longest = append(kept, f), max(longest, f.End()-f.Start())
			}
		}
		feats = kept
		sort.Slice(feats, func(i, j int) bool { return feats[i].Start() < feats[j].Start() })
		assignTable(feats, seq)
		for _, name := range names {
			ref.Seq[name] = seq
			ref.Features[name], ref.longest[name] = feats, longest
		}
		names, section, seq, feats = nil, "", nil, nil
	}

	for sc.Scan() {
		line := sc.Text()
		switch {
		case strings.HasPrefix(line, "LOCUS"), strings.HasPrefix(line, "ACCESSION"), strings.HasPrefix(line, "VERSION"):
			if fields := strings.Fields(line); len(fields) > 1 && fields[1] != "." {
				names = append(names, fields[1])
			}
		case strings.HasPrefix(line, "FEATURES"):
			section = "features"
		case strings.HasPrefix(line, "ORIGIN"):
			finishFeature()
			section = "origin"
		case strings.HasPrefix(line, "//"):
			finishRecord()
		case section == "origin":
			for _, r := range line {
				if r >= 'A' && r <= 'Z' || r >= 'a' && r <= 'z' {
					seq = append(seq, byte(r))
				}
			}
		case section == "features" && len(line) > 21:
			body := strings.TrimSpace(line[21:])
			if line[5] != ' ' { // a new feature key starts in column 6
				finishFeature()
				feats = append(feats, Feature{Kind: strings.TrimSpace(line[:21])})
				location = body
			} else if len(feats) == 0 {
				// ignore
			} else if strings.HasPrefix(body, "/") {
				finishFeature()
				key, value := splitQualifier(body)
				qualifier = key
				f := &feats[len(feats)-1]
				switch key {
				case "gene", "locus_tag", "label":
					if f.Name == "" || key == "gene" {
						f.Name = value
					}
				case "product":
					f.Product = value
				case "translation":
					f.Protein = value
				case "transl_table":
					f.Table, _ = strconv.Atoi(strings.TrimSpace(value))
					if f.Table != 0 && !builtInTable(f.Table) && ref.UnsupportedTable == 0 {
						ref.UnsupportedTable = f.Table
					}
				}
			} else if qualifier != "" {
				f := &feats[len(feats)-1]
				switch qualifier {
				case "gene", "locus_tag", "label":
					f.Name += " " + strings.Trim(body, `"`)
				case "product":
					f.Product += " " + strings.Trim(body, `"`)
				case "translation":
					f.Protein += body
				}
			} else if location != "" {
				location += body
			}
		}
	}
	finishRecord()
	return ref, sc.Err()
}

func splitQualifier(s string) (string, string) {
	key, value, _ := strings.Cut(strings.TrimPrefix(s, "/"), "=")
	return key, strings.Trim(value, `"`)
}

// parseLocation turns a GenBank location such as complement(join(1..9,20..30))
// into zero-based half-open ranges plus a strand.
func parseLocation(s string) ([]Range, bool) {
	s = strings.TrimSpace(s)
	reverse := false
	for {
		open := strings.IndexByte(s, '(')
		if open < 0 || !strings.HasSuffix(s, ")") {
			break
		}
		switch s[:open] {
		case "complement":
			reverse = !reverse
		case "join", "order":
		default:
			return nil, reverse
		}
		s = s[open+1 : len(s)-1]
	}
	var ranges []Range
	for _, part := range strings.Split(s, ",") {
		if strings.Contains(part, "(") { // a mixed-strand join; take the inner span
			part, _ = strings.CutPrefix(part, "complement(")
			part = strings.TrimSuffix(part, ")")
		}
		from, to, found := strings.Cut(part, "..")
		if !found {
			to = from
		}
		start, err1 := strconv.Atoi(strings.Trim(from, "<>"))
		end, err2 := strconv.Atoi(strings.Trim(to, "<>"))
		if err1 != nil || err2 != nil || end < start {
			continue
		}
		ranges = append(ranges, Range{start - 1, end})
	}
	return ranges, reverse
}

func (r *Reference) sequence(contig string) []byte {
	if r == nil {
		return nil
	}
	return r.Seq[contig]
}

// covers says whether a reference has a sequence for any of these contigs.
func (r *Reference) covers(contigs []Contig) bool {
	for _, c := range contigs {
		if r.Seq[c.Name] != nil {
			return true
		}
	}
	return false
}

// MatchContigs lets a reference whose contigs are named differently still line
// up with the alignment file, by matching on sequence length.
func (r *Reference) MatchContigs(contigs []Contig) {
	if r == nil {
		return
	}
	for _, c := range contigs {
		if r.Seq[c.Name] != nil {
			continue
		}
		for name, seq := range r.Seq {
			if len(seq) == c.Length {
				r.Seq[c.Name], r.Features[c.Name], r.longest[c.Name] = seq, r.Features[name], r.longest[name]
				break
			}
		}
	}
}
