package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

// Variant is one record from a VCF file. Start is zero-based.
type Variant struct {
	Contig string
	Start  int
	ID     string
	Ref    string
	Alt    string
	Qual   string
	Filter string
	Impact Impact
}

// End is the last reference position the variant touches, exclusive.
func (v Variant) End() int { return v.Start + max(1, len(v.Ref)) }

// Details are the VCF's own fields worth showing, skipping the empty ones.
func (v Variant) Details() []string {
	var out []string
	if present(v.ID) {
		out = append(out, v.ID)
	}
	if present(v.Qual) {
		out = append(out, "quality: "+v.Qual)
	}
	if present(v.Filter) {
		out = append(out, v.Filter)
	}
	return out
}

// present is false for a VCF field left empty or written as a full stop.
func present(field string) bool { return field != "" && field != "." }

func (v Variant) Label() string { return label(v.Ref, v.Alt) }

// label writes one allele changing into another, short enough for a status line.
func label(ref, alt string) string {
	if len(ref) > 6 {
		ref = ref[:5] + "…"
	}
	if len(alt) > 6 {
		alt = alt[:5] + "…"
	}
	return ref + "→" + alt
}

func LoadVariants(path string) ([]Variant, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var out []Variant
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		c := strings.Split(line, "\t")
		if len(c) < 5 {
			continue
		}
		pos, err := strconv.Atoi(c[1])
		if err != nil {
			continue
		}
		v := Variant{Contig: c[0], Start: pos - 1, ID: c[2], Ref: strings.ToUpper(c[3]), Alt: strings.ToUpper(c[4])}
		if len(c) > 5 {
			v.Qual = c[5]
		}
		if len(c) > 6 {
			v.Filter = c[6]
		}
		out = append(out, v)
	}
	return out, sc.Err()
}
