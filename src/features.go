package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// AddFeatures reads BED, GTF, or GFF3 gene annotations into a reference.
func (r *Reference) AddFeatures(path string) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	lower := strings.ToLower(path)
	var feats map[string][]Feature
	switch {
	case strings.HasSuffix(lower, ".bed"):
		feats, err = readBED(f)
	case strings.HasSuffix(lower, ".gtf"):
		feats, err = readGTF(f)
	case strings.HasSuffix(lower, ".gff"), strings.HasSuffix(lower, ".gff3"):
		feats, err = readGFF3(f)
	default:
		return fmt.Errorf("not a BED, GTF, or GFF3 file")
	}
	if err != nil {
		return err
	}
	r.addFeatures(feats)
	return nil
}

// addFeatures merges annotations and keeps each contig's list sorted.
func (r *Reference) addFeatures(feats map[string][]Feature) {
	for contig, list := range feats {
		longest := 0
		for _, f := range list {
			longest = max(longest, f.End()-f.Start())
			if f.Table != 0 && !builtInTable(f.Table) && r.UnsupportedTable == 0 {
				r.UnsupportedTable = f.Table
			}
		}
		r.mergeFeatures(contig, list, longest)
	}
}

func (r *Reference) mergeFeatures(contig string, list []Feature, longest int) {
	r.Features[contig] = append(r.Features[contig], list...)
	sort.SliceStable(r.Features[contig], func(i, j int) bool {
		return r.Features[contig][i].Start() < r.Features[contig][j].Start()
	})
	r.longest[contig] = max(r.longest[contig], longest)
}

func featureScanner(r io.Reader) *bufio.Scanner {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	return sc
}

func firstOf(values ...string) string {
	for _, v := range values {
		if v != "" && v != "." {
			return v
		}
	}
	return ""
}

func readBED(r io.Reader) (map[string][]Feature, error) {
	feats := map[string][]Feature{}
	sc := featureScanner(r)
	for sc.Scan() {
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") ||
			strings.HasPrefix(text, "track") || strings.HasPrefix(text, "browser") {
			continue
		}
		c := strings.Split(text, "\t")
		if len(c) < 3 {
			continue
		}
		start, err1 := strconv.Atoi(c[1])
		end, err2 := strconv.Atoi(c[2])
		if err1 != nil || err2 != nil || end <= start {
			continue
		}
		name := ""
		if len(c) > 3 && c[3] != "." {
			name = c[3]
		}
		feats[c[0]] = append(feats[c[0]], Feature{
			Kind: "gene", Name: name, Ranges: []Range{{start, end}},
			Reverse: len(c) > 5 && c[5] == "-",
		})
	}
	return feats, sc.Err()
}

func readGTF(r io.Reader) (map[string][]Feature, error) {
	feats := map[string][]Feature{}
	type transcriptKey struct{ contig, transcript string }
	cds := map[transcriptKey]*Feature{}
	var order []transcriptKey
	sc := featureScanner(r)
	for sc.Scan() {
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") {
			continue
		}
		c := strings.Split(text, "\t")
		if len(c) < 9 {
			continue
		}
		start, err1 := strconv.Atoi(c[3])
		end, err2 := strconv.Atoi(c[4])
		if err1 != nil || err2 != nil || end < start {
			continue
		}
		kind := c[2]
		attrs := parseGTFAttributes(c[8])
		name := firstOf(attrs["gene_name"], attrs["gene_id"], attrs["transcript_id"], kind)
		if kind == "CDS" {
			key := transcriptKey{c[0], firstOf(attrs["transcript_id"], attrs["gene_id"], name)}
			if f := cds[key]; f != nil {
				f.Ranges = append(f.Ranges, Range{start - 1, end})
			} else {
				f = &Feature{Kind: kind, Name: name, Product: attrs["product"], Reverse: c[6] == "-", Ranges: []Range{{start - 1, end}}}
				f.Table, _ = strconv.Atoi(attrs["transl_table"])
				cds[key] = f
				order = append(order, key)
			}
			continue
		}
		if kind == "source" || kind == "chromosome" {
			continue
		}
		feats[c[0]] = append(feats[c[0]], Feature{
			Kind: kind, Name: name, Product: attrs["product"],
			Reverse: c[6] == "-", Ranges: []Range{{start - 1, end}},
		})
	}
	for _, key := range order {
		f := cds[key]
		sort.Slice(f.Ranges, func(i, j int) bool { return f.Ranges[i].Start < f.Ranges[j].Start })
		feats[key.contig] = append(feats[key.contig], *f)
	}
	return feats, sc.Err()
}

func parseGTFAttributes(s string) map[string]string {
	attrs := map[string]string{}
	for _, field := range strings.Split(s, ";") {
		field = strings.TrimSpace(field)
		if field == "" {
			continue
		}
		key, value, ok := strings.Cut(field, " ")
		if ok {
			attrs[key] = strings.Trim(value, `"`)
		}
	}
	return attrs
}

// gffRow is one line of a GFF3 file, with its 1-based positions made
// zero-based and half-open.
type gffRow struct {
	contig  string
	kind    string
	start   int
	end     int
	reverse bool
	attrs   map[string]string
}

func readGFF3(r io.Reader) (map[string][]Feature, error) {
	var rows []gffRow
	byID := map[string]int{}
	sc := featureScanner(r)
	for sc.Scan() {
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") {
			continue
		}
		c := strings.Split(text, "\t")
		if len(c) < 9 {
			continue
		}
		start, err1 := strconv.Atoi(c[3])
		end, err2 := strconv.Atoi(c[4])
		if err1 != nil || err2 != nil || end < start {
			continue
		}
		at := gffRow{c[0], c[2], start - 1, end, c[6] == "-", parseGFFAttributes(c[8])}
		rows = append(rows, at)
		if id := at.attrs["ID"]; id != "" {
			byID[id] = len(rows) - 1
		}
	}

	feats := map[string][]Feature{}
	type transcriptKey struct{ contig, transcript string }
	cds := map[transcriptKey]*Feature{}
	var order []transcriptKey
	for i, at := range rows {
		if at.kind == "source" || at.kind == "region" || at.kind == "chromosome" {
			continue
		}
		if at.kind == "CDS" {
			key := transcriptKey{at.contig, firstOf(firstParent(at.attrs["Parent"]), at.attrs["ID"], at.attrs["gene"], gffName(at))}
			if f := cds[key]; f != nil {
				f.Ranges = append(f.Ranges, Range{at.start, at.end})
			} else {
				f = &Feature{
					Kind: "CDS", Name: gffGeneName(rows, byID, i), Product: at.attrs["product"],
					Reverse: at.reverse, Ranges: []Range{{at.start, at.end}},
				}
				f.Table, _ = strconv.Atoi(at.attrs["transl_table"])
				cds[key] = f
				order = append(order, key)
			}
			continue
		}
		feats[at.contig] = append(feats[at.contig], Feature{
			Kind: at.kind, Name: gffName(at), Product: at.attrs["product"],
			Reverse: at.reverse, Ranges: []Range{{at.start, at.end}},
		})
	}
	for _, key := range order {
		f := cds[key]
		sort.Slice(f.Ranges, func(i, j int) bool { return f.Ranges[i].Start < f.Ranges[j].Start })
		feats[key.contig] = append(feats[key.contig], *f)
	}
	return feats, sc.Err()
}

func parseGFFAttributes(s string) map[string]string {
	attrs := map[string]string{}
	for _, field := range strings.Split(s, ";") {
		key, value, ok := strings.Cut(strings.TrimSpace(field), "=")
		if ok {
			attrs[key] = value
		}
	}
	return attrs
}

func firstParent(s string) string {
	if i := strings.IndexByte(s, ','); i >= 0 {
		s = s[:i]
	}
	return strings.TrimSpace(s)
}

func gffName(at gffRow) string {
	return firstOf(at.attrs["Name"], at.attrs["gene"], at.attrs["ID"], at.attrs["Parent"], at.kind)
}

// gffGeneName walks up a CDS's parents to the gene it belongs to.
func gffGeneName(rows []gffRow, byID map[string]int, at int) string {
	name := gffName(rows[at])
	parent := firstParent(rows[at].attrs["Parent"])
	for depth := 0; parent != "" && depth < 8; depth++ {
		i, ok := byID[parent]
		if !ok {
			break
		}
		up := rows[i]
		if up.kind == "gene" {
			return gffName(up)
		}
		if gene := up.attrs["gene"]; gene != "" {
			return gene
		}
		if up.attrs["Name"] != "" {
			name = up.attrs["Name"]
		}
		parent = firstParent(up.attrs["Parent"])
	}
	return name
}
