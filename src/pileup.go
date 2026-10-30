package main

import "sort"

// packRows gives every read a row of the pileup, keeping the two reads of a
// pair on the same row so that a pair is easy to follow. Reads must arrive in
// coordinate order. It returns the number of rows used.
func packRows(reads []*Alignment, gapBases, windowBases int) int {
	var groups [][]*Alignment
	pairs := map[string]int{} // the group holding each paired read's name
	for _, a := range reads {
		if i, ok := pairs[a.Name]; ok && a.Paired() {
			groups[i] = append(groups[i], a)
			continue
		}
		if a.Paired() {
			pairs[a.Name] = len(groups)
		}
		groups = append(groups, []*Alignment{a})
	}

	var rowEnd []int
	place := func(start, end int) int {
		for row, e := range rowEnd {
			if start >= e {
				rowEnd[row] = end + gapBases
				return row
			}
		}
		rowEnd = append(rowEnd, end+gapBases)
		return len(rowEnd) - 1
	}

	for _, mates := range groups {
		start, end := mates[0].Start, mates[0].End
		for _, a := range mates[1:] {
			start, end = min(start, a.Start), max(end, a.End)
		}
		if end-start > max(windowBases, 2000) { // mates too far apart to share a row
			for _, a := range mates {
				a.Row = place(a.Start, a.End)
			}
			continue
		}
		row := place(start, end)
		for _, a := range mates {
			a.Row = row
		}
	}
	return len(rowEnd)
}

// site is what the reads say about one reference position.
type site struct {
	Depth  int
	Gaps   int // reads showing a deleted base here
	Counts [4]int
	Alt    byte
	Ratio  float64
}

// Count is how many reads carried one base at this position.
// tally counts what the reads show at one position.
func tally(reads []*Alignment, pos int) site {
	var s site
	for _, a := range reads {
		switch base := a.BaseAt(pos); {
		case base == 0 || base == ' ' || a.SoftClipped(pos):
		case base == '-':
			s.Depth, s.Gaps = s.Depth+1, s.Gaps+1
		default:
			s.Depth++
			if i := baseIndex(base); i >= 0 {
				s.Counts[i]++
			}
		}
	}
	return s
}

// ranked lists the alleles, deletions included, commonest first.
func (s site) ranked() []byte {
	bases := []byte("ACGT-")
	sort.SliceStable(bases, func(i, j int) bool { return s.Count(bases[i]) > s.Count(bases[j]) })
	return bases
}

func (s site) Count(base byte) int {
	if base == '-' {
		return s.Gaps
	}
	if i := baseIndex(base); i >= 0 {
		return s.Counts[i]
	}
	return 0
}

func baseIndex(base byte) int {
	return int(baseSlots[base])
}

// baseSlots is baseIndex as a table, so that counting a screenful of reads
// does not branch once per base.
var baseSlots = func() [256]int8 {
	var table [256]int8
	for i := range table {
		table[i] = -1
	}
	table['A'], table['C'], table['G'], table['T'] = 0, 1, 2, 3
	return table
}()

const (
	minVariantSupport = 3    // the second base must be seen three times
	minVariantRatio   = 0.05 // and be this share of the reads at the site
)

// baseCounts counts how many reads carried each base at every position in a
// region, and how deep the coverage was there. Deleted bases are depth but no
// base, and soft-clipped bases were never lined up so they count for nothing.
func baseCounts(reads []*Alignment, start, end int) ([][4]int, []int) {
	counts := make([][4]int, end-start)
	depth := make([]int, end-start)
	for _, a := range reads {
		if a.Bases == nil {
			continue // zoomed out, where no base is kept to count
		}
		bases := a.Bases
		from, to := max(a.Start+a.LeadClip, start), min(a.End-a.TailClip, end)
		for pos := from; pos < to; pos++ {
			i := pos - start
			depth[i]++
			if b := baseSlots[bases[pos-a.Start]]; b >= 0 {
				counts[i][b]++
			}
		}
	}
	return counts, depth
}

// consensusBase is the base most reads carry at one position, or 0 where there
// are no reads to say. A tie takes the first base in ACGT order, so that the
// same reads always give the same answer.
func consensusBase(counts [4]int) byte {
	best := 0
	for b := 1; b < 4; b++ {
		if counts[b] > counts[best] {
			best = b
		}
	}
	if counts[best] == 0 {
		return 0
	}
	return "ACGT"[best]
}

// consensusBases gives the consensus for every position in a region, one byte
// per position from start, so that drawing only has to look the answer up.
func consensusBases(counts [][4]int) []byte {
	out := make([]byte, len(counts))
	for i, c := range counts {
		out[i] = consensusBase(c)
	}
	return out
}

// consensusDeviations lists where a read carries a base other than the one
// most reads show. They are found as the region loads, because comparing every
// base against the consensus again on every frame is what made highlighting
// lag when zoomed out.
func consensusDeviations(a *Alignment, consensus []byte, start int) []int {
	if len(a.Bases) == 0 {
		return nil // zoomed out, where the bases that match are not kept
	}
	var out []int
	from, to := max(a.Start+a.LeadClip, start), min(a.End-a.TailClip, start+len(consensus))
	for pos := from; pos < to; pos++ {
		base := a.Bases[pos-a.Start]
		if baseSlots[base] >= 0 && consensus[pos-start] != base {
			out = append(out, pos)
		}
	}
	return out
}

// callSites finds positions where the reads disagree with the reference.
func callSites(counts [][4]int, depth []int, ref *Reference, contig string, start int) map[int]site {
	sites := map[int]site{}
	if ref == nil || ref.Seq[contig] == nil || len(counts) > 100000 {
		return sites
	}
	for i, c := range counts {
		refBase := ref.Base(contig, start+i)
		best, bestCount := byte(0), 0
		for b, n := range c {
			if base := "ACGT"[b]; base != refBase && n > bestCount {
				best, bestCount = base, n
			}
		}
		if bestCount < minVariantSupport {
			continue
		}
		if ratio := float64(bestCount) / float64(depth[i]); ratio >= minVariantRatio {
			sites[start+i] = site{Depth: depth[i], Counts: c, Alt: best, Ratio: ratio}
		}
	}
	return sites
}

// depthFromProfile reads the whole-contig coverage sketch back out at the
// resolution of the screen, for views too wide to load every read.
func depthFromProfile(profile []float64, contigLength, start, bpp, columns int) ([]int, int) {
	counts := make([]int, columns)
	if len(profile) == 0 || contigLength == 0 {
		return counts, 0
	}
	width := float64(contigLength) / float64(len(profile))
	peak := 0
	for col := range counts {
		from, to := float64(start+col*bpp)/width, float64(start+(col+1)*bpp)/width
		total, n := 0.0, 0
		for b := int(from); b <= int(to) && b < len(profile); b++ {
			if b >= 0 {
				total, n = total+profile[b], n+1
			}
		}
		if n > 0 {
			counts[col] = int(total/float64(n) + 0.5)
		}
		peak = max(peak, counts[col])
	}
	return counts, peak
}

// coverage counts reads per screen column.
func coverage(reads []*Alignment, cols []columnSpec, bpp int) ([]int, int) {
	if len(cols) == 0 {
		return nil, 0
	}
	start := cols[0].pos
	depth := make([]int, (cols[len(cols)-1].pos-start)/bpp+1)
	for _, a := range reads {
		first, last := (a.Start+a.LeadClip-start)/bpp, (a.End-a.TailClip-1-start)/bpp
		for i := max(first, 0); i <= last && i < len(depth); i++ {
			depth[i]++
		}
	}
	counts, peak := make([]int, len(cols)), 0
	for col, spec := range cols {
		counts[col] = depth[(spec.pos-start)/bpp]
		peak = max(peak, counts[col])
	}
	return counts, peak
}
