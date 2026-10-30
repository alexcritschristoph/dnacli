package main

import (
	"fmt"
	"math"
	"os"
	"runtime"
	"sort"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// Alignment is one read, laid out against the reference so that column
// Start+i of the reference shows Bases[i].
type Alignment struct {
	Name       string
	Start, End int
	MapQ       byte
	Flags      sam.Flags
	MateStart  int
	MateContig string
	InsertSize int
	Bases      []byte // '-' marks a deleted base, ' ' a skipped region; nil when zoomed out
	Sequence   string // the read's own sequence, as the file stores it
	Quals      []byte
	Insertions map[int]string // bases inserted just after a reference position
	MeanQual   int            // average PHRED score of the read's bases
	Mismatches int
	HasIndel   bool  // the read carries an insertion or a deletion
	MismatchAt []int // where the read differs from the reference, in order
	// The base and score at each MismatchAt position, kept instead of every
	// base of the read when the view is too wide to draw them.
	MismatchBase  []byte
	MismatchQual  []byte
	ConsensusAt   []int  // where the read differs from the consensus, filled in when drawn
	consensusDone bool   // whether ConsensusAt has been worked out yet
	LeadClip      int    // soft-clipped bases drawn before the alignment starts
	TailClip      int    // and after it ends
	Cigar         string // written out only when the view is close enough to show bases
	Row           int    // filled in by the pileup packer
}

func (a *Alignment) Reverse() bool { return a.Flags&sam.Reverse != 0 }

// SoftClipped is true for a base the aligner set aside rather than lined up
// against the reference.
func (a *Alignment) SoftClipped(pos int) bool {
	return pos < a.Start+a.LeadClip || pos >= a.End-a.TailClip
}
func (a *Alignment) Paired() bool { return a.Flags&sam.Paired != 0 }

// PairedWith is true when b is this read, or the other read of its pair. Only
// paired reads share a name; an unpaired read can be called anything at all.
func (a *Alignment) PairedWith(b *Alignment) bool {
	return b == a || b != nil && a.Paired() && a.Name == b.Name
}

// BaseAt returns the read base at a reference position, or 0 outside the read.
// When the view is zoomed out only the bases that differ from the reference are
// kept, so a matching base reads as 0, which means "nothing to show".
func (a *Alignment) BaseAt(pos int) byte {
	if pos < a.Start || pos >= a.End {
		return 0
	}
	if a.Bases != nil {
		return a.Bases[pos-a.Start]
	}
	if i := sort.SearchInts(a.MismatchAt, pos); i < len(a.MismatchAt) && a.MismatchAt[i] == pos {
		return a.MismatchBase[i]
	}
	return 0
}

// deletion is the span of deleted bases that a position is part of.
func (a *Alignment) deletion(pos int) (int, int) {
	from, to := pos, pos+1
	if a.Bases == nil {
		return from, to // zoomed out, where only the fact of the deletion is kept
	}
	for from > a.Start && a.Bases[from-1-a.Start] == '-' {
		from--
	}
	for to < a.End && a.Bases[to-a.Start] == '-' {
		to++
	}
	return from, to
}

func (a *Alignment) QualAt(pos int) byte {
	if pos < a.Start || pos >= a.End {
		return 0
	}
	if a.Quals != nil {
		return a.Quals[pos-a.Start]
	}
	if i := sort.SearchInts(a.MismatchAt, pos); i < len(a.MismatchAt) && a.MismatchAt[i] == pos {
		return a.MismatchQual[i]
	}
	return 0
}

// Contig is one reference sequence named in the alignment file header.
type Contig struct {
	Name    string
	Length  int
	Reads   int64
	Profile []float64 // coverage sketch, filled in the background
}

// Filter decides which reads are worth drawing.
type Filter struct {
	MinMapQ        int
	MinMeanQual    int // average PHRED score the read's bases must reach
	MaxMismatches  int // -1 means no limit
	HideSecondary  bool
	HideDuplicates bool
	HideUnpaired   bool
	MinLength      int  // shortest read worth drawing
	SoftClips      bool // draw the bases the aligner clipped off the ends
}

func defaultFilter() Filter {
	return Filter{MaxMismatches: -1, HideSecondary: true, SoftClips: true}
}

func (f Filter) keep(a *Alignment) bool {
	switch {
	case a.Flags&sam.Unmapped != 0:
		return false
	case int(a.MapQ) < f.MinMapQ:
		return false
	case a.MeanQual < f.MinMeanQual:
		return false
	case f.MaxMismatches >= 0 && a.Mismatches > f.MaxMismatches:
		return false
	case f.HideSecondary && a.Flags&(sam.Secondary|sam.Supplementary) != 0:
		return false
	case f.HideDuplicates && a.Flags&sam.Duplicate != 0:
		return false
	case f.HideUnpaired && a.Flags&sam.ProperPair == 0:
		return false
	case len(a.Sequence) < f.MinLength:
		return false
	}
	return true
}

// Source reads alignments from a BAM or SAM file. Indexed BAM files are read
// region by region; anything else is held in memory.
type Source struct {
	Path    string
	Contigs []Contig

	file    *os.File
	reader  *bam.Reader
	index   *bam.Index
	refs    map[string]*sam.Reference
	memory  map[string][]*sam.Record
	longest int // longest read seen, used to seek within in-memory records
}

const (
	maxAlignments = 20000 // as many reads as we will hold for one screen
	sampleSlices  = 64    // how many places we sample when a region holds far more
)

func OpenAlignments(path string) (*Source, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	s := &Source{Path: path, file: f, refs: map[string]*sam.Reference{}}
	var header *sam.Header
	if strings.HasSuffix(strings.ToLower(path), ".sam") {
		r, err := sam.NewReader(f)
		if err != nil {
			return nil, err
		}
		header = r.Header()
		if len(header.Refs()) == 0 {
			f.Close()
			return nil, fmt.Errorf("has no @SQ header lines, so the contig lengths are unknown")
		}
		s.readAll(recordReader(r.Read))
	} else {
		r, err := bam.NewReader(f, min(runtime.NumCPU(), 8)) // decompress blocks in parallel
		if err != nil {
			return nil, err
		}
		header, s.reader = r.Header(), r
		s.index = readIndex(path)
		if s.index == nil {
			s.readAll(recordReader(r.Read))
			s.reader = nil
		}
	}
	for _, ref := range header.Refs() {
		s.refs[ref.Name()] = ref
		c := Contig{Name: ref.Name(), Length: ref.Len()}
		if s.index != nil {
			if stats, ok := s.index.ReferenceStats(ref.ID()); ok {
				c.Reads = int64(stats.Mapped)
			}
		} else {
			c.Reads = int64(len(s.memory[ref.Name()]))
		}
		s.Contigs = append(s.Contigs, c)
	}
	return s, nil
}

type recordReader func() (*sam.Record, error)

func (s *Source) readAll(next recordReader) {
	s.memory = map[string][]*sam.Record{}
	for {
		rec, err := next()
		if err != nil {
			break
		}
		if rec.Ref != nil && rec.Flags&sam.Unmapped == 0 {
			s.memory[rec.Ref.Name()] = append(s.memory[rec.Ref.Name()], rec)
		}
	}
	for _, recs := range s.memory {
		sort.Slice(recs, func(i, j int) bool { return recs[i].Pos < recs[j].Pos })
		for _, rec := range recs {
			s.longest = max(s.longest, rec.End()-rec.Pos)
		}
	}
}

func readIndex(path string) *bam.Index {
	for _, name := range []string{path + ".bai", strings.TrimSuffix(path, ".bam") + ".bai"} {
		f, err := os.Open(name)
		if err != nil {
			continue
		}
		idx, err := bam.ReadIndex(f)
		f.Close()
		if err == nil {
			return idx
		}
	}
	return nil
}

func (s *Source) Close() {
	if s.file != nil {
		s.file.Close()
	}
}

// each walks every record overlapping a region, in coordinate order.
func (s *Source) each(contig string, start, end int, visit func(*sam.Record) bool) error {
	if s.memory != nil {
		recs := s.memory[contig]
		first := sort.Search(len(recs), func(i int) bool { return recs[i].Pos >= start-s.longest })
		for _, rec := range recs[first:] {
			if rec.Pos >= end {
				break
			}
			if rec.End() > start && !visit(rec) {
				break
			}
		}
		return nil
	}
	ref := s.refs[contig]
	if ref == nil {
		return nil
	}
	chunks, err := s.index.Chunks(ref, start, end)
	if err != nil {
		return nil // no chunks overlap the region
	}
	it, err := bam.NewIterator(s.reader, chunks)
	if err != nil {
		return err
	}
	defer it.Close()
	for it.Next() {
		rec := it.Record()
		if rec.Pos >= end {
			break
		}
		if rec.End() > start && !visit(rec) {
			break
		}
	}
	return it.Error()
}

// Reads returns the alignments overlapping a region that pass the filter, and
// how many the region held before the filter. A region holding far more reads
// than fit on screen is sampled at points spread evenly across it, so the
// picture still covers the whole region; the last result says when that
// happened. Detail says whether each read keeps a base for every position; when
// it is false only the bases that differ from the reference are kept, which is
// all a view too wide to draw single bases needs.
func (s *Source) Reads(contig string, start, end int, filter Filter, ref *Reference, detail bool) ([]*Alignment, int, bool) {
	if s.expected(contig, end-start) <= maxAlignments {
		out, seen := s.readRegion(contig, start, end, maxAlignments, filter, ref, detail)
		return out, seen, false
	}
	// Take an equal share from slices spread evenly across the region, in one
	// pass: reading each slice on its own would decompress the blocks between
	// them again and again.
	width, per := max(1, (end-start)/sampleSlices), max(1, maxAlignments/sampleSlices)
	taken := make([]int, sampleSlices)
	var out []*Alignment
	seen := 0
	s.each(contig, start, end, func(rec *sam.Record) bool {
		seen++
		slice := clamp((rec.Pos-start)/width, 0, sampleSlices-1)
		if taken[slice] >= per {
			return true
		}
		a := expand(rec, ref, filter.SoftClips, detail)
		if filter.keep(a) {
			out = append(out, a)
			taken[slice]++
		}
		return true
	})
	return out, seen, true
}

// expected guesses how many reads a region holds from the whole-contig count.
func (s *Source) expected(contig string, width int) int {
	c := s.contig(contig)
	if c == nil || c.Length == 0 {
		return 0
	}
	return int(float64(c.Reads) * float64(width) / float64(c.Length))
}

// readRegion loads up to limit of the reads overlapping a region.
func (s *Source) readRegion(contig string, start, end, limit int, filter Filter, ref *Reference, detail bool) ([]*Alignment, int) {
	var out []*Alignment
	seen := 0
	s.each(contig, start, end, func(rec *sam.Record) bool {
		seen++
		a := expand(rec, ref, filter.SoftClips, detail)
		if filter.keep(a) {
			out = append(out, a)
		}
		return len(out) < limit
	})
	return out, seen
}

// expand lays a record's bases out against the reference coordinates. Soft
// clips, when they are wanted, are laid out either side of the alignment. When
// detail is false the bases that match the reference are not stored, so a long
// read costs only the differences and a view far out stays quick.
func expand(rec *sam.Record, ref *Reference, clips, detail bool) *Alignment {
	a := &Alignment{
		Name: rec.Name, Start: rec.Pos, End: rec.End(), MapQ: rec.MapQ, Flags: rec.Flags,
		MateStart: rec.MatePos, InsertSize: rec.TempLen,
	}
	if detail {
		a.Cigar = rec.Cigar.String() // only the read line and the test ask for it
	}
	if clips {
		a.LeadClip, a.TailClip = softClip(rec.Cigar, 0), softClip(rec.Cigar, len(rec.Cigar)-1)
		a.Start, a.End = a.Start-a.LeadClip, a.End+a.TailClip
	}
	if c := rec.Cigar; len(c) > 0 && c[0].Type() == sam.CigarDeletion {
		a.Start += c[0].Len() // a deletion at either end is not part of the read
	}
	if c := rec.Cigar; len(c) > 0 && c[len(c)-1].Type() == sam.CigarDeletion {
		a.End -= c[len(c)-1].Len()
	}
	if rec.MateRef != nil {
		a.MateContig = rec.MateRef.Name()
	}
	seq, quals := rec.Seq.Expand(), rec.Qual // Expand gives the bases in capitals
	a.Sequence = string(seq)
	sum, n := 0, 0
	for _, q := range quals {
		if q != 0xff { // 255 stands in for a base the file gives no quality
			sum += int(q)
			n++
		}
	}
	if n > 0 {
		a.MeanQual = sum / n
	}
	if detail {
		a.Bases = make([]byte, a.End-a.Start)
		a.Quals = make([]byte, a.End-a.Start)
	}
	contig := ""
	if rec.Ref != nil {
		contig = rec.Ref.Name()
	}
	var refSeq []byte
	if ref != nil {
		refSeq = ref.Seq[contig] // read once, not through a map lookup per base
	}
	pos, query := rec.Pos, 0 // pos walks the reference, which the clipped bases sit outside of
	for _, op := range rec.Cigar {
		n := op.Len()
		if op.Type().Consumes().Query == 1 && query+n > len(seq) {
			break // the CIGAR and the sequence disagree; draw what we have
		}
		switch op.Type() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for i := 0; i < n; i++ {
				base, qual := seq[query+i], byte(0xff)
				if query+i < len(quals) {
					qual = quals[query+i]
				}
				if detail {
					a.Bases[pos-a.Start+i], a.Quals[pos-a.Start+i] = base, qual
				}
				r := byte(0)
				if at := pos + i; at >= 0 && at < len(refSeq) {
					r = upper(refSeq[at])
				}
				if r != 0 && base != 'N' && r != 'N' && base != r {
					a.Mismatches++
					a.MismatchAt = append(a.MismatchAt, pos+i)
					if !detail {
						a.MismatchBase, a.MismatchQual = append(a.MismatchBase, base), append(a.MismatchQual, qual)
					}
				}
			}
			pos, query = pos+n, query+n
		case sam.CigarDeletion:
			a.HasIndel = true
			if pos >= a.Start && pos < a.End {
				if detail {
					for i := 0; i < n; i++ {
						a.Bases[pos-a.Start+i] = '-'
					}
				}
				a.Mismatches, a.MismatchAt = a.Mismatches+1, append(a.MismatchAt, pos)
				if !detail {
					a.MismatchBase, a.MismatchQual = append(a.MismatchBase, '-'), append(a.MismatchQual, 0xff)
				}
			}
			pos += n
		case sam.CigarSkipped:
			if detail {
				for i := 0; i < n; i++ {
					a.Bases[pos-a.Start+i] = ' '
				}
			}
			pos += n
		case sam.CigarInsertion:
			a.HasIndel = true
			if detail {
				if a.Insertions == nil {
					a.Insertions = map[int]string{}
				}
				a.Insertions[max(pos-1, a.Start)] = strings.ToUpper(string(seq[query : query+n]))
			}
			query += n
			a.Mismatches++
		case sam.CigarSoftClipped:
			at := pos // a trailing clip sits after the alignment, a leading one before it
			if query == 0 {
				at = pos - n
			}
			if detail {
				for i := 0; clips && i < n; i++ {
					if k := at - a.Start + i; k >= 0 && k < len(a.Bases) {
						a.Bases[k] = seq[query+i]
						if query+i < len(quals) {
							a.Quals[k] = quals[query+i]
						}
					}
				}
			}
			query += n
		}
	}
	if ref == nil || ref.Seq[contig] == nil {
		if v, ok := rec.Tag([]byte("NM")); ok {
			if n, ok := auxInt(v); ok {
				a.Mismatches = n
			}
		}
	}
	return a
}

// softClip is the length of the soft clip at one end of a CIGAR, if there is one.
func softClip(cigar sam.Cigar, at int) int {
	if at >= 0 && at < len(cigar) && cigar[at].Type() == sam.CigarSoftClipped {
		return cigar[at].Len()
	}
	return 0
}

func auxInt(v sam.Aux) (int, bool) {
	switch n := v.Value().(type) {
	case int:
		return n, true
	case int8:
		return int(n), true
	case int16:
		return int(n), true
	case int32:
		return int(n), true
	case uint8:
		return int(n), true
	case uint16:
		return int(n), true
	case uint32:
		return int(n), true
	}
	return 0, false
}

// CoverageProfile measures mean read depth in equal bins across a whole
// contig. It reads every record, so it belongs in the background.
func (s *Source) CoverageProfile(contig string, bins int) []float64 {
	c := s.contig(contig)
	if c == nil || c.Length == 0 {
		return nil
	}
	profile := make([]float64, bins)
	width := float64(c.Length) / float64(bins)
	s.each(contig, 0, c.Length, func(rec *sam.Record) bool {
		from, to := float64(rec.Pos), float64(rec.End()) // End walks the CIGAR, so ask it once
		first, last := int(from/width), int((to-1)/width)
		for b := max(first, 0); b <= last && b < bins; b++ {
			lo, hi := float64(b)*width, float64(b+1)*width
			profile[b] += (math.Min(to, hi) - math.Max(from, lo)) / width
		}
		return true
	})
	return profile
}

func (s *Source) contig(name string) *Contig {
	for i := range s.Contigs {
		if s.Contigs[i].Name == name {
			return &s.Contigs[i]
		}
	}
	return nil
}
