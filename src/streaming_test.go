package main

import (
	"bufio"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/biogo/hts/bam"
)

func TestStreamedUnindexedBAMMatchesIndexed(t *testing.T) {
	path := "testdata/supertest.bam"
	want, err := OpenAlignments(path)
	if err != nil {
		t.Skip(err)
	}
	defer want.Close()
	if want.index == nil {
		t.Skip("supertest.bam has no index")
	}

	// A copy without its index is the file a user gets from an aligner that
	// forgets to run samtools index.
	copyPath := filepath.Join(t.TempDir(), "supertest.bam")
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(copyPath, data, 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := openAlignments(copyPath, 0) // force the streaming path
	if err != nil {
		t.Fatal(err)
	}
	defer got.Close()
	if !got.streaming || got.memory != nil || got.reader == nil {
		t.Fatalf("streaming=%v memory=%v reader=%v, want a streamed source", got.streaming, got.memory, got.reader)
	}
	if got.readyFor(got.Contigs[0].Name) {
		t.Fatal("a streamed source says it can answer regions before its sketch")
	}

	got.reader.Omit(bam.AllVariableLengthData)
	begin := time.Now()
	profiles, checkpoints, counts := got.SketchAll()
	t.Logf("streamed sketch of %d reads: %v", counts[got.Contigs[0].Name], time.Since(begin).Round(time.Microsecond))
	got.checkpoints = checkpoints
	for i, c := range got.Contigs {
		got.Contigs[i].Reads = counts[c.Name]
	}
	if !got.readyFor(got.Contigs[0].Name) {
		t.Fatal("the sketch did not make the streamed source ready")
	}
	if len(profiles) == 0 || len(checkpoints) == 0 {
		t.Fatalf("the sketch made %d profiles and %d checkpoint lists", len(profiles), len(checkpoints))
	}

	for i, c := range got.Contigs {
		if c.Reads != want.Contigs[i].Reads {
			t.Fatalf("%s has %d mapped reads, want %d", c.Name, c.Reads, want.Contigs[i].Reads)
		}
	}
	for _, index := range []int{0, len(got.Contigs) - 1} {
		c := got.Contigs[index]
		for _, start := range []int{0, c.Length / 3, c.Length / 2, max(0, c.Length-300)} {
			end := min(c.Length, start+300)
			a, seenA, _ := got.Reads(c.Name, start, end, defaultFilter(), nil, true)
			b, seenB, _ := want.Reads(c.Name, start, end, defaultFilter(), nil, true)
			if seenA != seenB || len(a) != len(b) {
				t.Fatalf("%s:%d-%d streamed %d of %d reads, indexed %d of %d",
					c.Name, start, end, len(a), seenA, len(b), seenB)
			}
			for j := range a {
				if a[j].Name != b[j].Name || a[j].Start != b[j].Start {
					t.Fatalf("%s:%d read %d is %s at %d, want %s at %d",
						c.Name, start, j, a[j].Name, a[j].Start, b[j].Name, b[j].Start)
				}
			}
		}
	}
}

func TestStreamedModelWaitsForItsSketch(t *testing.T) {
	path := filepath.Join(t.TempDir(), "supertest.bam")
	data, err := os.ReadFile("testdata/supertest.bam")
	if err != nil {
		t.Skip(err)
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	src, err := openAlignments(path, 0)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	m := newModel(src, nil, "", nil)
	m.width, m.height = 150, 44
	if m.showReads() {
		t.Fatal("a streamed model draws reads before its sketch is ready")
	}
	if summary := strings.Join(m.filterSummary(), " "); !strings.Contains(summary, "Measuring coverage") {
		t.Fatalf("before the sketch the filter line says %q", summary)
	}
	src.reader.Omit(bam.AllVariableLengthData)
	profiles, checkpoints, counts := src.SketchAll()
	src.checkpoints = checkpoints
	for i, c := range src.Contigs {
		src.Contigs[i].Reads = counts[c.Name]
	}
	m.Update(profileMsg{source: src, profiles: profiles, checkpoints: checkpoints, counts: counts})
	if !m.showReads() || len(m.reads) == 0 {
		t.Fatalf("after the sketch the model shows %d reads", len(m.reads))
	}
}

func TestStreamedSketchPublishesOneContigAtATime(t *testing.T) {
	path := filepath.Join(t.TempDir(), "supertest.bam")
	data, err := os.ReadFile("testdata/supertest.bam")
	if err != nil {
		t.Skip(err)
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	src, err := openAlignments(path, 0)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	want, err := OpenAlignments("testdata/supertest.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer want.Close()

	m := newModel(src, nil, "", nil)
	m.width, m.height = 150, 44
	begin := time.Now()
	msg, ok := sketchCoverageAt(src, 0)().(profileMsg)
	if !ok {
		t.Fatal("the streamed sketch returned no message")
	}
	first := time.Since(begin)
	_, next := m.Update(msg)
	if next == nil {
		t.Fatal("the whole stream was sketched in one message")
	}
	firstContig, secondContig := src.Contigs[0], src.Contigs[1]
	if !src.countsKnown(firstContig.Name) || src.countsKnown(secondContig.Name) {
		t.Fatalf("after one contig, known counts are %v, want only %s", src.known, firstContig.Name)
	}
	if !src.readyFor(firstContig.Name) || src.readyFor(secondContig.Name) {
		t.Fatal("only the first contig should be ready for region reads")
	}
	if src.sketched || len(m.reads) == 0 {
		t.Fatalf("after the first contig sketched=%v reads=%d", src.sketched, len(m.reads))
	}

	for next != nil {
		got := next()
		nextMsg, ok := got.(profileMsg)
		if !ok {
			t.Fatalf("the sketch returned %T", got)
		}
		_, next = m.Update(nextMsg)
	}
	if !src.sketched {
		t.Fatal("the stream never marked itself sketched")
	}
	for i, c := range src.Contigs {
		if c.Reads != want.Contigs[i].Reads {
			t.Fatalf("%s has %d reads, want %d", c.Name, c.Reads, want.Contigs[i].Reads)
		}
	}
	t.Logf("first contig in %s, whole sketch in %s", first.Round(time.Microsecond), time.Since(begin).Round(time.Millisecond))
}

func TestLargeUnindexedWarning(t *testing.T) {
	path := filepath.Join(t.TempDir(), "large.bam")
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := f.Truncate(largeAlignment + 1); err != nil {
		t.Fatal(err)
	}
	f.Close()

	var out strings.Builder
	if !confirmLargeUnindexed(path, bufio.NewReader(strings.NewReader("c\n")), &out) {
		t.Fatal("continue was refused")
	}
	warning := "Warning: you are opening a large unindexed SAM/BAM, which will be very slow. I recommend indexing the BAM first."
	if !strings.Contains(out.String(), warning) || !strings.Contains(out.String(), "[c]ontinue") {
		t.Fatalf("the warning says %q", out.String())
	}
	out.Reset()
	if confirmLargeUnindexed(path, bufio.NewReader(strings.NewReader("e\n")), &out) {
		t.Fatal("exit was accepted as continue")
	}
	if confirmLargeUnindexed(path, bufio.NewReader(strings.NewReader("continue\n")), &out) != true {
		t.Fatal("the word continue was refused")
	}
	if confirmLargeUnindexed(path, bufio.NewReader(strings.NewReader("exit\n")), &out) != false {
		t.Fatal("the word exit was accepted")
	}
	if confirmLargeUnindexed(path, bufio.NewReader(strings.NewReader("")), &out) {
		t.Fatal("an empty answer should exit")
	}

	// A large file with an index is not warned about, and nothing is read from
	// the prompt.
	indexed := filepath.Join(t.TempDir(), "large-indexed.bam")
	f, err = os.Create(indexed)
	if err != nil {
		t.Fatal(err)
	}
	if err := f.Truncate(largeAlignment + 1); err != nil {
		t.Fatal(err)
	}
	f.Close()
	data, err := os.ReadFile("testdata/supertest.bam.bai")
	if err != nil {
		t.Skip(err)
	}
	if err := os.WriteFile(indexed+".bai", data, 0o644); err != nil {
		t.Fatal(err)
	}
	out.Reset()
	if !confirmLargeUnindexed(indexed, bufio.NewReader(strings.NewReader("")), &out) || out.Len() != 0 {
		t.Fatalf("a large indexed file warned: %q", out.String())
	}
}

func TestSmallUnindexedFileStaysInMemory(t *testing.T) {
	path := filepath.Join(t.TempDir(), "supertest.bam")
	data, err := os.ReadFile("testdata/supertest.bam")
	if err != nil {
		t.Skip(err)
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	src, err := OpenAlignments(path)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	if src.memory == nil || src.streaming || src.reader != nil {
		t.Fatalf("small unindexed file memory=%v streaming=%v reader=%v", src.memory != nil, src.streaming, src.reader != nil)
	}
	if !src.readyFor(src.Contigs[0].Name) {
		t.Fatal("an in-memory file is not ready for regions")
	}
	msg, ok := sketchCoverage(src)().(profileMsg)
	if !ok || msg.profiles == nil || msg.checkpoints != nil {
		t.Fatalf("memory sketch returned %#v", msg)
	}
}
