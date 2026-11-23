package main

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeFeatureFile(t *testing.T, name, body string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func featureReference() *Reference {
	ref := newReference()
	ref.Seq["chr1"] = []byte(strings.Repeat("A", 200))
	return ref
}

func findFeature(ref *Reference, contig, kind string) *Feature {
	for i := range ref.Features[contig] {
		if ref.Features[contig][i].Kind == kind {
			return &ref.Features[contig][i]
		}
	}
	return nil
}

func TestBEDFeatures(t *testing.T) {
	path := writeFeatureFile(t, "genes.bed", "track name=test\n"+
		"chr1\t10\t20\tgeneA\t0\t+\n"+
		"chr1\t30\t40\tgeneB\t0\t-\n")
	ref := featureReference()
	if err := ref.AddFeatures(path); err != nil {
		t.Fatal(err)
	}
	feats := ref.Features["chr1"]
	if len(feats) != 2 {
		t.Fatalf("read %d BED features, want 2", len(feats))
	}
	if feats[0].Name != "geneA" || feats[0].Start() != 10 || feats[0].End() != 20 || feats[0].Reverse {
		t.Fatalf("first BED feature is %+v", feats[0])
	}
	if feats[1].Name != "geneB" || !feats[1].Reverse || feats[1].Start() != 30 {
		t.Fatalf("second BED feature is %+v", feats[1])
	}
}

func TestGTFFeaturesGroupCDSByTranscript(t *testing.T) {
	body := "chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\"; gene_name \"GeneA\";\n" +
		"chr1\ttest\tCDS\t10\t20\t.\t+\t0\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"GeneA\";\n" +
		"chr1\ttest\tCDS\t30\t40\t.\t+\t0\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"GeneA\";\n" +
		"chr1\ttest\texon\t10\t40\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"GeneA\";\n"
	path := writeFeatureFile(t, "genes.gtf", body)
	ref := featureReference()
	if err := ref.AddFeatures(path); err != nil {
		t.Fatal(err)
	}
	cds := findFeature(ref, "chr1", "CDS")
	if cds == nil {
		t.Fatal("no CDS feature")
	}
	if cds.Name != "GeneA" || len(cds.Ranges) != 2 {
		t.Fatalf("CDS is %+v", *cds)
	}
	if cds.Start() != 9 || cds.End() != 40 {
		t.Fatalf("CDS spans %d-%d, want 9-40", cds.Start(), cds.End())
	}
}

func TestGFF3FeaturesGroupCDSByTranscript(t *testing.T) {
	body := "##gff-version 3\n" +
		"chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name=GeneA\n" +
		"chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n" +
		"chr1\t.\tCDS\t10\t20\t.\t+\t0\tID=cds1;Parent=t1\n" +
		"chr1\t.\tCDS\t30\t40\t.\t+\t0\tID=cds2;Parent=t1\n"
	path := writeFeatureFile(t, "genes.gff3", body)
	ref := featureReference()
	if err := ref.AddFeatures(path); err != nil {
		t.Fatal(err)
	}
	cds := findFeature(ref, "chr1", "CDS")
	if cds == nil {
		t.Fatal("no CDS feature")
	}
	if cds.Name != "GeneA" || len(cds.Ranges) != 2 {
		t.Fatalf("CDS is %+v", *cds)
	}
	if cds.Start() != 9 || cds.End() != 40 {
		t.Fatalf("CDS spans %d-%d, want 9-40", cds.Start(), cds.End())
	}
	if gene := findFeature(ref, "chr1", "gene"); gene == nil || gene.Name != "GeneA" {
		t.Fatalf("gene is %+v", gene)
	}
}

func TestGzippedAnnotations(t *testing.T) {
	for _, c := range []struct{ name, body string }{
		{"genes.bed.gz", "chr1\t10\t20\tGeneA\t0\t+\n"},
		{"genes.gtf.gz", "chr1\ttest\tgene\t1\t20\t.\t+\t.\tgene_id \"g1\"; gene_name \"GeneA\";\n"},
		{"genes.gff3.gz", "##gff-version 3\nchr1\t.\tgene\t1\t20\t.\t+\t.\tID=g1;Name=GeneA\n"},
	} {
		path := filepath.Join(t.TempDir(), c.name)
		f, err := os.Create(path)
		if err != nil {
			t.Fatal(err)
		}
		gz := gzip.NewWriter(f)
		if _, err := gz.Write([]byte(c.body)); err != nil {
			t.Fatal(err)
		}
		if err := gz.Close(); err != nil {
			t.Fatal(err)
		}
		f.Close()
		ref := featureReference()
		if err := ref.AddFeatures(path); err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if gene := findFeature(ref, "chr1", "gene"); gene == nil || gene.Name != "GeneA" {
			t.Fatalf("%s: gene is %+v", c.name, gene)
		}
	}
}

func TestFeaturesFollowARenamedReference(t *testing.T) {
	fasta := writeFeatureFile(t, "ref.fa", ">chr1\n"+strings.Repeat("A", 100)+"\n")
	ref, err := LoadReference(fasta)
	if err != nil {
		t.Fatal(err)
	}
	gff := writeFeatureFile(t, "genes.gff3", "##gff-version 3\n"+
		"chr1\t.\tgene\t10\t20\t.\t+\t.\tID=g1;Name=GeneA\n")
	if err := ref.AddFeatures(gff); err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs([]Contig{{Name: "NC_TEST", Length: 100}})
	if gene := findFeature(ref, "NC_TEST", "gene"); gene == nil || gene.Name != "GeneA" {
		t.Fatalf("the renamed contig has %+v", ref.Features["NC_TEST"])
	}
}

func TestFeaturesAddToGenBank(t *testing.T) {
	ref, err := LoadReference("testdata/supertest.gb")
	if err != nil {
		t.Skip(err)
	}
	before := len(ref.Features["SUPER_A"])
	gff := writeFeatureFile(t, "extra.gff3", "##gff-version 3\n"+
		"SUPER_A\t.\tgene\t10\t20\t.\t+\t.\tID=extra;Name=Extra\n")
	if err := ref.AddFeatures(gff); err != nil {
		t.Fatal(err)
	}
	if len(ref.Features["SUPER_A"]) != before+1 {
		t.Fatalf("%d features, want %d", len(ref.Features["SUPER_A"]), before+1)
	}
}

func TestGzippedVCFMatchesPlain(t *testing.T) {
	plain, err := LoadVariants("testdata/bf_fragilis.vcf")
	if err != nil {
		t.Skip(err)
	}
	gz, err := LoadVariants("testdata/bf_fragilis.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	if len(gz) != len(plain) {
		t.Fatalf("gzipped VCF has %d variants, plain has %d", len(gz), len(plain))
	}
	for i := range plain {
		if gz[i].Contig != plain[i].Contig || gz[i].Start != plain[i].Start || gz[i].Ref != plain[i].Ref || gz[i].Alt != plain[i].Alt {
			t.Fatalf("variant %d differs: %+v vs %+v", i, gz[i], plain[i])
		}
	}
}

func TestCommandLineAcceptsFeaturesAndGzippedVCF(t *testing.T) {
	in, err := readInputs([]string{"reads.bam", "genome.fa", "genes.gff3.gz", "calls.vcf.gz"})
	if err != nil {
		t.Fatal(err)
	}
	if in.reference != "genome.fa" || in.features != "genes.gff3.gz" || in.variants != "calls.vcf.gz" {
		t.Fatalf("parsed %+v", in)
	}
	for _, args := range [][]string{
		{"genes.gff3"},
		{"reads.bam", "genes.gtf"},
		{"genome.fa", "genes.gff3", "more.bed"},
	} {
		if _, err := readInputs(args); err == nil {
			t.Fatalf("%v was accepted", args)
		}
	}
}

func TestOpenCommandAddsFeatures(t *testing.T) {
	m := ready(t, "supertest")
	before := len(m.ref.Features["SUPER_A"])
	path := writeFeatureFile(t, "extra.gff3", "##gff-version 3\n"+
		"SUPER_A\t.\tgene\t10\t20\t.\t+\t.\tID=extra;Name=Extra\n")
	run(m, "open "+path)
	if len(m.ref.Features["SUPER_A"]) != before+1 {
		t.Fatalf("%d features, want %d", len(m.ref.Features["SUPER_A"]), before+1)
	}
	if !strings.Contains(m.message, "gene features opened") {
		t.Fatalf("message is %q", m.message)
	}
}
