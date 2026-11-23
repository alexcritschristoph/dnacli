package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func namedFeature(ref *Reference, contig, kind, name string) *Feature {
	for i := range ref.Features[contig] {
		f := &ref.Features[contig][i]
		if f.Kind == kind && f.Name == name {
			return f
		}
	}
	return nil
}

// TestRealEcoliAnnotations uses downloaded Ensembl, NCBI, and UCSC files against
// the real E. coli K-12 BAM and FASTA in testdata. Set DNACLI_REALDATA to the
// download folder to run it.
func TestRealEcoliAnnotations(t *testing.T) {
	dir := os.Getenv("DNACLI_REALDATA")
	if dir == "" {
		t.Skip("set DNACLI_REALDATA")
	}
	src, err := OpenAlignments("testdata/ecoli_k12.bam")
	if err != nil {
		t.Skip(err)
	}
	defer src.Close()
	length := src.Contigs[0].Length

	ref, err := LoadReference("testdata/ecoli_k12.fa")
	if err != nil {
		t.Skip(err)
	}
	if err := ref.AddFeatures(filepath.Join(dir, "GCF_000005845.2_ASM584v2_genomic.gff.gz")); err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs(src.Contigs)
	gene := namedFeature(ref, "NC_000913.3", "gene", "thrL")
	if gene == nil || gene.Start() != 189 || gene.End() != 255 {
		t.Fatalf("NCBI GFF3 thrL is %+v", gene)
	}
	if cds := namedFeature(ref, "NC_000913.3", "CDS", "thrL"); cds == nil {
		t.Fatal("NCBI GFF3 has no thrL CDS")
	}
	t.Logf("NCBI GFF3: %d features on %s", len(ref.Features["NC_000913.3"]), src.Contigs[0].Name)

	m := newModel(src, ref, "testdata/ecoli_k12.fa", nil)
	m.width, m.height = 150, 44
	m.view = viewState{150, 1}
	m.reload()
	if frame := plain(m.View()); !strings.Contains(frame, "thrL") {
		t.Fatal("the real GFF3 did not put thrL in the gene track")
	}

	ref, err = LoadReference("testdata/ecoli_k12.fa")
	if err != nil {
		t.Fatal(err)
	}
	gtf := "Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.59.gtf.gz"
	if err := ref.AddFeatures(filepath.Join(dir, gtf)); err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs([]Contig{{Name: "Chromosome", Length: length}})
	if gene := namedFeature(ref, "Chromosome", "gene", "thrL"); gene == nil {
		t.Fatal("Ensembl GTF has no thrL gene")
	}
	if cds := namedFeature(ref, "Chromosome", "CDS", "thrL"); cds == nil {
		t.Fatal("Ensembl GTF has no thrL CDS")
	}
	t.Logf("Ensembl GTF: %d features on Chromosome", len(ref.Features["Chromosome"]))

	ref, err = LoadReference("testdata/ecoli_k12.fa")
	if err != nil {
		t.Fatal(err)
	}
	gff3 := "Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.59.gff3.gz"
	if err := ref.AddFeatures(filepath.Join(dir, gff3)); err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs([]Contig{{Name: "Chromosome", Length: length}})
	if gene := namedFeature(ref, "Chromosome", "gene", "thrL"); gene == nil {
		t.Fatal("Ensembl GFF3 has no thrL gene")
	}
	if cds := namedFeature(ref, "Chromosome", "CDS", "thrL"); cds == nil {
		t.Fatal("Ensembl GFF3 has no thrL CDS")
	}
	t.Logf("Ensembl GFF3: %d features on Chromosome", len(ref.Features["Chromosome"]))

	ref, err = LoadReference("testdata/ecoli_k12.fa")
	if err != nil {
		t.Fatal(err)
	}
	bed := "E_coli_K12-MG1655_Jan_2012_gb_complete.bed.gz"
	if err := ref.AddFeatures(filepath.Join(dir, bed)); err != nil {
		t.Fatal(err)
	}
	ref.MatchContigs([]Contig{{Name: "NC_000913", Length: length}})
	if gene := namedFeature(ref, "NC_000913", "gene", "thrL"); gene == nil || gene.Start() != 189 || gene.End() != 255 {
		t.Fatalf("BED thrL is %+v", gene)
	}
	t.Logf("UCSC BED: %d features on NC_000913", len(ref.Features["NC_000913"]))

	variants, err := LoadVariants(filepath.Join(dir, "ecoli_k12.vcf.gz"))
	if err != nil {
		t.Fatal(err)
	}
	if len(variants) != 592 || variants[0].Contig != "NC_000913.3" {
		t.Fatalf("gzipped E. coli VCF has %d variants, first %+v", len(variants), variants[0])
	}
	t.Logf("E. coli VCF.gz: %d variants", len(variants))
}
