package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	tea "github.com/charmbracelet/bubbletea"
)

// version is what ./build.sh stamps in, so a download can say which release it is.
var version = "0.9.2"

var usage = `dnacli: a modern terminal read alignment viewer.
version ` + version + `

Usage:
  dnacli reads.bam [reference.gb|reference.fa]

You can also load a reference genbank, fasta, and VCF:
  dnacli reads.bam [reference.gb|reference.fa] [variants.vcf|variants.vcf.gz]

You can also specify multiple BAMs at once, if they share the same reference:
  dnacli reads.bam [more.bam ...]

Examples:
  dnacli sample_a.bam sample_b.bam genome.gbk snvs.vcf
  dnacli samples*.bam

Files are recognised by their extensions, and can be in any order.

BED, GTF, and GFF3 files are also supported for gene features.
`

type inputs struct {
	alignments []string
	reference  string
	features   string
	variants   string
}

// readInputs sorts the files given on the command line by their extension.
func readInputs(args []string) (inputs, error) {
	var in inputs
	for _, arg := range args {
		switch {
		case hasExt(arg, ".bam", ".sam"):
			in.alignments = append(in.alignments, arg)
		case hasExt(arg, ".gb", ".gbk", ".genbank", ".fa", ".fasta", ".fna"):
			if in.reference != "" {
				return in, fmt.Errorf("Error: only one reference file can be provided.")
			}
			in.reference = arg
		case hasExt(arg, ".bed", ".gtf", ".gff", ".gff3"):
			if in.features != "" {
				return in, fmt.Errorf("Error: only one BED, GTF, or GFF3 file can be provided.")
			}
			in.features = arg
		case hasExt(arg, ".vcf", ".vcf.gz", ".vcf.bgz"):
			if in.variants != "" {
				return in, fmt.Errorf("Error: only one VCF can be provided.")
			}
			in.variants = arg
		default:
			return in, fmt.Errorf("Error, filetype not recognized %s", arg)
		}
	}
	if in.features != "" && in.reference == "" {
		return in, fmt.Errorf("Error: %s needs a FASTA or GenBank reference.", in.features)
	}
	return in, nil
}

// confirmLargeUnindexed warns before a large unindexed file is opened. It
// returns false when the user chooses to exit.
func confirmLargeUnindexed(path string, in *bufio.Reader, out io.Writer) bool {
	info, err := os.Stat(path)
	if err != nil || info.Size() <= largeAlignment {
		return true
	}
	if readIndex(path) != nil {
		return true
	}
	fmt.Fprintln(out, "Warning: you are opening a large unindexed SAM/BAM, which will be very slow. I recommend indexing the BAM first.")
	fmt.Fprint(out, "[c]ontinue  [e]xit: ")
	line, _ := in.ReadString('\n')
	switch strings.ToLower(strings.TrimSpace(line)) {
	case "c", "continue", "y", "yes":
		return true
	}
	return false
}

func main() {
	showVersion := flag.Bool("version", false, "print the version and exit")
	flag.Usage = func() { fmt.Fprint(os.Stderr, usage) }
	flag.Parse()
	if *showVersion {
		fmt.Println(version)
		return
	}

	in, err := readInputs(flag.Args())
	if err != nil {
		quit("%v", err)
	}
	if len(in.alignments) == 0 {
		flag.Usage()
		os.Exit(1)
	}
	input := bufio.NewReader(os.Stdin)
	for _, path := range in.alignments {
		if !confirmLargeUnindexed(path, input, os.Stderr) {
			return
		}
	}
	alignPath, refPath, vcfPath := in.alignments[0], in.reference, in.variants
	src, err := OpenAlignments(alignPath)
	if err != nil {
		quit("Cannot read %s: %v", alignPath, err)
	}
	defer src.Close()
	if len(src.Contigs) == 0 {
		quit("%s has no contigs in its header.", alignPath)
	}

	var ref *Reference
	if refPath != "" {
		if ref, err = LoadReference(refPath); err != nil {
			quit("Cannot read %s: %v", refPath, err)
		}
		if in.features != "" {
			if err := ref.AddFeatures(in.features); err != nil {
				quit("Cannot read %s: %v", in.features, err)
			}
		}
		ref.MatchContigs(src.Contigs)
	}
	var variants []Variant
	if vcfPath != "" {
		if variants, err = LoadVariants(vcfPath); err != nil {
			quit("Cannot read %s: %v", vcfPath, err)
		}
	}

	m := newModel(src, ref, refPath, variants)
	m.warnIfVariantsDoNotMatch(vcfPath)
	for _, path := range in.alignments[1:] {
		tabs := len(m.tabs)
		m.openAlignments(path)
		if len(m.tabs) == tabs {
			quit("%s", m.message)
		}
	}
	m.switchTab(0)
	m.async = true // panned and zoomed screens are built off the drawing thread

	if _, err := tea.NewProgram(m, tea.WithAltScreen(), tea.WithMouseAllMotion()).Run(); err != nil {
		quit("%v", err)
	}
	for _, tab := range m.tabs[1:] {
		tab.Close()
	}
}

func hasExt(path string, exts ...string) bool {
	lower := strings.ToLower(path)
	for _, e := range exts {
		if strings.HasSuffix(lower, e) {
			return true
		}
	}
	return false
}

func quit(format string, args ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", args...)
	os.Exit(1)
}
