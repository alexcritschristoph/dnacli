package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	tea "github.com/charmbracelet/bubbletea"
)

// version is what ./build.sh stamps in, so a download can say which release it is.
var version = "0.9.0"

const usage = `dnacli: a modern terminal read alignment viewer.

Usage:
  dnacli reads.bam [reference.gb|reference.fa]

You can also load a reference genbank, fasta, and VCF:
  dnacli reads.bam [reference.gb|reference.fa] [variants.vcf]

You can also specify multiple BAMs at once, if they share the same reference:
  dnacli reads.bam [more.bam ...]

Examples:
  dnacli sample_a.bam sample_b.bam genome.gbk snvs.vcf
  dnacli samples*.bam

Files are recognised by their extensions, and can be in any order.

dnacli --version prints which release this is.
`

type inputs struct {
	alignments []string
	reference  string
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
		case hasExt(arg, ".vcf"):
			if in.variants != "" {
				return in, fmt.Errorf("Error: only one VCF can be provided.")
			}
			in.variants = arg
		default:
			return in, fmt.Errorf("Error, filetype not recognized %s", arg)
		}
	}
	return in, nil
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
