package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// runOpen opens another file: an alignment file joins the tabs along the top,
// and a reference or a VCF is shared by every tab.
func (m *model) runOpen(arg string) {
	path := expandPath(strings.TrimSpace(arg))
	if path == "" {
		m.notify("Type a file to open after open.")
		return
	}
	switch {
	case hasExt(path, ".bam", ".sam"):
		m.openAlignments(path)
	case hasExt(path, ".gb", ".gbk", ".genbank", ".fa", ".fasta", ".fna"):
		m.openReference(path)
	case hasExt(path, ".vcf"):
		m.openVariants(path)
	default:
		m.notify("Error: wrong filetype. Please open a BAM, Genbank, or VCF.")
	}
}

// expandPath writes out a leading home folder, which the bar cannot.
func expandPath(path string) string {
	if path == "~" || strings.HasPrefix(path, "~/") {
		if home, err := os.UserHomeDir(); err == nil {
			return filepath.Join(home, path[1:])
		}
	}
	return path
}

func (m *model) openAlignments(path string) {
	src, err := OpenAlignments(path)
	if err != nil {
		m.notify("Cannot read %s: %v", path, err)
		return
	}
	if len(src.Contigs) == 0 {
		src.Close()
		m.notify("%s has no contigs in its header.", filepath.Base(path))
		return
	}
	first := m.tabs[0]
	conflict := !sameContigs(src.Contigs, first.Contigs)
	if !conflict && m.ref != nil {
		m.ref.MatchContigs(src.Contigs)
		conflict = !m.ref.covers(src.Contigs)
	}
	if conflict {
		src.Close()
		m.notify("Cannot open %s: does not share reference contigs with %s.",
			filepath.Base(path), filepath.Base(first.Path))
		return
	}
	m.tabs = append(m.tabs, src)
	m.switchTab(len(m.tabs) - 1)
	m.pending = sketchCoverage(src)
}

// sameContigs is true when every contig of a file is one the first tab has,
// so that only files built on the same reference join it.
func sameContigs(contigs, first []Contig) bool {
	for _, c := range contigs {
		found := false
		for _, f := range first {
			if f.Name == c.Name {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return true
}

func (m *model) openReference(path string) {
	ref, err := LoadReference(path)
	if err != nil {
		m.notify("Cannot read %s: %v", path, err)
		return
	}
	for _, src := range m.tabs {
		ref.MatchContigs(src.Contigs)
		if !ref.covers(src.Contigs) {
			m.notify("The reference %s does not match %s: no contig name or length is shared.",
				filepath.Base(path), filepath.Base(src.Path))
			return
		}
	}
	m.ref, m.refPath = ref, path
	m.measureImpacts()
	m.reload()
	m.notify("%s opened as the reference.", filepath.Base(path))
}

func (m *model) openVariants(path string) {
	variants, err := LoadVariants(path)
	if err != nil {
		m.notify("Cannot read %s: %v", path, err)
		return
	}
	if len(variants) == 0 {
		m.notify("%s has no variants.", filepath.Base(path))
		return
	}
	if err := variantsMatch(m.contigs(), variants); err != nil {
		m.notify("%s does not match %s: %v", filepath.Base(path), filepath.Base(m.src.Path), err)
		return
	}
	m.variants = variants
	m.variantIndex = clamp(m.variantIndex, 0, len(variants)-1)
	m.measureImpacts()
	if m.panel == "" {
		m.panel = "variants"
	}
	m.reload()
	m.notify("%d variants opened from %s.", len(variants), filepath.Base(path))
}

// warnIfVariantsDoNotMatch says so when a VCF given at startup names contigs
// or positions this alignment does not have. The variants still open, but the
// warning explains why they may not line up.
func (m *model) warnIfVariantsDoNotMatch(path string) {
	if path == "" {
		return
	}
	if err := variantsMatch(m.contigs(), m.variants); err != nil {
		m.notify("Warning: %s does not match this alignment: %v", filepath.Base(path), err)
	}
}

// variantsMatch says whether every variant names a contig this file has, at a
// position it reaches. It reads the contigs into a map, so a VCF with a
// million records against a thousand contigs costs one pass, not a billion.
func variantsMatch(contigs []Contig, variants []Variant) error {
	byName := make(map[string]Contig, len(contigs))
	for _, c := range contigs {
		byName[c.Name] = c
	}
	for _, v := range variants {
		c, found := byName[v.Contig]
		if !found {
			return fmt.Errorf("VCF has %s, not found in BAM and/or reference", v.Contig)
		}
		if v.Start >= c.Length {
			return fmt.Errorf("position %s in VCF is past the end of %s", comma(v.Start+1), c.Name)
		}
	}
	return nil
}
