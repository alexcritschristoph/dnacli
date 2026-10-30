# dnacli
A modern BAM/SAM terminal viewer.

`dnacli` is an interactive read mapping alignment viewer for the terminal. 

<img width="1770" height="930" alt="ecoli-genes-three-quarters" src="https://alexcritschristoph.github.io/dnacli/ecoli-genes-three-quarters.gif" />


`dnacli` is inspired by `samtools tview`, but has many improvements and quality of life features: full pan and zoom; can open multiple bams for comparison; visualizes reads on genes and amino acid translations; highlights mutations and their effects; quickly jumps to mutations in a vcf; switching between contigs, and more.

Use dnacli to open and inspect BAM or SAM files. Run it like this: 

`dnacli reads.bam [more.bam ...] [reference.gb|reference.fa] [variants.vcf]`

for example: `dnacli reads.bam reference.gb`

---

### Installation

Instant install into `~/.local/bin`:

```
curl -fsSL https://raw.githubusercontent.com/alexcritschristoph/dnacli/main/install.sh | sh
```

Or download the folder for your computer from the
[latest release](https://github.com/alexcritschristoph/dnacli/releases/latest),
unpack it, and run the binary inside:

| System | Download |
| --- | --- |
| Linux, Intel or AMD | [dnacli-linux-amd64.tar.gz](https://github.com/alexcritschristoph/dnacli/releases/latest/download/dnacli-linux-amd64.tar.gz) |
| Linux, ARM | [dnacli-linux-arm64.tar.gz](https://github.com/alexcritschristoph/dnacli/releases/latest/download/dnacli-linux-arm64.tar.gz) |
| Mac, Apple silicon | [dnacli-darwin-arm64.tar.gz](https://github.com/alexcritschristoph/dnacli/releases/latest/download/dnacli-darwin-arm64.tar.gz) |
| Mac, Intel | [dnacli-darwin-amd64.tar.gz](https://github.com/alexcritschristoph/dnacli/releases/latest/download/dnacli-darwin-amd64.tar.gz) |

```
tar xzf ~/Downloads/dnacli-linux-amd64.tar.gz
~/Downloads/dnacli-linux-amd64/dnacli reads.bam reference.gb
```

Note: A file fetched in a browser is quarantined by macOS the first time it runs, so open it once with a right click Open, or use the installer above.


### Features

#### 1. Compare multiple samples

You can open up several BAMs of samples mapped to the same reference at once. You can then navigate to a given position, and then quickly view that position across each BAM:

<img width="1180" height="620" alt="compare-samples-half" src="https://alexcritschristoph.github.io/dnacli/compare-samples-half.gif" />

#### 2. Browse sites in a VCF

You can include a vcf file with SNPs of interest. The VCF is shown in a table view, and you can select variants to instantly jump to their pileups for inspection:

<img width="1180" height="620" alt="vcholerae-variants-half" src="https://alexcritschristoph.github.io/dnacli/vcholerae-variants-half.gif" />

#### 3. Filter reads, copy reads, view mutational effects

You can filter the reads in the view by many metrics: mismatches, mapq, phred, length. You can copy the sequence of any read to clipboard. 

Alleles at >5% frequency (>=3x reads) are flagged, and if they occur in protein coding genes the codon and AA changes are shown:

<img width="1770" height="930" alt="filters-three-quarters" src="https://alexcritschristoph.github.io/dnacli/filters-three-quarters.gif" />

#### 4. Export the view to a  SVG

You can export any dnacli view as a beautiful SVG screenshot:

<img width="1584" height="424" alt="test6" src="https://github.com/user-attachments/assets/5828fabb-a882-4665-b513-752ab2ac550f" />

### Building from source

```
git clone https://github.com/alexcritschristoph/dnacli
cd dnacli
go build -o dnacli ./src
```

Use Go 1.25 or newer.

### Notes 

`dnacli` is primarily intended for use in microbial genomics (bacterial, viral, plasmids, synthetic biology). It is not especially designed for the human genome, and I have not tested it with the human genome.

This program was made with agentic programming, with ample human planning and guidance. I apologize if you run into any bugs or issues that may result from that.

### Contributing

If you're interested in contributing to dnacli, submit a pull request. I want to maintain the program's performance so I will extensively test any PR's impact on performance.


### License

Copyright (C) 2026 Alex Crits-Christoph.

dnacli is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, version 3. It comes with no warranty. See [LICENSE](LICENSE).
