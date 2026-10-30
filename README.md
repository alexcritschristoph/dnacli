# dnacli
A modern BAM/SAM terminal viewer.

`dnacli` is an interactive read mapping alignment viewer entirely within the terminal. 

It is inspired by `samtools tview`, but has many improvements and quality of life features: full pan and zoom; can open multiple bams for comparison; visualizes reads on genes and amino acid translations; highlights mutations and their effects; quickly jumps to mutations in a vcf; switching between contigs, and more.

Use dnacli to open and inspect BAM or SAM files. Run it like this: `dnacli reads.bam [more.bam ...] [reference.gb|reference.fa] [variants.vcf]`, e.g. `dnacli reads.bam reference.gb`.
<img width="2360" height="1240" alt="ecoli-genes" src="https://github.com/user-attachments/assets/6ff11889-13b0-49b9-8a5a-30102642f62a" />


---

### Installation

### Features

1. Compare multiple samples
<img width="2360" height="1240" alt="contigs" src="https://github.com/user-attachments/assets/4007f633-34be-4a5c-8b25-e1a76f6775a2" />

2. Investigate sites in a VCF
<img width="2360" height="1240" alt="vcholerae-variants" src="https://github.com/user-attachments/assets/db37d39c-d858-4953-8aba-4f9d0a80d247" />

4. Browse by contig
<img width="2360" height="1240" alt="compare-samples" src="https://github.com/user-attachments/assets/c2b54cf0-2893-4f2c-a8ef-05e352489ef9" />

5. Filter reads, copy reads, view mutational effects
<img width="2360" height="1240" alt="filters" src="https://github.com/user-attachments/assets/bddaa98b-2f3b-4955-b549-8800e506fdd8" />

6. Export the view to a  SVG
You can export any dnacli view as a beautiful SVG screenshot:
<img width="1584" height="424" alt="test6" src="https://github.com/user-attachments/assets/5828fabb-a882-4665-b513-752ab2ac550f" />

### Notes 

`dnacli` is primarily intended for use in microbial genomics (bacterial, viral, plasmids, synthetic biology). It is not especially designed for the human genome, and I have not tested it with the human genome.

### Contributing

If you're interested in contributing to dnacli, submit a pull request.
