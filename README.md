# dnacli
A modern BAM/SAM terminal viewer.

`dnacli` is an interactive read mapping alignment viewer entirely within the terminal. 

It is inspired by `samtools tview`, but has many improvements and quality of life features: full pan and zoom; can open multiple bams for comparison; visualizes reads on genes and amino acid translations; highlights mutations and their effects; quickly jumps to mutations in a vcf; switching between contigs, and more.

Use dnacli to open and inspect BAM or SAM files. Run it like this: `dnacli reads.bam [more.bam ...] [reference.gb|reference.fa] [variants.vcf]`, e.g. `dnacli reads.bam reference.gb`.


---

### Installation


### Notes 

`dnacli` is primarily intended for use in microbial genomics (bacterial, viral, plasmids, synthetic biology). It is not especially designed for the human genome, and I have not tested it with the human genome.

### Contributing

If you're interested in contributing to dnacli, submit a pull request.
