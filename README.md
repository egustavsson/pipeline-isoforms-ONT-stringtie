# pipeline-isoforms-ONT-stringtie

<!-- badges: start -->
![Maintainer](https://img.shields.io/badge/maintainer-egustavsson-blue)
![Generic badge](https://img.shields.io/badge/WMS-snakemake-blue.svg)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
<!-- badges: end -->

This is a `Snakemake` pipeline for processing Oxford Nanopore Technologies
(ONT) long-read RNA-sequencing data. Read preprocessing is selected in
`config.yml` and can be `restrander`, `pychopper` or `none`. The pipeline
generates read statistics using `SeqKit`, aligns reads directly to the genome
FASTA using `minimap2`, generates alignment statistics using `samtools`, and
uses `StringTie` to assemble and quantify transcripts.

# Getting Started

## Input

- ONT cDNA or direct RNA reads in FASTQ format
- Reference genome assembly in FASTA format
- For SMART-Seq mRNA Long Read data, Dorado-demultiplexed FASTQ files and the
  Takara Bio `SMART-Seq_mRNA_LR.json` Restrander configuration
- Optional reference annotation in GTF format, for example a matching
  [GENCODE GTF](https://www.gencodegenes.org/human/)

## Dependencies

- [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
- The remaining dependencies, including Snakemake, Pychopper and Restrander,
  are installed through `environment.yml`

## Installation

Clone the repository:

```bash
git clone https://github.com/egustavsson/pipeline-isoforms-ONT-stringtie.git
cd pipeline-isoforms-ONT-stringtie
```

Create and activate the Conda environment:

```bash
conda config --set channel_priority strict
conda env create -f environment.yml
conda activate ont_stringtie
```

`environment.yml` installs Restrander 1.1.1 from the third-party `genomedk`
Conda channel because Restrander is not currently distributed through
Bioconda. Its upstream repository documents source installation with `make`
as an alternative.

## Configuration

Edit `config.yml` to specify the working directory, sample name, input reads,
reference genome and analysis options.

For ONT cDNA reads:

```yaml
library_type: "cdna"
read_preprocessing: "restrander"
restrander_config: "/path/to/config/SMART-Seq_mRNA_LR.json"
```

For a compatible ONT cDNA kit using Pychopper:

```yaml
library_type: "cdna"
read_preprocessing: "pychopper"
pychopper_kit: "PCS111"
```

For ONT direct RNA reads:

```yaml
library_type: "direct_rna"
read_preprocessing: "none"
```

Only the rule file for the selected preprocessing method is loaded. With
`none`, neither Restrander nor Pychopper runs.

If `reads_fastq` points to a directory containing multiple FASTQ files, set:

```yaml
concatenate_reads: true
```

To guide StringTie assembly with an existing annotation:

```yaml
use_guide_annotation: true
annotation: "/path/to/annotation.gtf"
```

## Usage

Snakemake commands should be issued from within the pipeline directory after
activating the environment.

It is a good idea to perform a dry run before executing the pipeline:

```bash
snakemake --dry-run --cores <num_cores>
```

Run the pipeline:

```bash
snakemake \
    --software-deployment-method conda \
    --cores <num_cores> \
    --rerun-incomplete
```

If a previous run was interrupted, the same command with
`--rerun-incomplete` resumes completed outputs and reruns only incomplete
jobs.

You can visualise the processes in a directed acyclic graph:

```bash
snakemake --dag | dot -Tpng > dag.png
```

To stop a running Snakemake process, press `Ctrl+C`. If Snakemake is running in
the background, a `TERM` signal will stop the scheduling of new jobs and allow
currently running jobs to finish:

```bash
killall -TERM snakemake
```

To deactivate the environment:

```bash
conda deactivate
```

## Output

Outputs are written to `<workdir>/<sample_name>/`:

```text
<workdir>/
`-- <sample_name>/
    |-- config.yml
    |-- Restrander/
    |   `-- <sample>_restranded.fastq.gz
    |-- Pychopper/
    |   |-- <sample>_full_length_reads.fastq.gz
    |   |-- <sample>_pychopper_report.pdf
    |   |-- <sample>_pychopper_stats.tsv
    |   |-- <sample>_rescued.fastq.gz
    |   `-- <sample>_unclassified.fastq.gz
    |-- QC/
    |   `-- <sample>_seqkit_stats.tsv
    |-- Mapping/
    |   |-- <sample>_minimap2.sorted.bam
    |   |-- <sample>_minimap2.sorted.bam.bai
    |   |-- <sample>_alignment_stats.txt
    |   `-- <sample>_flagstat.txt
    |-- StringTie/
    |   |-- <sample>_stringtie.gtf
    |   `-- <sample>_gene_abundance.tsv
    |-- logs/
    `-- benchmarks/
```

Only the directory selected by `read_preprocessing` is produced.
Temporary concatenated reads are removed after successful completion.

## Workflow structure

```text
Snakefile
config.yml
environment.yml
schema/
`-- config_schema.yaml
rules/
|-- common.smk
|-- restrander.smk
|-- pychopper.smk
|-- alignment.smk
`-- stringtie.smk
```
