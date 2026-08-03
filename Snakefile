from pathlib import Path
import shlex

from snakemake.utils import min_version, validate


min_version("9.23")


# -----------------------------------------------------------------------------
# Configuration and paths

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")

SNAKEDIR = Path(workflow.basedir).resolve()


def resolve_path(value):
    """Resolve user-supplied paths relative to the pipeline directory."""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = SNAKEDIR / path
    return path.resolve()


WORKDIR = resolve_path(config["workdir"])
sample = config.get("sample_name", "sample")
sample_dir = WORKDIR / sample


def output_path(*parts):
    """Construct an output path inside the sample results directory."""
    return str(sample_dir.joinpath(*parts))


# -----------------------------------------------------------------------------
# Input reads

in_fastq_path = resolve_path(config["reads_fastq"])

if not in_fastq_path.exists():
    raise FileNotFoundError(f"FASTQ input does not exist: {in_fastq_path}")

fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

if in_fastq_path.is_dir():
    in_fastq = sorted(
        str(path)
        for path in in_fastq_path.rglob("*")
        if path.is_file() and path.name.lower().endswith(fastq_extensions)
    )
else:
    if not in_fastq_path.name.lower().endswith(fastq_extensions):
        raise ValueError(
            "reads_fastq must end in .fastq, .fq, .fastq.gz or .fq.gz"
        )
    in_fastq = [str(in_fastq_path)]

if not in_fastq:
    raise ValueError(f"No FASTQ files found in: {in_fastq_path}")

if len(in_fastq) > 1 and not config.get("concatenate_reads", False):
    raise ValueError(
        "Multiple FASTQ files were found, but concatenate_reads is false. "
        "Set concatenate_reads: true or provide a single FASTQ file."
    )


# -----------------------------------------------------------------------------
# Reference files

in_genome = resolve_path(config["genome"])
if not in_genome.is_file():
    raise FileNotFoundError(f"Reference genome does not exist: {in_genome}")

use_guide_annotation = config.get("use_guide_annotation", False)
in_annotation = None

if use_guide_annotation:
    annotation = config.get("annotation", "")
    if not annotation:
        raise ValueError(
            "annotation is required when use_guide_annotation is true"
        )
    in_annotation = resolve_path(annotation)
    if not in_annotation.is_file():
        raise FileNotFoundError(
            f"Reference annotation does not exist: {in_annotation}"
        )


# -----------------------------------------------------------------------------
# Library type and read preprocessing

library_type = config.get("library_type", "cdna")
if library_type not in {"cdna", "direct_rna"}:
    raise ValueError("library_type must be one of: cdna, direct_rna")

read_preprocessing = config.get("read_preprocessing", "none")
if read_preprocessing not in {"restrander", "pychopper", "none"}:
    raise ValueError(
        "read_preprocessing must be one of: restrander, pychopper, none"
    )

if library_type == "direct_rna" and read_preprocessing != "none":
    raise ValueError(
        "Direct RNA reads must use read_preprocessing: none"
    )

restrander_config = None
if read_preprocessing == "restrander":
    configured_restrander = config.get("restrander_config", "")
    if not configured_restrander:
        raise ValueError(
            "restrander_config is required when read_preprocessing is restrander"
        )
    restrander_config = resolve_path(configured_restrander)
    if not restrander_config.is_file():
        raise FileNotFoundError(
            f"Restrander configuration does not exist: {restrander_config}"
        )

if read_preprocessing == "pychopper" and not config.get("pychopper_kit", ""):
    raise ValueError(
        "pychopper_kit is required when read_preprocessing is pychopper"
    )


# Restrander and Pychopper orient cDNA reads. Unprocessed cDNA may contain
# reads in either orientation. Direct RNA uses minimap2's direct-RNA settings.
if library_type == "direct_rna":
    minimap2_mapping_opts = "-uf -k14"
elif read_preprocessing in {"restrander", "pychopper"}:
    minimap2_mapping_opts = "-uf"
else:
    minimap2_mapping_opts = "-ub"

minimap2_extra_opts = config.get("minimap2_extra_opts", "--MD")

stringtie_strand = config.get("stringtie_strand", "")
if stringtie_strand not in {"", "--rf", "--fr"}:
    raise ValueError('stringtie_strand must be "", "--rf" or "--fr"')


# -----------------------------------------------------------------------------
# Thread allocation

max_threads = int(config.get("threads", 30))
if max_threads < 1:
    raise ValueError("threads must be at least 1")


def configured_threads(name, default):
    """Read and validate a per-tool thread setting."""
    value = int(config.get(name, default))
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
    return min(value, max_threads)


seqkit_threads = configured_threads("seqkit_threads", 8)
pychopper_threads = configured_threads("pychopper_threads", 16)
requested_minimap2_threads = configured_threads("minimap2_threads", 24)
requested_sort_threads = configured_threads("samtools_sort_threads", 4)
samtools_threads = configured_threads("samtools_threads", 8)
stringtie_threads = configured_threads("stringtie_threads", 16)

if requested_minimap2_threads + requested_sort_threads > max_threads:
    raise ValueError(
        "minimap2_threads + samtools_sort_threads must not exceed threads"
    )

mapping_threads = requested_minimap2_threads + requested_sort_threads


def mapping_minimap2_threads(allocated_threads):
    """Scale minimap2 threads if Snakemake reduces the mapping allocation."""
    if allocated_threads <= 1:
        return 1
    fraction = requested_minimap2_threads / mapping_threads
    return max(1, min(allocated_threads - 1, round(allocated_threads * fraction)))


def mapping_sort_threads(allocated_threads):
    """Assign the remaining mapping-job threads to samtools sort."""
    if allocated_threads <= 1:
        return 1
    return max(1, allocated_threads - mapping_minimap2_threads(allocated_threads))


# -----------------------------------------------------------------------------
# Workflow paths

prepared_fastq = output_path(
    "processed_reads",
    f"{sample}_reads.fastq.gz",
)
restranded_fastq = output_path(
    "Restrander",
    f"{sample}_restranded.fastq.gz",
)
pychopper_fastq = output_path(
    "Pychopper",
    f"{sample}_full_length_reads.fastq.gz",
)
pychopper_report = output_path(
    "Pychopper",
    f"{sample}_pychopper_report.pdf",
)
pychopper_stats = output_path(
    "Pychopper",
    f"{sample}_pychopper_stats.tsv",
)
pychopper_unclassified = output_path(
    "Pychopper",
    f"{sample}_unclassified.fastq.gz",
)
pychopper_rescued = output_path(
    "Pychopper",
    f"{sample}_rescued.fastq.gz",
)

if read_preprocessing == "restrander":
    analysis_fastq = restranded_fastq
elif read_preprocessing == "pychopper":
    analysis_fastq = pychopper_fastq
else:
    analysis_fastq = prepared_fastq

config_copy = output_path("config.yml")
read_stats = output_path("QC", f"{sample}_seqkit_stats.tsv")
bam = output_path("Mapping", f"{sample}_minimap2.sorted.bam")
bai = f"{bam}.bai"
alignment_stats = output_path(
    "Mapping",
    f"{sample}_alignment_stats.txt",
)
flagstat = output_path("Mapping", f"{sample}_flagstat.txt")
stringtie_gtf = output_path(
    "StringTie",
    f"{sample}_stringtie.gtf",
)
gene_abundance = output_path(
    "StringTie",
    f"{sample}_gene_abundance.tsv",
)


# -----------------------------------------------------------------------------
# Default targets

target_list = [
    config_copy,
    read_stats,
    bam,
    bai,
    alignment_stats,
    flagstat,
    stringtie_gtf,
    gene_abundance,
]

if read_preprocessing == "pychopper":
    target_list.extend([pychopper_report, pychopper_stats])


rule all:
    input:
        target_list


# -----------------------------------------------------------------------------
# Modular rules

include: "rules/common.smk"

if read_preprocessing == "restrander":
    include: "rules/restrander.smk"
elif read_preprocessing == "pychopper":
    include: "rules/pychopper.smk"

include: "rules/alignment.smk"
include: "rules/stringtie.smk"
