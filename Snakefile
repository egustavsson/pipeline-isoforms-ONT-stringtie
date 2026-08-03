from pathlib import Path
import shlex

from snakemake.utils import min_version
from snakemake.utils import validate

min_version("9.23")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")

SNAKEDIR = Path(workflow.basedir).resolve()


def resolve_path(value):
    """Resolve paths relative to the pipeline directory."""
    value = Path(value).expanduser()
    if not value.is_absolute():
        value = SNAKEDIR / value
    return value.resolve()


WORKDIR = resolve_path(config["workdir"])
sample = config.get("sample_name", "sample")
sample_dir = WORKDIR / sample


def output_path(*parts):
    """Construct an output path inside the sample directory."""
    return str(sample_dir.joinpath(*parts))


# ----------------------------------------------------------------

in_fastq_path = resolve_path(config["reads_fastq"])

if not in_fastq_path.exists():
    raise FileNotFoundError(f"FASTQ input does not exist: {in_fastq_path}")

fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

if in_fastq_path.is_dir():
    in_fastq = sorted(
        str(file)
        for file in in_fastq_path.rglob("*")
        if file.is_file() and file.name.lower().endswith(fastq_extensions)
    )
else:
    if not in_fastq_path.name.lower().endswith(fastq_extensions):
        raise ValueError(
            "reads_fastq must be a FASTQ file ending in "
            ".fastq, .fq, .fastq.gz or .fq.gz"
        )
    in_fastq = [str(in_fastq_path)]

if not in_fastq:
    raise ValueError(f"No FASTQ files found in: {in_fastq_path}")

if len(in_fastq) > 1 and not config.get("concatenate_reads", False):
    raise ValueError(
        "Multiple FASTQ files were found, but concatenate_reads is false. "
        "Set concatenate_reads: true or provide a single FASTQ file."
    )


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


run_pychopper = config.get("run_pychopper", False)
library_type = config.get(
    "library_type",
    "cdna" if run_pychopper else "direct_rna",
)

valid_library_types = {"cdna", "direct_rna"}
if library_type not in valid_library_types:
    raise ValueError(
        "library_type must be one of: cdna, direct_rna"
    )

if run_pychopper and library_type != "cdna":
    raise ValueError(
        "Pychopper should only be enabled for ONT cDNA libraries"
    )


minimap2_presets = {
    "cdna": {
        "preset": "splice",
        "index_opts": "",
        "mapping_opts": "-ub",
    },
    "direct_rna": {
        "preset": "splice",
        "index_opts": "-k14",
        "mapping_opts": "-uf -k14",
    },
}

minimap2_preset = minimap2_presets[library_type]["preset"]
minimap2_index_opts = minimap2_presets[library_type]["index_opts"]
minimap2_mapping_opts = minimap2_presets[library_type]["mapping_opts"]
minimap2_extra_opts = config.get("minimap2_extra_opts", "--MD")


stringtie_strand = config.get("stringtie_strand", "")
if stringtie_strand not in {"", "--rf", "--fr"}:
    raise ValueError(
        'stringtie_strand must be "", "--rf" or "--fr"'
    )


max_threads = config.get("threads", 30)
mapping_threads = min(max_threads, 32)
sort_threads = min(4, max(1, mapping_threads // 4))
minimap2_threads = max(1, mapping_threads - sort_threads)
stringtie_threads = min(max_threads, 16)
qc_threads = min(max_threads, 8)


prepared_fastq = output_path(
    "processed_reads",
    f"{sample}_reads.fastq.gz",
)

if run_pychopper:
    analysis_fastq = output_path(
        "Pychopper",
        f"{sample}_full_length_reads.fastq.gz",
    )
else:
    analysis_fastq = prepared_fastq


genome_index = str(
    WORKDIR
    / "reference"
    / f"{in_genome.name}.{library_type}.mmi"
)


# ----------------------------------------------------------------

target_list = [
    output_path("config.yml"),
    output_path("QC", f"{sample}_NanoStat.tsv"),
    output_path("Mapping", f"{sample}_minimap2.sorted.bam"),
    output_path("Mapping", f"{sample}_minimap2.sorted.bam.bai"),
    output_path("Mapping", f"{sample}_alignment_stats.txt"),
    output_path("Mapping", f"{sample}_flagstat.txt"),
    output_path("StringTie", f"{sample}_stringtie.gtf"),
    output_path("StringTie", f"{sample}_gene_abundance.tsv"),
]

if run_pychopper:
    target_list.extend(
        [
            output_path("Pychopper", f"{sample}_pychopper_report.pdf"),
            output_path("Pychopper", f"{sample}_pychopper_stats.tsv"),
        ]
    )


rule all:
    input:
        target_list


# ----------------------------------------------------------------

rule save_config:
    input:
        config = str(SNAKEDIR / "config.yml")

    output:
        config = output_path("config.yml")

    log:
        output_path(
            "logs",
            "save_config",
            f"{sample}.log",
        )

    shell:
        """
        mkdir -p "$(dirname {output.config:q})"
        mkdir -p "$(dirname {log:q})"

        cp \
            {input.config:q} \
            {output.config:q} \
            > {log:q} \
            2>&1
        """


# ----------------------------------------------------------------

rule concatenate_reads:
    input:
        fq = in_fastq

    output:
        fq_concat = temp(prepared_fastq)

    params:
        number_of_files = len(in_fastq)

    log:
        output_path(
            "logs",
            "concatenate_reads",
            f"{sample}.log",
        )

    threads: 1

    shell:
        """
        mkdir -p "$(dirname {output.fq_concat:q})"
        mkdir -p "$(dirname {log:q})"

        if [[ {params.number_of_files} -eq 1 && {input.fq[0]:q} == *.gz ]]; then
            ln -s "$(realpath {input.fq[0]:q})" {output.fq_concat:q}
        else
            {{
                (
                    for fq in {input.fq:q}; do
                        case "$fq" in
                            *.gz)
                                gzip -cd -- "$fq"
                                ;;
                            *)
                                cat -- "$fq"
                                ;;
                        esac
                    done
                ) | gzip -c > {output.fq_concat:q}
            }} 2> {log:q}
        fi
        """


# ----------------------------------------------------------------

if run_pychopper:

    rule pychopper:
        input:
            fq = rules.concatenate_reads.output.fq_concat

        output:
            pyfq = analysis_fastq,
            report = output_path(
                "Pychopper",
                f"{sample}_pychopper_report.pdf",
            ),
            stats = output_path(
                "Pychopper",
                f"{sample}_pychopper_stats.tsv",
            ),
            unclassified = output_path(
                "Pychopper",
                f"{sample}_unclassified.fastq.gz",
            ),
            rescued = output_path(
                "Pychopper",
                f"{sample}_rescued.fastq.gz",
            )

        params:
            opts = config.get("pychopper_opts", ""),
            kit = config.get("kit", "PCS109")

        log:
            output_path(
                "logs",
                "pychopper",
                f"{sample}.log",
            )

        threads: min(max_threads, 16)

        shell:
            """
            mkdir -p "$(dirname {output.pyfq:q})"
            mkdir -p "$(dirname {log:q})"

            pychopper \
                {params.opts} \
                -t {threads} \
                -k {params.kit} \
                -r {output.report:q} \
                -S {output.stats:q} \
                -u {output.unclassified:q} \
                -w {output.rescued:q} \
                {input.fq:q} \
                {output.pyfq:q} \
                2> {log:q}
            """


# ----------------------------------------------------------------

rule nanostat:
    input:
        fq = analysis_fastq

    output:
        stats = output_path(
            "QC",
            f"{sample}_NanoStat.tsv",
        )

    log:
        output_path(
            "logs",
            "nanostat",
            f"{sample}.log",
        )

    threads: qc_threads

    shell:
        """
        mkdir -p "$(dirname {output.stats:q})"
        mkdir -p "$(dirname {log:q})"

        NanoStat \
            --fastq {input.fq:q} \
            --threads {threads} \
            --tsv \
            > {output.stats:q} \
            2> {log:q}
        """


# ----------------------------------------------------------------

rule minimap2_index:
    input:
        genome = str(in_genome)

    output:
        mmi = genome_index

    params:
        preset = minimap2_preset,
        opts = minimap2_index_opts

    log:
        output_path(
            "logs",
            "minimap2_index",
            f"{sample}.log",
        )

    threads: 1

    shell:
        """
        mkdir -p "$(dirname {output.mmi:q})"
        mkdir -p "$(dirname {log:q})"

        minimap2 \
            -x {params.preset} \
            {params.opts} \
            -d {output.mmi:q} \
            {input.genome:q} \
            2> {log:q}
        """


# ----------------------------------------------------------------

rule minimap_mapping:
    input:
        index = rules.minimap2_index.output.mmi,
        fq = analysis_fastq

    output:
        bam = output_path(
            "Mapping",
            f"{sample}_minimap2.sorted.bam",
        )

    params:
        preset = minimap2_preset,
        protocol_opts = minimap2_mapping_opts,
        extra_opts = minimap2_extra_opts,
        minimap2_threads = minimap2_threads,
        sort_threads = sort_threads,
        sort_memory = config.get("samtools_sort_memory", "1G"),
        benchmark_dir = output_path(
            "benchmarks",
            "minimap2",
        )

    log:
        minimap2 = output_path(
            "logs",
            "minimap2",
            f"{sample}.log",
        ),
        samtools = output_path(
            "logs",
            "samtools_sort",
            f"{sample}.log",
        )

    benchmark:
        output_path(
            "benchmarks",
            "minimap2",
            f"{sample}.tsv",
        )

    threads: mapping_threads

    shell:
        """
        mkdir -p "$(dirname {output.bam:q})"
        mkdir -p "$(dirname {log.minimap2:q})"
        mkdir -p "$(dirname {log.samtools:q})"
        mkdir -p {params.benchmark_dir:q}

        minimap2 \
            -ax {params.preset} \
            {params.protocol_opts} \
            {params.extra_opts} \
            --secondary=no \
            -t {params.minimap2_threads} \
            {input.index:q} \
            {input.fq:q} \
            2> {log.minimap2:q} \
        | samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_memory} \
            -O BAM \
            -o {output.bam:q} \
            - \
            2> {log.samtools:q}
        """


# ----------------------------------------------------------------

rule index_alignment:
    input:
        bam = rules.minimap_mapping.output.bam

    output:
        bai = output_path(
            "Mapping",
            f"{sample}_minimap2.sorted.bam.bai",
        )

    log:
        output_path(
            "logs",
            "samtools_index",
            f"{sample}.log",
        )

    threads: qc_threads

    shell:
        """
        mkdir -p "$(dirname {log:q})"

        samtools index \
            -@ {threads} \
            {input.bam:q} \
            {output.bai:q} \
            2> {log:q}
        """


# ----------------------------------------------------------------

rule aln_stats:
    input:
        bam = rules.minimap_mapping.output.bam,
        bai = rules.index_alignment.output.bai

    output:
        stats = output_path(
            "Mapping",
            f"{sample}_alignment_stats.txt",
        ),
        flagstat = output_path(
            "Mapping",
            f"{sample}_flagstat.txt",
        )

    log:
        output_path(
            "logs",
            "alignment_stats",
            f"{sample}.log",
        )

    threads: qc_threads

    shell:
        """
        mkdir -p "$(dirname {output.stats:q})"
        mkdir -p "$(dirname {log:q})"

        samtools stats \
            -@ {threads} \
            {input.bam:q} \
            > {output.stats:q} \
            2> {log:q}

        samtools flagstat \
            -@ {threads} \
            {input.bam:q} \
            > {output.flagstat:q} \
            2>> {log:q}
        """


# ----------------------------------------------------------------

rule run_stringtie:
    input:
        bam = rules.minimap_mapping.output.bam

    output:
        gtf = output_path(
            "StringTie",
            f"{sample}_stringtie.gtf",
        ),
        abundance = output_path(
            "StringTie",
            f"{sample}_gene_abundance.tsv",
        )

    params:
        opts = config.get("stringtie_opts", ""),
        strand = stringtie_strand,
        annotation = (
            f"-G {shlex.quote(str(in_annotation))}"
            if use_guide_annotation
            else ""
        ),
        benchmark_dir = output_path(
            "benchmarks",
            "stringtie",
        )

    log:
        output_path(
            "logs",
            "stringtie",
            f"{sample}.log",
        )

    benchmark:
        output_path(
            "benchmarks",
            "stringtie",
            f"{sample}.tsv",
        )

    threads: stringtie_threads

    shell:
        """
        mkdir -p "$(dirname {output.gtf:q})"
        mkdir -p "$(dirname {log:q})"
        mkdir -p {params.benchmark_dir:q}

        stringtie \
            -L \
            -v \
            -p {threads} \
            {params.strand} \
            {params.annotation} \
            {params.opts} \
            -A {output.abundance:q} \
            -o {output.gtf:q} \
            {input.bam:q} \
            2> {log:q}
        """
