import os
import Bio
import shutil
from os import path
from re import search
from pathlib import Path

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("7.18")

# ----------------------------------------------------------------

configfile: "config.yml"
validate(config, schema="schema/config_schema.yaml")
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

shutil.copy2(SNAKEDIR + "/config.yml", WORKDIR)

sample = config["sample_name"]


in_fastq = config["reads_fastq"]
if not path.isabs(in_fastq):
    in_fastq = path.join(SNAKEDIR, in_fastq)
    assert os.path.exists(in_fastq)

in_genome = config["genome"]
if not path.isabs(in_genome):
    in_genome = path.join(SNAKEDIR, in_genome)
    assert os.path.exists(in_genome)

in_annotation = config["annotation"]
if not path.isabs(in_annotation) and in_annotation != "":
    in_annotation = path.join(SNAKEDIR, in_annotation)


# ----------------------------------------------------------------

target_list = [
    "Nanostat/stat_out.txt",
    path.join("StringTie", f"{sample}_stringtie.gff")
]


rule all:
    input:
        target_list

# ----------------------------------------------------------------

rule nanostat:
    input:
        fq=config["reads_fastq"]

    output:
        ns="Nanostat/stat_out.txt"

    threads: config["threads"]

    shell:
        """
        NanoStat -n {output.ns} -t {threads} --tsv --fastq {input.fq}
        """


# ----------------------------------------------------------------

rule pychopper:
    input:
        fq=config["reads_fastq"]

    output:
        pyfq=path.join("Pychopper", f"{sample}_full_length_reads.fastq")

    params:
        outpath="Pychopper",
        prefix=f"{sample}_full_length_reads.fastq",
        pc="True" if config["run_pychopper"] else "False",
        pc_opts=config["pychopper_opts"],
        kit=config["kit"]

    log: path.join(WORKDIR, "Pychopper", f"{sample}_pychopped.log")

    threads: config["threads"]

    run:
        if params.pc == "True":
            shell(
                """
                cd {params.outpath};
                pychopper {params.pc_opts} -t {threads} -k {params.kit} -r report.pdf -u unclassified.fq -w rescued.fq {input.fq} {params.prefix} 2> {log}
                """
            )
        else:
            shell(
                """
                ln -s `realpath {input.fq}` {output.pyfq}
                """
            )

# ----------------------------------------------------------------

rule minimap_mapping:
    input:
        genome=config["genome"],
        fq=rules.pychopper.output.pyfq

    output:
        sam=path.join("Mapping", f"{sample}_minimap.sam")

    params:
        opts=config["minimap2_opts"]

    log: path.join("Mapping", f"{sample}_minimap.log")

    threads: config["threads"]

    shell:
        """
        minimap2 -ax splice {params.opts} --secondary=no -t {threads} {input.genome} {input.fq} \
        | samtools sort -@ {threads} -o {output.sam}
        """
# ----------------------------------------------------------------

rule aln_stats:
    input:
        sam=rules.minimap_mapping.output.sam

    output:
        tsv=path.join("Mapping", f"{sample}_alignment_stats.tsv")

    threads: config["threads"]

    shell:
        """
        samtools stats -@ {threads} {input.sam} > {output.tsv}
        """
# ----------------------------------------------------------------

str_threads = 10
if config["threads"] < str_threads:
    str_threads = config["threads"]

rule run_stringtie:
    input:
        sam=rules.minimap_mapping.output.sam,
        fa=config["genome"]

    output:
        gff = path.join("StringTie", f"{sample}_stringtie.gff")

    params:
        opts = config["stringtie_opts"],
        ann = in_annotation if config["use_guide_annotation"] else ""

    log: path.join("StringTie", f"{sample}_StringTie.log")
    
    threads: config["threads"]

    run:
        g_flag = f"-G {params.ann}" if params.ann else ""
        shell(
            f"""
            stringtie --rf {g_flag} -L -v -p {threads} {params.opts} -o {output.gff} {input.sam} 2> {log}
            """
        )


