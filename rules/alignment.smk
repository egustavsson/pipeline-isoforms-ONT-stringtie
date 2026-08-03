# ----------------------------------------------------------------
# Genome alignment and alignment QC

rule minimap_mapping:
    input:
        genome=str(in_genome),
        fq=analysis_fastq
    output:
        bam
    params:
        protocol_opts=minimap2_mapping_opts,
        extra_opts=minimap2_extra_opts,
        minimap2_threads=lambda wildcards, threads: (
            mapping_minimap2_threads(threads)
        ),
        sort_threads=lambda wildcards, threads: (
            mapping_sort_threads(threads)
        ),
        sort_memory=config.get("samtools_sort_memory", "1G"),
        benchmark_dir=output_path("benchmarks", "minimap2")
    log:
        minimap2=output_path("logs", "minimap2", f"{sample}.log"),
        samtools=output_path("logs", "samtools_sort", f"{sample}.log")
    benchmark:
        output_path("benchmarks", "minimap2", f"{sample}.tsv")
    threads: mapping_threads
    conda:
        "../envs/alignment.yml"
    shell:
        """
        mkdir -p \
            "$(dirname {output:q})" \
            "$(dirname {log.minimap2:q})" \
            "$(dirname {log.samtools:q})" \
            {params.benchmark_dir:q}

        minimap2 \
            -ax splice \
            {params.protocol_opts} \
            {params.extra_opts} \
            --secondary=no \
            -t {params.minimap2_threads} \
            {input.genome:q} \
            {input.fq:q} \
            2> {log.minimap2:q} \
        | samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_memory} \
            -O BAM \
            -o {output:q} \
            - \
            2> {log.samtools:q}
        """


rule index_alignment:
    input:
        rules.minimap_mapping.output
    output:
        bai
    log:
        output_path("logs", "samtools_index", f"{sample}.log")
    threads: samtools_threads
    conda:
        "../envs/alignment.yml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        samtools index \
            -@ {threads} \
            {input:q} \
            {output:q} \
            2> {log:q}
        """


rule aln_stats:
    input:
        bam=rules.minimap_mapping.output,
        bai=rules.index_alignment.output
    output:
        stats=alignment_stats,
        flagstat=flagstat
    log:
        output_path("logs", "alignment_stats", f"{sample}.log")
    threads: samtools_threads
    conda:
        "../envs/alignment.yml"
    shell:
        """
        mkdir -p "$(dirname {output.stats:q})" "$(dirname {log:q})"

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
