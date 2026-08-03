# ----------------------------------------------------------------
# Transcript assembly and abundance estimation

rule run_stringtie:
    input:
        rules.minimap_mapping.output
    output:
        gtf=stringtie_gtf,
        abundance=gene_abundance
    params:
        opts=config.get("stringtie_opts", ""),
        strand=stringtie_strand,
        annotation=(
            f"-G {shlex.quote(str(in_annotation))}"
            if use_guide_annotation
            else ""
        ),
        benchmark_dir=output_path("benchmarks", "stringtie")
    log:
        output_path("logs", "stringtie", f"{sample}.log")
    benchmark:
        output_path("benchmarks", "stringtie", f"{sample}.tsv")
    threads: stringtie_threads
    conda:
        "../envs/stringtie.yml"
    shell:
        """
        mkdir -p \
            "$(dirname {output.gtf:q})" \
            "$(dirname {log:q})" \
            {params.benchmark_dir:q}

        stringtie \
            -L \
            -v \
            -p {threads} \
            {params.strand} \
            {params.annotation} \
            {params.opts} \
            -A {output.abundance:q} \
            -o {output.gtf:q} \
            {input:q} \
            2> {log:q}
        """
