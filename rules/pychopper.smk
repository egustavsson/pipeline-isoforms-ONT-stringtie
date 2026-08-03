# ----------------------------------------------------------------
# Pychopper preprocessing
#
# Loaded only when config contains:
# read_preprocessing: "pychopper"

rule pychopper:
    input:
        rules.concatenate_reads.output
    output:
        fq=temp(pychopper_fastq),
        report=pychopper_report,
        stats=pychopper_stats,
        unclassified=temp(pychopper_unclassified),
        rescued=temp(pychopper_rescued)
    params:
        kit=config["pychopper_kit"],
        opts=config.get("pychopper_opts", "")
    log:
        output_path("logs", "pychopper", f"{sample}.log")
    threads: pychopper_threads
    conda:
        "../envs/pychopper.yml"
    shell:
        """
        mkdir -p "$(dirname {output.fq:q})" "$(dirname {log:q})"

        pychopper \
            {params.opts} \
            -t {threads} \
            -k {params.kit} \
            -r {output.report:q} \
            -S {output.stats:q} \
            -u {output.unclassified:q} \
            -w {output.rescued:q} \
            {input:q} \
            {output.fq:q} \
            2> {log:q}
        """
