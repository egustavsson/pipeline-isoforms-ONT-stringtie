# ----------------------------------------------------------------
# Restrander preprocessing
#
# Loaded only when config contains:
# read_preprocessing: "restrander"

rule restrander:
    input:
        fq=rules.concatenate_reads.output,
        config=str(restrander_config)
    output:
        fq=temp(restranded_fastq)
    params:
        executable=config.get("restrander_executable", "restrander")
    log:
        output_path("logs", "restrander", f"{sample}.log")
    threads: 1
    conda:
        "../envs/restrander.yml"
    shell:
        """
        mkdir -p "$(dirname {output.fq:q})" "$(dirname {log:q})"

        {params.executable:q} \
            {input.fq:q} \
            {output.fq:q} \
            {input.config:q} \
            > {log:q} \
            2>&1
        """
