# ----------------------------------------------------------------
# Common read preparation and QC

rule save_config:
    input:
        str(SNAKEDIR / "config.yml")
    output:
        config_copy
    log:
        output_path("logs", "save_config", f"{sample}.log")
    conda:
        "../envs/common.yml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        cp {input:q} {output:q} > {log:q} 2>&1
        """


rule concatenate_reads:
    input:
        in_fastq
    output:
        temp(prepared_fastq)
    params:
        number_of_files=len(in_fastq)
    log:
        output_path("logs", "concatenate_reads", f"{sample}.log")
    threads: 1
    conda:
        "../envs/common.yml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        if [[ {params.number_of_files} -eq 1 && {input[0]:q} == *.gz ]]; then
            ln -s "$(realpath {input[0]:q})" {output:q} 2> {log:q}
        else
            {{
                for fq in {input:q}; do
                    case "$fq" in
                        *.gz) gzip -cd -- "$fq" ;;
                        *)    cat -- "$fq" ;;
                    esac
                done
            }} 2> {log:q} \
            | gzip -c > {output:q} 2>> {log:q}
        fi
        """


rule seqkit_stats:
    input:
        analysis_fastq
    output:
        read_stats
    log:
        output_path("logs", "seqkit_stats", f"{sample}.log")
    threads: seqkit_threads
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        seqkit stats \
            --all \
            --tabular \
            --threads {threads} \
            {input:q} \
            > {output:q} \
            2> {log:q}
        """
