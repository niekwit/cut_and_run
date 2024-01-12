rule fastqc:
    input:
        "reads/{sample}{end}_001.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}{end}.html",
        zip="results/qc/fastqc/{sample}{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"],
        mem_mb=2048
    wrapper:
        "v3.3.3/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1","_R2"])
    output:
        "results/qc/multiqc.html",
        "results/qc/multiqc_data/multiqc_general_stats.txt",
    params:
        extra="",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "multiqc " 
        "--force "
        "-o results/qc "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "2> {log}"

