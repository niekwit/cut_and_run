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
        "v3.3.6/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1","_R2"])
    output:
        "results/qc/multiqc.html",
        "results/qc/multiqc_data/multiqc_general_stats.txt",
    params:
        extra="",  # Optional: extra parameters for multiqc
    threads: 1
    resources:
        runtime=15,
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


rule calculate_effective_genome_sizes:
        input:
            **calculate_effective_genome_sizes_input()
        output:
            egs="results/effective_genome_sizes/effective_genome_sizes.csv",
        params:
            genome=genome,
            remove_MT_seqs=config["remove_MT_seqs"],
        threads: 1
        resources:
            runtime=10
        conda:
            "../envs/deeptools.yaml"
        log:
            "logs/effective_genome_sizes/effective_genome_sizes.log"
        script:
            "../scripts/effective_genome_size.py"
