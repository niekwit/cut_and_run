rule fastqc:
    input:
        "reads/{sample}_R{end}_001.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}_{end}.html",
        zip="results/qc/fastqc/{sample}_{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"],
        mem_mb = 2048,
    wrapper:
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
        input:
            expand("results/qc/fastqc/{sample}_{end}_fastqc.zip", sample=SAMPLES, end=["1","2"])
        output:
            "results/qc/multiqc/multiqc.html",
            "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",
        params:
            dir=lambda wildcards, output: os.path.dirname(output[0]),
            extra="",  # Optional: extra parameters for multiqc
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        log:
            "logs/multiqc/multiqc.log"
        conda:
            "../envs/mapping.yaml"
        shell:
            "multiqc " 
            "--force "
            "--outdir {params.dir} "
            "-n multiqc.html "
            "{params.extra} "
            "{input} "
            "> {log} 2>&1"

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