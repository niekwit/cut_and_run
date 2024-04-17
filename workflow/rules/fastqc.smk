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
        f"{wrapper_version}/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1","_R2"])
    output:
        "results/qc/pre_trim/multiqc.html",
        "results/qc/pre_trim/multiqc_data/multiqc_general_stats.txt",
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
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
        "-o {params.dir} "
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


use rule fastqc as fastqc_post_trim with:
    input:
        "results/trimmed/{sample}_val_{end}.fq.gz"
    output:
        html="results/qc/fastqc/{sample}{end}_post_trim.html",
        zip="results/qc/fastqc/{sample}{end}_fastqc_post_trim.zip",
    log:
        "logs/fastqc/{sample}{end}_post_trim.log"


rule multiqc_post_trim:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc_post_trim.zip", sample=SAMPLES, end=["1","2"])
    output:
        "results/qc/post_trim/multiqc.html",
        "results/qc/post_trim/multiqc_data/multiqc_general_stats.txt",
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
        extra="",  # Optional: extra parameters for multiqc
    threads: 1
    resources:
        runtime=15,
    log:
        "logs/multiqc/multiqc_post_trim.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "multiqc " 
        "--force "
        "-o {params.dir} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "2> {log}"