rule plot_mapping_rates:
    input:
        bam=expand("results/mapped/bowtie2_genome/{sample}.bam", sample=SAMPLES),
        log=expand("logs/bowtie2_align_genome/{sample}.log", sample=SAMPLES),
    output:
        rates="results/plots/qc/mapping_rates.pdf",
        counts="results/plots/qc/mapping_read_number.pdf",
    threads: config["resources"]["plotting"]["cpu"],
    resources:
        runtime=config["resources"]["plotting"]["time"],
    log:
        "logs/plots/mapping_rates.log",
    conda:
        "../envs/R.yaml",
    script:
        "../scripts/plot_mapping_rates.R"


rule plot_rsequence_lengths:
    input:
        txt="results/read_lengths/{sample}.txt",
    output:
        pdf="results/read_lengths/{sample}.pdf",
    threads: config["resources"]["samtools"]["cpu"],
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/plots/read_lengths/{sample}.log",
    conda:
        "../envs/R.yaml",
    script:
        "../scripts/plot_sequence_lengths.R"


rule plot_PCA:
    input:
        "results/deeptools/PCA.tab",
    output:
        pca="results/plots/qc/PCA.pdf",
        scree="results/plots/qc/scree.pdf",
    params:
        extra=""
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/plotPCA.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"



