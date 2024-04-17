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


rule plot_heatmap:
    input:
        mat="results/deeptools/matrix.gz",
    output:
        pdf="results/plots/heatmap.pdf",
        mat="results/deeptools/heatmap_matrix.gz",
    params:
        im = config["deeptools"]["plotHeatmap"]["interpolationMethod"],
        pt = config["deeptools"]["plotHeatmap"]["plotType"],
        cm = config["deeptools"]["plotHeatmap"]["colorMap"],
        a = config["deeptools"]["plotHeatmap"]["alpha"],
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/plotHeatmap.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotHeatmap "
        "--matrixFile {input.mat} "
        "--outFileNameMatrix {output.mat} "
        "--outFileName {output.pdf} "
        "--perGroup "
        "--colorMap {params.cm} "
        "--alpha {params.a} "
        "{params.extra} "
        "> {log} 2>&1"


if config["peak_calling"]["htseq_count"]["use_htseq_count"]:
    rule plot_peaks_volcano:
        input:
            xlsx="results/peaks/DESeq2/differential_peaks.xlsx",
            dir="results/peaks/DESeq2/",
        output:
            directory("results/plots/differential_peaks/volcano_plots"),
        params:
            fdr=config["peak_calling"]["htseq_count"]["DESeq2"]["alpha"],
            fc=config["peak_calling"]["htseq_count"]["DESeq2"]["fc"],
            extra="",
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"]
        log:
            "logs/plots/volcano_peaks.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/plot_peaks_volcano.R"