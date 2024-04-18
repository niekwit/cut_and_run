rule multiBigwigSummary:
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
    output:
        "results/deeptools/scores_per_bin.npz",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/multiBigwigSummary.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins "
        "--bwfiles {input} "
        "--smartLabels "
        "--outFileName {output} "
        "--numberOfProcessors {threads} "
        "{params.extra} "
        "> {log} 2>&1"


rule PCA:
    input:
        "results/deeptools/scores_per_bin.npz",
    output:
        "results/deeptools/PCA.tab",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/PCA.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--outFileNameData {output} "
        "--transpose "
        "{params.extra} "
        "> {log} 2>&1"


rule BAM_fragment_sizes:
    input:
        bam=expand("results/mapped/{sample}.bl.bam", sample=SAMPLES),
        bai=expand("results/mapped/{sample}.bl.bam.bai", sample=SAMPLES),
    output:
        hist="results/plots/qc/fragment_lengths.pdf",
        table="results/qc/fragment_lengths.tsv",
        raw="results/qc/fragment_lengths_raw.tsv",
    params:
        names= " ".join(SAMPLES),
        max_len=config["bowtie2"]["max_length"],
    threads: config["resources"]["deeptools"]["cpu"],
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/deeptools/BAM_fragment_sizes.log",
    conda:
        "../envs/deeptools.yaml",
    shell:
        "bamPEFragmentSize "
        "--numberOfProcessors {threads} " 
        "--maxFragmentLength {params.max_len} "
        "--bamfiles {input.bam} "
        "--histogram {output.hist} " 
        "--table {output.table} "
        "--outRawFragmentLengths {output.raw} " 
        "--samplesLabel {params.names} " 
        "--plotTitle 'Fragment size of PE data' "
        "> {log} 2>&1"


rule computeMatrix:
    input:
        bw=expand("results/bigwig/average_bw/{conditions}.bw", conditions=CONDITIONS_NO_CONTROL),
        gtf=resources.gtf,
    output:
        mat="results/deeptools/matrix.gz",
    params:
        args=computematrix_args(),
    threads: config["resources"]["deeptools"]["cpu"] * 4 # Otherwise it will take very long
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/computeMatrix.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "computeMatrix "
        "{params.args} "
        "--numberOfProcessors {threads} "
        "--smartLabels "
        "--scoreFileName {input.bw} "
        "--outFileName {output.mat} "
        "> {log} 2>&1"


