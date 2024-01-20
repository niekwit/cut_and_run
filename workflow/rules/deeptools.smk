if config["deduplication"]:
    use rule bam_index as bam_index2 with:
        input:
            "results/mapped/{bw_input_dir}/{sample}.bam",
        output:
            "results/mapped/{bw_input_dir}/{sample}.bam.bai",
        log:
            "logs/samtools_index/{sample}_{bw_input_dir}.log",


rule multiBigwigSummary:
    input:
        expand("results/bigwig/{bw_input_dir}/{sample}.bw", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
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
        bam=expand("results/mapped/{bw_input_dir}/{sample}.bam", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
        bai=expand("results/mapped/{bw_input_dir}/{sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
    output:
        hist="results/plots/qc/fragment_lengths.pdf",
        table="results/qc/fragment_lengths.tsv",
        raw="results/qc/fragment_lengths_raw.tsv",
    params:
        names= " ".join(SAMPLES),
        dir=BW_INPUT_DIR[0],
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
        "--plotTitle 'Fragment size of PE data ({params.dir})' "
        "> {log} 2>&1"


rule computeMatrix:
    input:
        bw=expand("results/bigwig/average_bw/{condition}.bw", condition=CONDITIONS),
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


