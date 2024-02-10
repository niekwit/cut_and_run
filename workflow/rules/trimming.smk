rule cutadapt:
    input:
        ["reads/{sample}_R1_001.fastq.gz","reads/{sample}_R2_001.fastq.gz"],
    output:
        fastq1="results/trimmed/{sample}_val_1.fq.gz",
        fastq2="results/trimmed/{sample}_val_2.fq.gz",
        qc="logs/cutadapt/{sample}.qc.txt",
    params:
        adapters=cutadapt_args(config,"adapters"),
        extra=cutadapt_args(config,"extra"),
    log:
        "logs/cutadapt/{sample}.log",
    threads: config["resources"]["trim"]["cpu"]
    resources:
        runtime=config["resources"]["trim"]["time"],
    wrapper:
        "v3.3.6/bio/cutadapt/pe"

