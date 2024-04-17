if config["use_trim_galore"]:
    rule trim_galore:
        input:
            ["reads/{sample}_R1_001.fastq.gz","reads/{sample}_R2_001.fastq.gz"],
        output:
            fasta_fwd="results/trimmed/{sample}_val_1.fq.gz",
            report_fwd="logs/trim_galore/{sample}_R1_trimming_report.txt",
            fasta_rev="results/trimmed/{sample}_val_2.fq.gz",
            report_rev="logs/trim_galore/{sample}_R2_trimming_report.txt",
        params:
            extra="--illumina -q 20",
        log:
            "logs/trim_galore/{sample}.log",
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        wrapper:
            f"{wrapper_version}/bio/trim_galore/pe"
else:
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
            f"{wrapper_version}/bio/cutadapt/pe"

