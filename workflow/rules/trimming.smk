if PAIRED_END:
    if config["use_trim_galore"]:
        rule trim_galore_pe:
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
                "logs/trim_galore_pe/{sample}.log",
            threads: config["resources"]["trim"]["cpu"]
            resources:
                runtime=config["resources"]["trim"]["time"],
            wrapper:
                "v7.0.0/bio/trim_galore/pe"
    else:
        rule cutadapt_pe:
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
                "logs/cutadapt_pe/{sample}.log",
            threads: config["resources"]["trim"]["cpu"]
            resources:
                runtime=config["resources"]["trim"]["time"],
            wrapper:
                "v7.0.0/bio/cutadapt/pe"
else:
    if config["use_trim_galore"]:
        rule trim_galore_se:
            input:
                "reads/{sample}.fastq.gz",
            output:
                fasta="results/trimmed/{sample}.fq.gz",
                report="results/trimmed/report/{sample}_trimming_report.txt",
            params:
                extra="--illumina -q 20",
            log:
                "logs/trim_galore_se/{sample}.log",
            threads: config["resources"]["trim"]["cpu"]
            resources:
                runtime=config["resources"]["trim"]["time"],
            wrapper:
                "v7.0.0/bio/trim_galore/se"
    else:
        rule cutadapt_se:
            input:
                "reads/{sample}.fastq.gz",
            output:
                fastq="results/trimmed/{sample}.fq.gz",
                qc="results/trimmed/{sample}.qc.txt",
            params:
                adapters=cutadapt_args(config,"adapters"),
                extra=cutadapt_args(config,"extra"),
            log:
                "logs/cutadapt_se/{sample}.log",
            threads: 4  # set desired number of threads here
            wrapper:
                "v7.0.0/bio/cutadapt/se"

