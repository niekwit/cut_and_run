# Downsampling of reads with seqtk (only for samples set in config/downsample.csv)
# Using the same seed for R1 and R2 keeps read pairs in sync
if DOWNSAMPLE:
    if PAIRED_END:

        rule downsample_pe:
            input:
                r1="reads/{sample}_R1_001.fastq.gz",
                r2="reads/{sample}_R2_001.fastq.gz",
            output:
                r1="results/downsampled/{sample}_R1_001.fastq.gz",
                r2="results/downsampled/{sample}_R2_001.fastq.gz",
            params:
                fraction=lambda wildcards: DOWNSAMPLE[wildcards.sample] / 100,
                seed=config["downsample_seed"],
            wildcard_constraints:
                sample="|".join(re.escape(s) for s in DOWNSAMPLE),
            log:
                "logs/downsample/{sample}.log",
            threads: config["resources"]["trim"]["cpu"]
            resources:
                runtime=config["resources"]["trim"]["time"],
            conda:
                "../envs/seqtk.yaml"
            shell:
                "(seqtk sample -s {params.seed} {input.r1} {params.fraction} | "
                "pigz -p {threads} > {output.r1}; "
                "seqtk sample -s {params.seed} {input.r2} {params.fraction} | "
                "pigz -p {threads} > {output.r2}) 2> {log}"

    else:

        rule downsample_se:
            input:
                "reads/{sample}.fastq.gz",
            output:
                "results/downsampled/{sample}.fastq.gz",
            params:
                fraction=lambda wildcards: DOWNSAMPLE[wildcards.sample] / 100,
                seed=config["downsample_seed"],
            wildcard_constraints:
                sample="|".join(re.escape(s) for s in DOWNSAMPLE),
            log:
                "logs/downsample/{sample}.log",
            threads: config["resources"]["trim"]["cpu"]
            resources:
                runtime=config["resources"]["trim"]["time"],
            conda:
                "../envs/seqtk.yaml"
            shell:
                "seqtk sample -s {params.seed} {input} {params.fraction} 2> {log} | "
                "pigz -p {threads} > {output}"
