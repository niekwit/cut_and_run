rule bowtie2_build:
    input:
        ref=fasta(config, resources),
    output:
        multiext(
            f"resources/bowtie2_index/{genome}/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build_genome/build.log",
    params:
        extra="",  # optional parameters
    threads: config["resources"]["index"]["cpu"]
    resources:
        runtime=config["resources"]["index"]["time"],
    wrapper:
        "v3.3.6/bio/bowtie2/build"


rule bowtie2_align:
    input:
        idx=multiext(
            f"resources/bowtie2_index/{genome}/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        r1="results/trimmed/{sample}_val_1.fq.gz",
        r2="results/trimmed/{sample}_val_2.fq.gz",
    output:
        "results/mapped/bowtie2_genome/{sample}.bam",
    params:
        idx=f"resources/bowtie2_index/{genome}/index",
        min_len=config["bowtie2"]["min_length"],
        max_len=config["bowtie2"]["max_length"],
        min_mq=config["bowtie2"]["MAPQ_cutoff"],
        extra=config["bowtie2"]["extra"],
        samtools_extra=config["bowtie2"]["samtools_extra"],
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"],
    log:
        "logs/bowtie2_align_genome/{sample}.log",
    conda:
        "../envs/mapping.yaml",
    # Mapping based on:
    # https://elifesciences.org/articles/21856 Henikoff Cut&Run paper
    # https://github.com/peteskene/py_bowtie_fastq_2_sam/blob/master/py_bowtie_fastq_2_sam.py
    # and Marta's email
    shell:
        "bowtie2 "
        #"--no-hd " # Suppress SAM header lines (starting with @).
        "--local " # read characters from one or both ends of the alignment might be trimmed to maximize the alignment score
        "--very-sensitive-local "
        "--soft-clipped-unmapped-tlen " # Consider soft-clipped bases unmapped when calculating TLEN (observed Template LENgth, is SAM field)
        "--dovetail " # Considers cases where the mate alignments dovetail as concordant
        "--no-unal " # Suppress SAM records for reads that failed to align
        "--no-mixed " # Bowtie 2 runs a little faster , but will only consider alignment status of pairs per se, not individual mates.
        "--no-discordant " # Disables searching for discordant alignments (does not satisfy the paired-end constraints)
        "--phred33 " # Input qualities are ASCII chars equal to the Phred quality plus 33 (data dependent?)
        "-I {params.min_len} " # The minimum fragment length for valid paired-end alignments
        "-X {params.max_len} " # The maximum fragment length for valid paired-end alignments
        "--threads {threads} "
        "{params.extra} "
        "-x {params.idx} "
        "-1 {input.r1} "
        "-2 {input.r2} 2> {log} | "
        "samtools view " # Convert SAM to BAM
        "--min-MQ {params.min_mq} " # Minimum mapping quality
        "{params.extra} "
        "-bhS > {output} "


if config["spike-in"]["apply_spike_in"]:
    rule bowtie2_build_spike_in:
        input:
            ref=resources_spike_in.fasta
        output:
            multiext(
                "resources/bowtie2_index/spike_in/index",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        log:
            "logs/bowtie2_build_spike_in/build.log",
        params:
            extra="",  # optional parameters
        threads: config["resources"]["index"]["cpu"]
        resources:
            runtime=config["resources"]["index"]["time"],
        wrapper:
            "v3.3.6/bio/bowtie2/build"


    rule bowtie2_align_spike_in:
        input:
            idx=multiext(
                "resources/bowtie2_index/spike_in/index",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
            r1="results/trimmed/{sample}_val_1.fq.gz",
            r2="results/trimmed/{sample}_val_2.fq.gz",
        output:
            "results/mapped/bowtie2_spike_in/{sample}.bam",
        params:
            idx="resources/bowtie2_index/spike_in/index",
            extra="",  # optional parameters
        threads: config["resources"]["mapping"]["cpu"]
        resources:
            runtime=config["resources"]["mapping"]["time"],
        log:
            "logs/bowtie2_align_spike_in/{sample}.log",
        conda:
            "../envs/mapping.yaml",
        shell:
            "bowtie2 "
            #" --no-hd " # Suppress SAM header lines (starting with @).
            "--local " # read characters from one or both ends of the alignment might be trimmed to maximize the alignment score
            "--very-sensitive-local "
            "--no-unal " # Suppress SAM records for reads that failed to align
            "--no-mixed " # Bowtie 2 runs a little faster , but will only consider alignment status of pairs per se, not individual mates.
            "--no-discordant " # Disables searching for discordant alignments (does not satisfy the paired-end constraints)
            "--phred33 " # Input qualities are ASCII chars equal to the Phred quality plus 33 (data dependent?)
            "--threads {threads} "
            "-x {params.idx} "
            "{params.extra} "
            "-1 {input.r1} "
            "-2 {input.r2} 2> {log} | "
            "samtools view " # Convert SAM to BAM
            "-bhS > {output} "


rule bam_sort:
    input:
        "results/mapped/bowtie2_genome/{sample}.bam",
    output:
        "results/mapped/sorted/{sample}.bam",
    params:
        extra="-m 4G",  # optional params string
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/bam_sort/{sample}.log",
    wrapper:
        "v3.3.6/bio/samtools/sort"


rule remove_blacklisted_regions:
    input:
        left="results/mapped/sorted/{sample}.bam",
        right=resources.blacklist,
    output:
        "results/mapped/bl_removed/{sample}.bam",
    params:
        extra="-v ", # only keeps regions in bam file that are not in bed file
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"],
    log:
        "logs/bedtools_bl/{sample}.log",
    wrapper:
        "v3.3.6/bio/bedtools/intersect"


rule bam_index:
    input:
        "results/mapped/bl_removed/{sample}.bam",
    output:
        "results/mapped/bl_removed/{sample}.bam.bai",
    params:
        extra="",  # optional params string
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/samtools_index/{sample}.log",
    wrapper:
        "v3.3.6/bio/samtools/index"


