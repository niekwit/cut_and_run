if PAIRED_END:

    rule bowtie2_align_pe:
        input:
            idx=multiext(
                f"resources/bowtie2_index/{genome}_{resources.build}/index",
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
            temp("results/mapped/bowtie2_genome/{sample}.bam"),
        params:
            idx=lambda wc, input: input[0].replace(".1.bt2", ""),
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
            "../envs/mapping.yaml"
        # Mapping based on:
        # https://elifesciences.org/articles/21856 Henikoff Cut&Run paper
        # https://github.com/peteskene/py_bowtie_fastq_2_sam/blob/master/py_bowtie_fastq_2_sam.py
        shell:
            "bowtie2 "
            "--local "
            "--very-sensitive-local "
            "--soft-clipped-unmapped-tlen "
            "--dovetail "
            "--no-unal "
            "--no-mixed "
            "--no-discordant "
            "--phred33 "
            "-I {params.min_len} "
            "-X {params.max_len} "
            "--threads {threads} "
            "{params.extra} "
            "-x {params.idx} "
            "-1 {input.r1} "
            "-2 {input.r2} 2> {log} | "
            "samtools view "
            "--min-MQ {params.min_mq} "
            "{params.samtools_extra} "
            "-bhS > {output}"
            # read characters from one or both ends of the alignment might be trimmed to maximize the alignment score
            # Consider soft-clipped bases unmapped when calculating TLEN (observed Template LENgth, is SAM field)
            # Considers cases where the mate alignments dovetail as concordant
            # Suppress SAM records for reads that failed to align
            # Bowtie 2 runs a little faster , but will only consider alignment status of pairs per se, not individual mates.
            # Disables searching for discordant alignments (does not satisfy the paired-end constraints)
            # Input qualities are ASCII chars equal to the Phred quality plus 33 (data dependent?)
            # The minimum fragment length for valid paired-end alignments
            # The maximum fragment length for valid paired-end alignments
            # Convert SAM to BAM
            # Minimum mapping quality

else:

    rule bowtie2_align_se:
        input:
            idx=multiext(
                f"resources/bowtie2_index/{genome}_{resources.build}/index",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
            sample="results/trimmed/{sample}.fq.gz",
        output:
            temp("results/mapped/bowtie2_genome/{sample}.bam"),
        params:
            idx=lambda wc, input: input[0].replace(".1.bt2", ""),
            min_mq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["bowtie2"]["extra"],
            samtools_extra=config["bowtie2"]["samtools_extra"],
        threads: config["resources"]["mapping"]["cpu"]
        resources:
            runtime=config["resources"]["mapping"]["time"],
        log:
            "logs/bowtie2_align_genome/{sample}.log",
        conda:
            "../envs/mapping.yaml"
        shell:
            "bowtie2 "
            "--threads {threads} "
            "{params.extra} "
            "-x {params.idx} "
            "-U {input.sample} 2> {log} | "
            "samtools view "
            "--min-MQ {params.min_mq} "
            "{params.samtools_extra} "
            "-bhS > {output}"
            # Convert SAM to BAM
            # Minimum mapping quality


if config["spike_in"]["apply_spike_in"]:
    # TODO: single-end data support

    rule bowtie2_build_spike_in:
        input:
            ref=resources_spike_in.fasta,
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
        cache: True
        threads: config["resources"]["index"]["cpu"]
        resources:
            runtime=config["resources"]["index"]["time"],
        wrapper:
            "v7.0.0/bio/bowtie2/build"

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
            temp("results/mapped/bowtie2_spike_in/{sample}.bam"),
        params:
            idx="resources/bowtie2_index/spike_in/index",
            extra="",  # optional parameters
        threads: config["resources"]["mapping"]["cpu"]
        resources:
            runtime=config["resources"]["mapping"]["time"],
        log:
            "logs/bowtie2_align_spike_in/{sample}.log",
        conda:
            "../envs/mapping.yaml"
        shell:
            "bowtie2 "
            "--local "
            "--very-sensitive-local "
            "--no-unal "
            "--no-mixed "
            "--no-discordant "
            "--phred33 "
            "--threads {threads} "
            "-x {params.idx} "
            "{params.samtools_extra} "
            "-1 {input.r1} "
            "-2 {input.r2} 2> {log} | "
            "samtools view "
            "-bhS > {output} "
            # read characters from one or both ends of the alignment might be trimmed to maximize the alignment score
            # Suppress SAM records for reads that failed to align
            # Bowtie 2 runs a little faster , but will only consider alignment status of pairs per se, not individual mates.
            # Disables searching for discordant alignments (does not satisfy the paired-end constraints)
            # Input qualities are ASCII chars equal to the Phred quality plus 33 (data dependent?)
            # Convert SAM to BAM


rule bam_sort:
    input:
        "results/mapped/bowtie2_genome/{sample}.bam",
    output:
        temp("results/mapped/sorted/{sample}.bam"),
    params:
        extra="-m 4G",  # optional params string
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/bam_sort/{sample}.log",
    wrapper:
        "v7.0.0/bio/samtools/sort"


rule remove_blacklisted_regions:
    input:
        left="results/mapped/sorted/{sample}.bam",
        right=resources.ensembl_blacklist,
    output:
        "results/mapped/{sample}.bl.bam",
    params:
        extra="-v ",  # Only keeps regions in bam file that are not in bed file
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"],
    log:
        "logs/bedtools_bl/{sample}.log",
    wrapper:
        "v7.0.0/bio/bedtools/intersect"


rule bam_index:
    input:
        "results/mapped/{sample}.bl.bam",
    output:
        "results/mapped/{sample}.bl.bam.bai",
    params:
        extra="",  # Optional params string
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/samtools_index/{sample}.log",
    wrapper:
        "v7.0.0/bio/samtools/index"
