if config["deduplication"]:
    rule remove_duplicates_bam:
        input:
            bams="results/mapped/bl_removed/{sample}.bam",
            bai="results/mapped/bl_removed/{sample}.bam.bai",
        output:
            bam="results/mapped/deduplicated/{sample}.bam",
            metrics="results/mapped/deduplicated/{sample}.metrics.txt",
        log:
            "logs/dedup_bam/{sample}.log",
        params:
            extra="--REMOVE_DUPLICATES true",
        threads:
            config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["cpu"],
            mem_mb=2048,
        wrapper:
            "v3.3.6/bio/picard/markduplicates"

    rule dedup_bam_index:
        input:
            "results/mapped/deduplicated/{sample}.bam",
        output:
            "results/mapped/deduplicated/{sample}.bam.bai",
        params:
            extra="",  # optional params string
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"],
        log:
            "logs/samtools_index/{sample}.log",
        wrapper:
            "v3.3.6/bio/samtools/index"