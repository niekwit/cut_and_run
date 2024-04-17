rule annotate_peaks:
    input:
        bed="results/peaks/overlapping_peaks/{bg_sample}.extended.bed",
        gtf=resources.gtf,
    output:
        txt="results/peaks/overlapping_peaks/{bg_sample}.annotated.txt",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotate_peaks/{bg_sample}.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/annotate_peaks.R"