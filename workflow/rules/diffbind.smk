if run_diffbind():

    rule diffbind:
        input:
            xls=expand(
                f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls",
                conditions=CONDITIONS,
            ),
            bam=expand("results/mapped/{sample}.bl.bam", sample=SAMPLES),
            bai=expand("results/mapped/{sample}.bl.bam.bai", sample=SAMPLES),
        output:
            dba=f"results/diffbind/{PEAK_MODE}/fdr{fdr}/dba.RData",
            pca=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/PCA.pdf",
            sc=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/SampleCorrelation.pdf",
            pp=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/ProfilePLot.pdf",
            sh=f"results/diffbind/{PEAK_MODE}/fdr{fdr}/samplesheet.csv",
        params:
            genome=genome,
        threads: config["resources"]["diffbind"]["cpu"]
        resources:
            runtime=config["resources"]["diffbind"]["time"],
            mem_mb=120000,
        conda:
            "../envs/diffbind.yaml"
        log:
            "logs/diffbind/{peak_mode}/fdr{fdr}/diffbind.log",
        script:
            "../scripts/diffbind.R"
