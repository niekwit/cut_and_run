if config["spike-in"]["apply_spike_in"]:
    rule calculate_scale_factors:
        input:
            bam=expand("results/mapped/{bw_input_dir}/{sample}.bam", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
            si_bam=expand("results/mapped/bowtie2_spike_in/{sample}.bam", sample=SAMPLES),
        output:
            "results/scale_factors/scale_factors.csv",
        params: 
            extra="",
        threads: config["resources"]["samtools"]["cpu"]
        resources:
            runtime=config["resources"]["samtools"]["time"]
        log:
            "logs/scale_factors/calculate_scale_factors.log"
        conda:
            "../envs/deeptools.yaml"
        script:
            "../scripts/calculate_scale_factors.py"

