rule bigwig:
    input:
        unpack(bw_input_data),
    output:
        "results/bigwig/{bw_input_dir}/{sample}.bw",
    params:
        genome=config["genome"],
        binsize=config["bigwig"]["binsize"],
        norm=config["bigwig"]["normalisation"],
        extra=config["bigwig"]["extra"],
        apply_spike_in=config["apply_spike_in"],
        remove_MT_seqs=config["remove_MT_seqs"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/bigwig/{bw_input_dir}/{sample}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/bigwig.py"


rule average_bigwigs:
    input:
        expand("results/bigwig/{bw_input_dir}/{sample}.bw", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES)