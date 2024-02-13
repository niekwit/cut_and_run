rule bigwig:
    input:
        unpack(bw_input),
    output:
        "results/bigwig/{bw_input_dir}/{sample}.bw",
    params:
        genome=config["genome"],
        binsize=config["deeptools"]["bigwig"]["binsize"],
        norm=config["deeptools"]["bigwig"]["normalisation"],
        extra=config["deeptools"]["bigwig"]["extra"],
        apply_spike_in=config["spike-in"]["apply_spike_in"],
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


rule average_wig:
    input:
        expand("results/bigwig/{bw_input_dir}/{sample}.bw", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
    output:
        wig=temp("results/bigwig/average_bw/{bw_input_dir}/{condition}.wig"),
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wiggletools/{bw_input_dir}/wig_average_{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig="results/bigwig/average_bw/{bw_input_dir}/{condition}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig/{bw_input_dir}/average_bw/{condition}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wigToBigWig/{bw_input_dir}/{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"

'''
rule average_bigwigs:
    input:
        expand("results/bigwig/{bw_input_dir}/{sample}.bw", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
    output:
        bw="results/bigwig/average_bw/{condition}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/bw_average_{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_bigwig.py"
'''

