rule bigwig:
    input:
        unpack(bw_input),
    output:
        "results/bigwig/{sample}.bw",
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
        "logs/bigwig/{sample}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/bigwig.py"


rule average_wig:
    input:
        expand("results/bigwig/{sample}.bw",  sample=SAMPLES),
    output:
        wig=temp("results/bigwig/average_bw/{condition}.wig"),
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wiggletools/wig_average_{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig="results/bigwig/average_bw/{condition}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig/average_bw/{condition}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wigToBigWig/{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output} > {log} 2>&1"
