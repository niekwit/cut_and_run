rule compute_sequence_lengths:
    input:
        "results/mapped/{bw_input_dir}/{sample}.bam",
    output:
        "results/sequence_lengths/{bw_input_dir}/{sample}.txt",
    threads: config["resources"]["samtools"]["cpu"],
    resources:
        runtime=config["resources"]["samtools"]["time"],
    log:
        "logs/sequence_lengths/{sample}_{bw_input_dir}_compute_sequence_lengths.log",
    conda:
        "../envs/sequence_lengths.yaml",
    script:
        "scripts/compute_sequence_lengths.sh"


