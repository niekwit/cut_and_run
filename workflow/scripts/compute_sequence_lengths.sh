#!/usr/bin/env bash

mkdir -p results/read_lengths/${snakemake_wildcards["bw_input_dir"]}

samtools view -@ {threads} ${snakemake_input[0]} | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort > ${snakemake_output[0]} 2> ${snakemake_log[0]}
