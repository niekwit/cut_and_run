"""
Calculates the effective genome size for each
sample based on the read length and genome.
Required for bigwig generation and peak calling.
"""

import pandas as pd
import re

# Load Snakemake variables
genome = snakemake.params["genome"]
MT_seqs_removed = snakemake.params["remove_MT_seqs"]
multiqc = snakemake.input["multiqc"]

# Default effective genome sizes
default_effective_genome_size = {
    "dm6": {
        "50": 125464728,
        "75": 127324632,
        "100": 129789873,
        "150": 129941135,
        "200": 132509163,
    },
    "hg19": {
        "50": 2685511504,
        "75": 2736124973,
        "100": 2776919808,
        "150": 2827437033,
        "200": 2855464000,
    },
    "hg38": {
        "50": 2701495761,
        "75": 2747877777,
        "100": 2805636331,
        "150": 2862010578,
        "200": 2887553303,
    },
    "mm37": {
        "50": 2304947926,
        "75": 2404646224,
        "100": 2462481010,
        "150": 2489384235,
        "200": 2513019276,
    },
    "mm38": {
        "50": 2308125349,
        "75": 2407883318,
        "100": 2467481108,
        "150": 2494787188,
        "200": 2520869189,
    },
}
# Load multiqc file
multiqc = pd.read_csv(multiqc, sep="\t")

# Data frame to store effective genome size per sample
df = pd.DataFrame(columns=["sample", "effective_genome_size"])

# Iterate over rows to get read length and effective genome size and add to df
for index, row in multiqc.iterrows():
    sample = row[0]
    sample = re.sub(r"_R[12]_001", "", sample)
    read_length = row[4]
    
    # If read length is not in dict, pick the closest read length
    if read_length not in default_effective_genome_size[genome]:
        read_length = min(default_effective_genome_size[genome], key=lambda x: abs(int(x) - int(read_length)))
    
    effective_genome_size = default_effective_genome_size[genome][read_length]
    
    # Remove MT genome size from effective genome size
    if MT_seqs_removed:
        with open(snakemake.input["mgs"]) as f:
            mt_genome_size = int(f.read().strip())
        effective_genome_size -= mt_genome_size
    
    # Add to df
    df.loc[index] = [sample, effective_genome_size]
    
# Save df to file
df.to_csv(snakemake.output["egs"], index=False)


