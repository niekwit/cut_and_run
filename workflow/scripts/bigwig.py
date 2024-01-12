from snakemake.shell import shell
import pandas as pd

# Input/output files
bam = snakemake.input["bam"]
multiqc = snakemake.input["multiqc"]
output = snakemake.output[0]

# Parameters
genome = snakemake.params["genome"]
binsize = snakemake.params["binsize"]
norm = snakemake.params["norm"]
extra = snakemake.params["extra"]
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
sample = snakemake.wildcards["sample"]

# Load scale factors for sample if spike in was applied
if snakemake.params["apply_spike_in"]:
    scale_factors = snakemake.input["sf"]
    scale_factors = pd.read_csv(scale_factors)
    scale_factor = scale_factors[scale_factors["sample"] == sample]["scale_factor"].values[0]
    
    # Set normalisation parameter to None 
    norm = "None"
    
else: # Apply no scale factor (i.e. scale factor = 1, set normalisation parameter in config to not None)
    if norm in ["None", None]:
        print("WARNING: No scale factor applied and no normalisation parameter set.\nThis will result in a bigwig file with raw read counts.")
    scale_factor = 1

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

# Read multiqc file to get read length
# WHY?: Some experiments have multiple read lengths (not ideal...), so we need to get the read length from the multiqc file for each file
with open(multiqc) as f:
    lines = f.readlines()
    line = [line.strip() for line in lines if sample in line][0]

read_length = line.split("\t")[4]

# Get effective genome size
effective_genome_size = default_effective_genome_size[genome][read_length]

# Correct effective genome size for removal of MT sequences, if applied
if snakemake.params["remove_MT_seqs"]:
    # Load MT genome size
    with open(snakemake.input["mgs"]) as f:
        mt_genome_size = int(f.read().strip())
    effective_genome_size -= mt_genome_size

# Run deeptools
shell(
    "bamCoverage "
    "--numberOfProcessors {threads} "
    "--effectiveGenomeSize {effective_genome_size} "
    "--bam {bam} "
    "--scaleFactor {scale_factor} "
    "--binSize {binsize} "
    "--normalizeUsing {norm} "
    "--extendReads "
    "--outFileName {output} "
    "{log} "
)


