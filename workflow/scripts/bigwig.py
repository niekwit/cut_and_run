from snakemake.shell import shell
import pandas as pd
import csv

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
multiqc = snakemake.input["multiqc"]
egs = snakemake.input["egs"]

# Load scale factors for sample if spike in was applied
if snakemake.params["apply_spike_in"]:
    scale_factors = snakemake.input["sf"]
    scale_factors = pd.read_csv(scale_factors)
    scale_factor = scale_factors[scale_factors["sample"] == sample]["scale_factor"].values[0]
    scaling = f"--scaleFactor {scale_factor}"
    # Set normalisation parameter to None 
    normalisation = ""
    
else: # Apply no scale factor (set normalisation parameter in config to not None)
    if norm in ["None", None]:
        print("WARNING: No scale factor applied and no normalisation parameter set.\nThis will result in a bigwig file with raw read counts.")
    scaling = ""
    normalisation = f"--normalizeUsing {norm}"

# Retrieve effective genome size for sample
with open(egs, 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)  # Skip header row
    for row in csv_reader:
        sample_name = row[0]
        if sample_name == sample:
            effective_genome_size = int(row[1])
            break

# Run deeptools
shell(
    "bamCoverage "
    "--numberOfProcessors {threads} "
    "--effectiveGenomeSize {effective_genome_size} "
    "--bam {bam} "
    "{scaling} "
    "--binSize {binsize} "
    "{normalisation} "
    #"--extendReads "
    "--outFileName {output} "
    "{extra} "
    "{log} "
)

