import pandas as pd
from snakemake.shell import shell
import csv
import statistics
import re
import os

def get_egs(sample):
    # Retrieve effective genome size for sample
    with open(egs, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header row
        for row in csv_reader:
            sample_name = row[0]
            if sample_name == sample:
                return int(row[1])
 
# Load sample information
_csv = pd.read_csv("config/samples.csv")

# Load Snakemake variables
control_available = snakemake.params["control"]
mode = snakemake.params["mode"]
tbams = snakemake.input["bams"]
qvalue = snakemake.params["qvalue"]
extra = snakemake.params["extra"]
egs = snakemake.input["egs"]
xls_out = snakemake.output["xls"]
logs = snakemake.log

# Set peak calling mode
if mode == "broad":
    broad_cutoff = str(snakemake.params["bc"])
    mode = f"--broad --broad-cutoff {broad_cutoff} "
else:
    mode = ""

# Run MACS2 for each condition and match controls
conditions = list(set(re.sub("_[\d]$", "", x) for x in _csv["sample"].tolist()))
for c in conditions:
    #log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    
    bams = [x for x in tbams if c in x]
    if not control_available:
        control = ""
    else:
        # Get matching control bam files for c
        cbams = snakemake.input["cbams"]
        control_names = _csv[_csv["sample"].str.contains(c)]["control"].tolist()
        control_bams = []
        for i in cbams:
            for j in control_names:
                if j in i:
                    control_bams.append(i)
        control_bams = list(set(control_bams))
        control = f"--control {' '.join(control_bams)}"

    # Get average effective genome size for all samples matching condition
    effective_genome_sizes = []
    samples = _csv[_csv["sample"].str.contains(c)]["sample"].tolist()
    for s in samples:
        effective_genome_sizes.append(get_egs(s))
    mean_egs = int(statistics.fmean(effective_genome_sizes))
    
    # Create outdir name
    outdir = f"{os.path.dirname(os.path.dirname(xls_out[0]))}/{c}"

    # Run MACS2
    shell(
        "macs2 callpeak "
        "--treatment {bams} "
        "{control} "
        "--gsize {mean_egs} "
        "--qvalue {qvalue} "
        "--format BAMPE "
        "--keep-dup all "
        "{mode}"
        "--outdir {outdir} "
        "--name {c} "
        "{extra} "
        #"{log} "
    )

