import pandas as pd
from snakemake.shell import shell
import csv


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Load sample information
sample_info = pd.read_csv("config/samples.csv")

# Load Snakemake variables
sample = snakemake.wildcards["ip_sample"]
control_available = snakemake.params["control"]
mode = snakemake.params["mode"]
outdir = snakemake.params["outdir"]
ip_bam = snakemake.input["bam"]
qvalue = snakemake.params["qvalue"]
extra = snakemake.params["extra"]
egs = snakemake.input["egs"]

# Check if control is available
if not control_available:
    control = ""
else:
    control_bam = snakemake.input["cbam"]
    control = f"--control {control_bam}"

if mode == "broad":
    broad_cutoff = str(snakemake.params["bc"])
    mode = f"--broad --broad-cutoff {broad_cutoff}"
else:
    mode = ""
 
# Retrieve effective genome size for sample
with open(egs, 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)  # Skip header row
    for row in csv_reader:
        sample_name = row[0]
        if sample_name == sample:
            effective_genome_size = int(row[1])
            break
 
# Run MACS2
shell(
    "macs2 callpeak "
    "--treatment {ip_bam} "
    "{control} "
    "--gsize {effective_genome_size} "
    "--qvalue {qvalue} "
    "--format BAMPE "
    "--keep-dup all "
    "{mode} "
    "--outdir {outdir} "
    "--name {sample} "
    "{extra} "
    "{log} "
)
