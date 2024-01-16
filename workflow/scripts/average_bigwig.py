from snakemake.shell import shell

# Load Snakemake variables
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads

all_bw = snakemake.input
condition = snakemake.wildcards["condition"]
out = snakemake.output["bw"]

# Get all samples in condition
bw = [x for x in all_bw if condition in x] # use input lambda function instead? (just shell command in rule)

# Create average bigwig file
shell(
    "bigwigAverage "
    "--bigwigs {bw} "
    "--outFileName {out} "
    "--numberOfProcessors {threads} "
    "{log}"
    )

