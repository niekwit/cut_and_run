from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads

bam_files = " ".join(snakemake.input["bam"])
histogram = snakemake.output["hist"]
table = snakemake.output["table"]

names = snakemake.params["names"]
dir = snakemake.params["dir"]
max_length = snakemake.params["max_len"]

shell(
    "bamPEFragmentSize "
    "--numberOfProcessors {threads} " 
    "--maxFragmentLength {max_length} "
    "--bamfiles {bam_files} "
    "--histogram {histogram} " 
    "--table {table} " 
    "--samplesLabel {names} " 
    "--plotTitle 'Fragment size of PE data ({dir})' "
    "{log}"
)



