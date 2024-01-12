import pandas as pd
import os

def samples():
    """Returns sample wildcard values for Snakemake
    """
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"]
    
    # check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        r1= f"reads/{sample}_R1_001.fastq.gz"
        r2= f"reads/{sample}_R2_001.fastq.gz"
        if not os.path.isfile(r1):
            not_found.append(r1)
        if not os.path.isfile(r2):
            not_found.append(r2)
    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"ERROR: following files not found:\n{not_found}")
        
    return SAMPLES         


def cutadapt_args(config, param):
    """Returns cutadapt adapter or extra arguments as string read from config file
    """
    if param == "adapters":
        a_arg = config["cutadapt"]["a"]
        A_arg = config["cutadapt"]["A"]
        return f'-a "{a_arg}" -A "{A_arg}"'
    elif param == "extra":
        return f"--minimum-length {config['cutadapt']['min_length']} {config['cutadapt']['extra']}"


def bw_input_dir(config):
    """Input function for bigwig rule.
    Determines which bam files to use for bigwig generation: deduplicated or not.
    """
    if config["deduplication"]:
        return ["deduplicated"]
    else:
        return ["bl_removed"]
    

def fasta(config, resources):
    if config["remove_MT_seqs"]:
        return resources.nomt_fasta
    else:
        return resources.fasta


def bw_input_data(wildcards):
    """Returns named input files as dictionary for bigwig rule.
    Bigwig rule input can change depdending whether spike-in normalization is applied or not and mitochondrial sequences are omited from the analysis.
    """
    # Create base input dictionary (these input files are always required)
    dict = {
        "bam": "results/mapped/{wildcards.bw_input_dir}/{wildcards.sample}.bam".format(wildcards=wildcards),
        "bai": "results/mapped/{wildcards.bw_input_dir}/{wildcards.sample}.bam.bai".format(wildcards=wildcards),
        "multiqc": "results/qc/multiqc_data/multiqc_general_stats.txt",
    }
    # Add additional input files depending on config file
    if config["apply_spike_in"]:
        dict["sf"] = "results/scale_factors/scale_factors.csv"
    if config["remove_MT_seqs"]:
        dict["mgs"] = f"resources/{genome}_mt_genome_size.txt"
    
    return dict


