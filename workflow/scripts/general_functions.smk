import pandas as pd
import os
import re
from snakemake.logging import logger

def targets():
    """
    Returns file targets for rule all
    """
    # Base targets
    TARGETS = [
        "results/plots/qc/mapping_rates.pdf",
        "results/plots/qc/mapping_read_number.pdf",
        "results/plots/qc/PCA.pdf",
        "results/plots/qc/scree.pdf",
        "results/plots/qc/bam_fragment_lengths.pdf",
        "results/plots/heatmap.pdf",
        "results/deeptools/heatmap_matrix.gz",
    ]

    ### Add conditional targets
    # Peak calling: 
    # there is a lot of redundancy in the peak calling targets code but
    # it was very tricky to make a general function for this, optimise later?
    if config["peak_calling"]["macs2"]["use_macs2"]:
        if config["peak_calling"]["macs2"]["broad"]:
            TARGETS.extend([
                f"results/plots/macs2_broad/fdr{fdr}/peaks_distance_to_TSS.pdf",
                f"results/plots/macs2_broad/fdr{fdr}/peak_distributions.pdf",
            ]) 
            if control_available():
                TARGETS.extend([
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.broadPeak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.gappedPeak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                ])
            else:
                TARGETS.extend([
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.xls", ip_sample=IP_SAMPLES),
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.broadPeak", ip_sample=IP_SAMPLES),
                    expand(f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.gappedPeak", ip_sample=IP_SAMPLES),
                ])
        else:
            TARGETS.extend([
                f"results/plots/macs2_narrow/fdr{fdr}/peaks_distance_to_TSS.pdf",
                f"results/plots/macs2_narrow/fdr{fdr}/peak_distributions.pdf",
            ])
            if control_available():
                TARGETS.extend([
                    expand(f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                    expand(f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.narrowPeak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                    expand(f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_summits.bed", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
                ]) 
            else:
                TARGETS.extend([
                    expand(f"results/macs2_narrow/{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.xls", ip_sample=IP_SAMPLES),
                    expand(f"results/macs2_narrow/{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.narrowPeak", ip_sample=IP_SAMPLES),
                    expand(f"results/macs2_narrow/{fdr}/{{ip_sample}}/{{ip_sample}}_summits.bed", ip_sample=IP_SAMPLES),
                ])
        TARGETS.extend([
            expand(f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.bed", conditions= CONDITIONS_NO_CONTROL),
            expand(f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_annotated.peaks.txt", conditions= CONDITIONS_NO_CONTROL),
        ])
    elif config["peak_calling"]["htseq_count"]["use_htseq_count"]:
        TARGETS.extend([
            "results/htseq_count/DESeq2/differential_peaks.xlsx",
            "results/htseq_count/DESeq2/dds.RData"
                        ])
    
    return TARGETS


def check_config():
    """
    Checks if config file values are valid entries
    REPLACE WITH SCHEMA VALIDATION
    """
    # List to store parameters with wrong value
    wrong = []
    
    # Check boolean values
    check = {
        "apply_spike_in": config["spike-in"]["apply_spike_in"],
        "remove_MT_seqs": config["remove_MT_seqs"],
    }
    for key, value in check.items():
        if not isinstance(value, bool):
            wrong.append(key)
    
    # Check numerical values
    check = {
        "cutadapt_min_length": config["cutadapt"]["min_length"],
        "bowtie2_min_length": config["bowtie2"]["min_length"],
        "bowtie2_max_length": config["bowtie2"]["max_length"],
        "bowtie2_min_mapq": config["bowtie2"]["MAPQ_cutoff"],
        "bigwig_binsize": config["deeptools"]["bigwig"]["binsize"],
    }
    for key, value in check.items():
        if not isinstance(value, int):
            wrong.append(key)

    # Raise error if any invalid values were found
    if len(wrong) != 0:
        wrong = "\n".join(wrong)
        raise ValueError(f"ERROR: following config parameters have invalid values:\n{wrong}")


def samples():
    """
    Checks sample names/files and returns sample wildcard values for Snakemake
    """
    SAMPLES = csv["sample"].tolist() # this gets just IP samples
    try:
        # This adds control samples
        # Some IP samples might share the same control sample 
        # so unique control samples are extracted
        SAMPLES.extend(list(set(csv["control"]))) 
    except KeyError:
        # No control samples
        pass
    
    # Check if sample names contain any characters that are not alphanumeric or underscore
    illegal = []
    for sample in SAMPLES:
        if not re.match("^[a-zA-Z0-9_]*$", sample):
            illegal.append(sample)
    if len(illegal) != 0:
        illegal = "\n".join(illegal)
        raise ValueError(f"Following samples contain illegal characters:\n{illegal}")

    # Check if each sample name ends with _[0-9]
    wrong = []
    for sample in SAMPLES:
        if not re.match(".*_[0-9]$", sample):
            wrong.append(sample)
    if len(wrong) != 0:
        wrong = "\n".join(wrong)
        raise ValueError(f"Following samples do not end with _[0-9]:\n{wrong}")

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
        raise ValueError(f"Following files not found:\n{not_found}")

    return SAMPLES         


def control_available():
    """
    Check if control samples are available in csv
    """
    try:
        test = csv["control"]
        return True
    except KeyError:
        return False


def cutadapt_args(config, param):
    """
    Returns cutadapt adapter or extra arguments as string read from config file
    """
    if param == "adapters":
        a_arg = config["cutadapt"]["a"]
        A_arg = config["cutadapt"]["A"]
        return f'-a "{a_arg}" -A "{A_arg}"'
    elif param == "extra":
        return f"--minimum-length {config['cutadapt']['min_length']} {config['cutadapt']['extra']}"


def conditions(include_controls=False):
    """
    Returns condition wildcard values 
    """
    # Get unique conditions from sample names
    conditions = list(set(re.sub("_[0-9]$", "", x) for x in csv["sample"].tolist()))

    # The same but for control samples and add to conditions
    if include_controls:
        if control_available():
            conditions.extend(list(set(re.sub("_[0-9]$", "", x) for x in csv["control"].tolist())))

    return conditions


def ip_samples():
    """Returns lists of paired IP/control samples for peak calling
    """
    if config["peak_calling"]["macs2"]["use_macs2"] or config["peak_calling"]["htseq_count"]["use_htseq_count"]:
        ip_samples = csv["sample"].tolist()
        
        if control_available():
            input_samples = csv["control"].tolist()
        else:
            input_samples = []
            
            logger.info("WARNING: No input/IgG/control samples applied for peak calling...")
            logger.info("Peak calling will continue without control samples\n...")
    else:
        ip_samples = []
        input_samples = []
        logger.info("WARNING: Skipping peak calling (no peak calling method selected)...")

    return ip_samples, input_samples


def fasta(config, resources):
    if config["remove_MT_seqs"]:
        return resources.nomt_fasta
    else:
        return resources.fasta


def bw_input(wildcards):
    """
    Returns named input files as dictionary for bigwig rule.
    Bigwig rule input can change depdending whether spike-in normalization is applied or not and mitochondrial sequences are omited from the analysis.
    """
    # Create base input dictionary (these input files are always required)
    _dict = {
        "bam": "results/mapped/{wildcards.sample}.bl.bam".format(wildcards=wildcards),
        "bai": "results/mapped/{wildcards.sample}.bl.bam.bai".format(wildcards=wildcards),
        "multiqc": "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",#"results/qc/pre_trim/multiqc_data/multiqc_general_stats.txt",
        "egs": "results/effective_genome_sizes/effective_genome_sizes.csv",
    }
    # Add additional input files depending on config file
    if config["spike-in"]["apply_spike_in"]:
        _dict["sf"] = "results/scale_factors/scale_factors.csv"
    if config["remove_MT_seqs"]:
        _dict["mgs"] = f"resources/{genome}_mt_genome_size.txt"
    
    return _dict


def computematrix_args():
    """
    Returns computeMatrix arguments as string based on config file
    """
    # Add mode argument
    mode = config["deeptools"]["matrix"]["mode"]
    if mode == "scale-regions".lower():
        rbl = config["deeptools"]["matrix"]["regionBodyLength"]
        args = f"scale-regions --regionBodyLength {rbl} "
    elif mode == "reference-point".lower():
        rp = config["deeptools"]["matrix"]["referencePoint"]
        args = f"reference-point --referencePoint {rp} "
    else:
        raise ValueError(f"ERROR: deeptools matrix mode {mode} not supported")
    
    
    # Add common arguments
    b = config["deeptools"]["matrix"]["upstream"]
    a = config["deeptools"]["matrix"]["downstream"]
    bs = config["deeptools"]["matrix"]["binSize"]
    atb = config["deeptools"]["matrix"]["averageTypeBins"]
        
    args = f"{args} --upstream {b} --downstream {a} --binSize {bs} --averageTypeBins {atb} "

    # Add region argument
    r   = config["deeptools"]["matrix"]["regionsFileName"]
    no_whole_genome = config["deeptools"]["matrix"]["no_whole_genome"]

    if no_whole_genome and r:
        args = f"{args} --regionsFileName {r} "
    elif not no_whole_genome and r:
        args = f"{args} --regionsFileName {resources.gtf} {r} "
    else: 
        args = f"{args} --regionsFileName {resources.gtf} "
    
    return args


def macs2_mode():
    """
    Returns macs2 peak calling mode as string based on config file
    """
    broad = config["peak_calling"]["macs2"]["broad"]
    if not isinstance(broad, bool):
        raise ValueError(f"ERROR: config.peak_calling.macs2.broad must be True or False")
    
    if broad:
        return "broad"
    else:
        return "narrow"


def run_diffbind():
    """
    Returns True if diffbind rule should be run
    """
    if config["peak_calling"]["macs2"]["use_macs2"]:
        references = list(set(csv["reference"]))
        if len(references) > 1 and "yes" in references:
            return True
        else:
            False
    else:
        False


def diffbind_input(wildcards):
    """
    Returns named input files as dictionary for diffbind rule.
    """
    mode = peak_mode()

    # Base input
    _dict = {
            "ip_bam": expand("results/mapped/{wildcards.ip_sample}.bam".format(wildcards=wildcards)),
            "bai": expand("results/mapped/{wildcards.ip_sample}.bam.bai".format(wildcards=wildcards)),
            "xls": expand(f"results/{mode}/fdr{fdr}/{{wildcard.ip_sample}}/{{wildcard.ip_sample}}_vs_{{wildcard.control_sample}}_peaks.xls".format(wildcards=wildcards)),
        }
    # Add control bam files if available
    if control_available():
        _dict["control_bam"] = expand("results/mapped/{wildcards.control_sample}.bam".format(wildcards=wildcards)),
        _dict["control_bai"] = expand("results/mapped/{wildcards.control_sample}.bam.bai".format(wildcards=wildcards))
    
    return _dict


def mitochondrial_genome_name():
    """
    Returns mitochondrial genome name as string based species.
    Needed for removing mitochondrial sequences from fasta file.
    """
    if "hg" in genome or "mm" in genome:
        return "MT"
    elif "dm" in genome:
        return "mitochondrion_genome"


def calculate_effective_genome_sizes_input():
    """
    Returns input files for calculate_effective_genome_sizes rule
    """
    _dict = {
            "multiqc": "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt", #"results/qc/pre_trim/multiqc_data/multiqc_general_stats.txt",
        }

    if config["remove_MT_seqs"]:
        _dict["mgs"] = f"resources/{genome}_mt_genome_size.txt"

    return _dict


def peak_mode():
    """
    Returns MACS2 peak calling mode as string based on config file
    """
    if config["peak_calling"]["macs2"]["broad"]:
        return "macs2_broad"
    else:
        return "macs2_narrow"


def peak_fdr(type_):
    if type_ == "macs2_narrow":
        return config["peak_calling"]["macs2"]["qvalue"]
    elif type_ == "macs2_broad":
        return config["peak_calling"]["macs2"]["broad_cutoff"]