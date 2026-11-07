import pandas as pd
import os, re
import glob
import datetime
from scripts.resources import GenomeResources
from snakemake.utils import min_version, validate
from snakemake.logging import logger


def targets():
    """
    Returns file targets for rule all
    """
    # Base targets
    TARGETS = [
        # "results/plots/qc/mapping_rates.pdf",
        # "results/plots/qc/mapping_read_number.pdf",
        "results/plots/qc/PCA.pdf",
        "results/plots/qc/scree.pdf",
        "results/plots/qc/bam_fragment_lengths.pdf",
        "results/plots/heatmap.pdf",
        "results/deeptools/heatmap_matrix.gz",
    ]

    ### Add conditional targets
    # Peak calling (per peak calling mode and FDR cutoff)
    if config["peak_calling"]["macs2"]["use_macs2"]:
        for peak_mode in PEAK_MODES:
            if peak_mode == "macs2_narrow":
                suffixes = ["peaks.xls", "peaks.narrowPeak", "summits.bed"]
            else:
                suffixes = ["peaks.xls", "peaks.broadPeak", "peaks.gappedPeak"]
            for fdr in FDRS[peak_mode]:
                for suffix in suffixes:
                    TARGETS.extend(macs2_sample_files(peak_mode, fdr, suffix))
                TARGETS.extend(
                    [
                        f"results/plots/{peak_mode}/fdr{fdr}/peaks_distance_to_TSS.pdf",
                        f"results/plots/{peak_mode}/fdr{fdr}/peak_distributions.pdf",
                    ]
                )
                TARGETS.extend(
                    expand(
                        f"results/{peak_mode}/fdr{fdr}/consensus_peaks/{{condition}}.annotated.peaks.txt",
                        condition=mode_conditions(peak_mode),
                    )
                )
    elif config["peak_calling"]["htseq_count"]["use_htseq_count"]:
        TARGETS.extend(
            [
                "results/htseq_count/DESeq2/differential_peaks.xlsx",
                "results/htseq_count/DESeq2/dds.RData",
            ]
        )

    return TARGETS


def samples(PAIRED_END):
    """
    Checks sample names/files and returns sample wildcard values for Snakemake
    """
    SAMPLES = csv["sample"].tolist()  # this gets just IP samples
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
    if PAIRED_END:
        for sample in SAMPLES:
            r1 = f"reads/{sample}_R1_001.fastq.gz"
            r2 = f"reads/{sample}_R2_001.fastq.gz"
            if not os.path.isfile(r1):
                not_found.append(r1)
            if not os.path.isfile(r2):
                not_found.append(r2)
    else:
        for sample in SAMPLES:
            r = f"reads/{sample}.fastq.gz"
            if not os.path.isfile(r):
                not_found.append(r)

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


def downsample_percentages():
    """
    Reads optional config/downsample.csv and returns a dictionary with the
    percentage of reads to keep for each sample that should be downsampled.
    Samples not in the file or set to 100 are not downsampled.
    """
    file = "config/downsample.csv"
    if not os.path.isfile(file):
        return {}

    df = pd.read_csv(file)
    for col in ["sample", "percentage"]:
        if col not in df.columns:
            raise ValueError(f"{file} must contain a '{col}' column")

    # Check if all samples are IP or control samples
    unknown = [s for s in df["sample"] if s not in SAMPLES]
    if len(unknown) != 0:
        unknown = "\n".join(unknown)
        raise ValueError(
            f"Following samples in {file} are not in samples.csv:\n{unknown}"
        )

    # Check for duplicate samples
    duplicated = df.loc[df["sample"].duplicated(), "sample"].tolist()
    if len(duplicated) != 0:
        duplicated = "\n".join(duplicated)
        raise ValueError(
            f"Following samples are listed more than once in {file}:\n{duplicated}"
        )

    # Check if percentages are numbers between 0 (exclusive) and 100 (inclusive)
    percentages = pd.to_numeric(df["percentage"], errors="coerce")
    wrong = df.loc[
        percentages.isna() | (percentages <= 0) | (percentages > 100), "sample"
    ].tolist()
    if len(wrong) != 0:
        wrong = "\n".join(wrong)
        raise ValueError(
            f"Following samples in {file} do not have a percentage between 0 and 100:\n{wrong}"
        )

    return {s: float(p) for s, p in zip(df["sample"], percentages) if p < 100}


def trim_input(wildcards):
    """
    Returns input fastq file(s) for trimming rules: downsampled reads if the
    sample is set in config/downsample.csv, otherwise the raw reads.
    """
    if wildcards.sample in DOWNSAMPLE:
        dir_ = "results/downsampled"
    else:
        dir_ = "reads"

    if PAIRED_END:
        return [
            f"{dir_}/{wildcards.sample}_R1_001.fastq.gz",
            f"{dir_}/{wildcards.sample}_R2_001.fastq.gz",
        ]
    return f"{dir_}/{wildcards.sample}.fastq.gz"


def regex(values):
    """
    Returns a regular expression that only matches the given values exactly
    (for wildcard constraints)
    """
    return "|".join(re.escape(str(v)) for v in values)


def mode_samples(peak_mode):
    """
    Returns the IP samples and their matching control samples (empty list
    without controls) that use the given MACS2 peak calling mode
    """
    ips = [s for s in IP_SAMPLES if SAMPLE_PEAK_MODES[s] == peak_mode]
    controls = [
        c
        for s, c in zip(IP_SAMPLES, CONTROL_SAMPLES)
        if SAMPLE_PEAK_MODES[s] == peak_mode
    ]
    return ips, controls


def mode_conditions(peak_mode):
    """
    Returns the consensus peak conditions of the IP samples that use the given
    MACS2 peak calling mode
    """
    ips, _ = mode_samples(peak_mode)
    return sorted(set(SAMPLE_CONDITIONS[s] for s in ips))


def macs2_sample_files(peak_mode, fdr, suffix):
    """
    Returns per-sample MACS2 output files with the given suffix
    (e.g. peaks.narrowPeak) for the given peak calling mode and FDR cutoff
    """
    ips, controls = mode_samples(peak_mode)
    if control_available():
        return expand(
            f"results/{peak_mode}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_{suffix}",
            zip,
            ip_sample=ips,
            control_sample=controls,
        )
    return expand(
        f"results/{peak_mode}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_{suffix}",
        ip_sample=ips,
    )


def macs2_sample_peaks(peak_mode, fdr):
    """
    Returns per-sample MACS2 peak files for the given peak calling mode
    and FDR cutoff
    """
    ext = "narrowPeak" if peak_mode == "macs2_narrow" else "broadPeak"
    return macs2_sample_files(peak_mode, fdr, f"peaks.{ext}")


def sample_peak_modes():
    """
    Returns a dictionary with the MACS2 peak calling mode of each IP sample
    (macs2_narrow or macs2_broad), read from the optional macs2 column
    (narrow/broad) of samples.csv. Without this column, the
    peak_calling:macs2:broad setting in config.yaml is used for all samples.
    """
    df = pd.read_csv("config/samples.csv", dtype=str, keep_default_na=False)

    if "macs2" not in df.columns:
        default = "broad" if config["peak_calling"]["macs2"]["broad"] else "narrow"
        return {s: f"macs2_{default}" for s in df["sample"]}

    modes = {}
    wrong = []
    for row in df.itertuples():
        mode = row.macs2.strip().lower()
        if mode not in ["narrow", "broad"]:
            wrong.append(f"{row.sample}: '{row.macs2}'")
        modes[row.sample] = f"macs2_{mode}"
    if wrong:
        wrong = "\n".join(wrong)
        raise ValueError(
            f"macs2 column in samples.csv must be narrow or broad:\n{wrong}"
        )

    # All samples of a factor must use the same mode
    df["mode"] = [modes[s] for s in df["sample"]]
    mixed = [
        factor for factor, group in df.groupby("factor") if group["mode"].nunique() > 1
    ]
    if mixed:
        raise ValueError(
            "Following factors have both narrow and broad in the macs2 column "
            f"of samples.csv: {', '.join(mixed)}"
        )

    return modes


def sample_conditions():
    """
    Returns a dictionary with the condition of each IP sample,
    i.e. the genotype, treatment and factor columns of samples.csv joined by "_"
    """
    # Read as strings so that values such as "None" are kept as is
    df = pd.read_csv("config/samples.csv", dtype=str, keep_default_na=False)

    columns = ["genotype", "treatment", "factor"]
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"Following columns are missing in samples.csv: {', '.join(missing)}"
        )

    conditions = {}
    for row in df.itertuples():
        values = [getattr(row, c).strip() for c in columns]
        if "" in values:
            raise ValueError(
                f"Sample {row.sample} has an empty genotype, treatment or factor "
                "in samples.csv"
            )
        conditions[row.sample] = "_".join(values)

    # Conditions are used in file names
    illegal = sorted(
        set(c for c in conditions.values() if not re.match(r"^[a-zA-Z0-9_.-]+$", c))
    )
    if illegal:
        illegal = "\n".join(illegal)
        raise ValueError(
            "Following genotype_treatment_factor conditions in samples.csv contain "
            f"illegal characters:\n{illegal}"
        )

    return conditions


def consensus_beds(wildcards):
    """
    Returns the MACS2 peak files of all IP samples that belong to the
    condition (exact match on genotype_treatment_factor) and peak calling mode
    of the consensus peaks rule
    """
    ips, _ = mode_samples(wildcards.peak_mode)
    beds = macs2_sample_peaks(wildcards.peak_mode, wildcards.fdr)
    return [
        bed
        for ip_sample, bed in zip(ips, beds)
        if SAMPLE_CONDITIONS[ip_sample] == wildcards.condition
    ]


def annotation_plot_beds(wildcards):
    """
    Returns the consensus peak files of all conditions of the peak calling
    mode and FDR cutoff of the peak annotation plots rule
    """
    return expand(
        f"results/{wildcards.peak_mode}/fdr{wildcards.fdr}/consensus_peaks/{{condition}}.bed",
        condition=mode_conditions(wildcards.peak_mode),
    )


def cutadapt_args(config, param):
    """
    Returns cutadapt adapter or extra arguments as string read from config file
    """
    if param == "adapters":
        args = ""
        args_list = ["a", "A", "g", "G", "b", "B"]
        for letter in args_list:
            seqs = config["cutadapt"][letter]
            if len(seqs) > 0:
                for s in seqs:
                    args = f'{args} -{letter} "{s}" '

        # Check if any adapter sequences were found
        if len(args) == 0:
            raise ValueError(
                "No cutadapt adapter sequences found in config file. "
                "Please add at least one adapter sequence."
            )
        return args

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
            conditions.extend(
                list(set(re.sub("_[0-9]$", "", x) for x in csv["control"].tolist()))
            )

    return conditions


def ip_samples():
    """Returns lists of paired IP/control samples for peak calling"""
    if (
        config["peak_calling"]["macs2"]["use_macs2"]
        or config["peak_calling"]["htseq_count"]["use_htseq_count"]
    ):
        ip_samples = csv["sample"].tolist()

        if control_available():
            input_samples = csv["control"].tolist()
        else:
            input_samples = []

            logger.info(
                "WARNING: No input/IgG/control samples applied for peak calling..."
            )
            logger.info("Peak calling will continue without control samples\n...")
    else:
        ip_samples = []
        input_samples = []
        logger.info(
            "WARNING: Skipping peak calling (no peak calling method selected)..."
        )

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
        "bai": "results/mapped/{wildcards.sample}.bl.bam.bai".format(
            wildcards=wildcards
        ),
        "multiqc": "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",  # "results/qc/pre_trim/multiqc_data/multiqc_general_stats.txt",
        "egs": "results/effective_genome_sizes/effective_genome_sizes.csv",
    }
    # Add additional input files depending on config file
    if config["spike_in"]["apply_spike_in"]:
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
    r = config["deeptools"]["matrix"]["regionsFileName"]
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
        raise ValueError(
            f"ERROR: config.peak_calling.macs2.broad must be True or False"
        )

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
            return False
    else:
        return False


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
        "multiqc": "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",  # "results/qc/pre_trim/multiqc_data/multiqc_general_stats.txt",
    }

    if config["remove_MT_seqs"]:
        _dict["mgs"] = f"resources/{genome}_mt_genome_size.txt"

    return _dict


def peak_fdrs(peak_mode):
    """
    Returns the MACS2 FDR cutoffs (fdr wildcard values) as strings:
    qvalue for narrow peaks and broad_cutoff for broad peaks.
    Each can be a single number or a list of numbers in config.yaml.
    """
    key = "qvalue" if peak_mode == "macs2_narrow" else "broad_cutoff"
    values = config["peak_calling"]["macs2"][key]
    if not isinstance(values, list):
        values = [values]

    fdrs = [str(v) for v in values]
    duplicates = sorted(set(f for f in fdrs if fdrs.count(f) > 1))
    if duplicates:
        raise ValueError(f"Duplicate MACS2 {key} cutoffs in config: {duplicates}")
    return fdrs


def broad_qvalue():
    """
    Returns the MACS2 qvalue used together with the broad cutoff for broad
    peaks: the strictest (smallest) value if qvalue is a list
    """
    qvalue = config["peak_calling"]["macs2"]["qvalue"]
    return min(qvalue) if isinstance(qvalue, list) else qvalue


def paired_end():
    """
    Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*fastq.gz")
    if len(reads) == 0:
        reads = glob.glob("reads/*fastq.gz")
    assert len(reads) != 0, "No fastq files found..."

    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith("_R1_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    elif fastq.endswith("_R2_001.fastq.gz"):
        logger.info("Paired-end reads detected...")
        return True
    else:
        logger.info("Single-end reads detected...")
        return False
