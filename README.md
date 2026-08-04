# Snakemake workflow: `cut_and_run`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/cut_and_run/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/cut_and_run/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/739444701.svg)](https://zenodo.org/doi/10.5281/zenodo.10667876)

A Snakemake workflow for CUT&RUN sequencing data analysis, supporting paired-end and single-end reads, spike-in normalisation, and multiple peak calling strategies.

---

## Workflow overview

```
FASTQ reads
    │
    ▼
FastQC + MultiQC (QC)
    │
    ▼
Adapter trimming (Cutadapt or Trim Galore)
    │
    ▼
Reference genome download + Bowtie2 index
    │
    ▼
Bowtie2 alignment → SAM → BAM (sorted, filtered by MAPQ)
    │
    ▼
Blacklist filtering (samtools)
    │
    ├── [optional] Spike-in scale factor calculation (E. coli MG1655 or custom)
    │
    ▼
BigWig generation (deepTools bamCoverage; RPKM/CPM/BPM/RPGC or spike-in normalised)
    │
    ├── computeMatrix + plotHeatmap (deepTools)
    │
    ▼
Peak calling
    ├── MACS2 narrow or broad peaks (per-sample and replicate-merged)
    │       └── Consensus peaks → ChIPseeker annotation
    └── HTSeq-count + DESeq2 (gene-level differential analysis)
    │
    ▼
QC plots (PCA, scree plot, fragment length distribution)
```

---

## Deployment

### 1. Install Snakemake

Install Snakemake (≥ 8.25.0) using conda/mamba:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### 2. Clone the workflow

```bash
git clone https://github.com/niekwit/cut_and_run.git
cd cut_and_run
```

Alternatively, use [Snakedeploy](https://snakedeploy.readthedocs.io) to deploy the workflow into a new project directory:

```bash
pip install snakedeploy
snakedeploy deploy-workflow https://github.com/niekwit/cut_and_run . --tag main
```

### 3. Place FASTQ files

Copy or symlink your FASTQ files into the `reads/` directory.

**Paired-end** (auto-detected):

```
reads/{sample}_R1_001.fastq.gz
reads/{sample}_R2_001.fastq.gz
```

**Single-end**:

```
reads/{sample}.fastq.gz
```

> **Note:** Sample names must contain only alphanumeric characters and underscores, and must end with `_` followed by a single digit (e.g. `WT_H3K4me3_1`, `KO_IgG_2`).

---

## Configuration

### `config/samples.csv`

Describe your samples in `config/samples.csv`. Required columns:

| Column      | Description                                                                             |
| ----------- | --------------------------------------------------------------------------------------- |
| `sample`    | Sample name matching the FASTQ file name (without extension)                            |
| `genotype`  | Genotype of the sample                                                                  |
| `treatment` | Treatment condition (`None` if untreated)                                               |
| `factor`    | Antibody target or `IgG`/`input` for controls                                           |
| `control`   | _(optional)_ Name of the matched control sample (IgG/input) for peak calling            |
| `batch`     | _(optional)_ Numerical batch identifier for DiffBind                                    |
| `reference` | _(optional)_ Set to `yes` for the reference condition in DiffBind differential analysis |

Example:

| sample       | genotype | treatment | factor  | control  |
| ------------ | -------- | --------- | ------- | -------- |
| WT_H3K4me3_1 | WT       | None      | H3K4me3 | WT_IgG_1 |
| WT_H3K4me3_2 | WT       | None      | H3K4me3 | WT_IgG_2 |
| WT_IgG_1     | WT       | None      | IgG     |          |
| WT_IgG_2     | WT       | None      | IgG     |          |

### `config/config.yaml`

Key settings to adjust before running:

| Setting                                    | Description                                                     | Example           |
| ------------------------------------------ | --------------------------------------------------------------- | ----------------- |
| `genome`                                   | Reference genome (`hg38`, `hg19`, `mm39`, `mm38`, `dm3`, `dm6`) | `hg38`            |
| `ensembl_genome_build`                     | Ensembl build number                                            | `110`             |
| `spike_in.apply_spike_in`                  | Enable spike-in normalisation                                   | `False`           |
| `spike_in.genome`                          | Spike-in genome (e.g. `MG1655` for _E. coli_)                   | `MG1655`          |
| `spike_in.normalisation`                   | Spike-in normalisation method: `divide` or `multiply`           | `divide`          |
| `use_trim_galore`                          | Use Trim Galore instead of Cutadapt                             | `False`           |
| `remove_MT_seqs`                           | Exclude mitochondrial sequences from the genome index           | `True`            |
| `bowtie2.min_length`                       | Minimum fragment length for paired-end alignments               | `10`              |
| `bowtie2.max_length`                       | Maximum fragment length for paired-end alignments               | `700`             |
| `bowtie2.MAPQ_cutoff`                      | Minimum MAPQ score to retain an alignment                       | `5`               |
| `deeptools.bigwig.normalisation`           | BigWig normalisation (`RPKM`, `CPM`, `BPM`, `RPGC`, `None`)     | `RPKM`            |
| `deeptools.matrix.mode`                    | computeMatrix mode: `scale-regions` or `reference-point`        | `scale-regions`   |
| `peak_calling.macs2.use_macs2`             | Enable MACS2 peak calling                                       | `True`            |
| `peak_calling.macs2.broad`                 | Call broad peaks (e.g. for H3K27me3)                            | `False`           |
| `peak_calling.macs2.qvalue`                | MACS2 q-value cutoff for narrow peaks                           | `0.05`            |
| `peak_calling.htseq_count.use_htseq_count` | Enable HTSeq-count + DESeq2 instead of MACS2                    | `False`           |
| `resources.*`                              | CPUs and walltime (minutes) per rule                            | see `config.yaml` |

---

## Usage

Activate the Snakemake environment, then run from the workflow root directory.

**Dry run** (check what will be executed):

```bash
snakemake -n
```

**Local execution** (e.g. 8 cores):

```bash
snakemake --use-conda --cores 8
```

**HPC execution with SLURM** (using a Snakemake profile):

```bash
snakemake --use-conda --profile slurm
```

> All software dependencies are managed automatically via conda environments (`workflow/envs/`). The `--use-conda` flag is required.

### Spike-in normalisation

Set `spike_in.apply_spike_in: True` and provide the spike-in genome identifier and build number in `config.yaml`. Reads will be aligned to both the experimental and spike-in genomes, and scale factors calculated from spike-in alignment rates will be applied during BigWig generation.

### Peak calling modes

- **MACS2 narrow** (`broad: False`): for transcription factors.
- **MACS2 broad** (`broad: True`): for broad histone marks (e.g. H3K4me3).
- **HTSeq-count + DESeq2** (`use_htseq_count: True`): gene-level read counting followed by differential analysis with DESeq2. Use this when treating CUT&RUN data as a counting experiment rather than peak-based analysis.

---

## Outputs

| Path                                              | Description                                                       |
| ------------------------------------------------- | ----------------------------------------------------------------- |
| `results/qc/`                                     | FastQC and MultiQC reports                                        |
| `results/trimmed/`                                | Adapter-trimmed FASTQ files                                       |
| `results/mapped/`                                 | Blacklist-filtered, sorted BAM files                              |
| `results/bigwig/`                                 | Normalised BigWig files for genome browser                        |
| `results/deeptools/`                              | computeMatrix output                                              |
| `results/macs2_narrow/` or `results/macs2_broad/` | MACS2 peak files per sample and consensus peaks                   |
| `results/plots/`                                  | QC plots (PCA, scree, fragment lengths, heatmap, peak annotation) |
| `results/htseq_count/`                            | HTSeq-count tables and DESeq2 results (if enabled)                |
| `logs/snakemake/`                                 | Timestamped Snakemake log files                                   |

---

## Citation

If you use this workflow in your research, please cite the original repository and its Zenodo DOI:

> Wit, N. (2024). niekwit/cut_and_run. Zenodo. https://doi.org/10.5281/zenodo.10667876

Also cite the key tools used by the workflow:

- **Snakemake**: Mölder et al. (2021). _Sustainable data analysis with Snakemake_. F1000Research. https://doi.org/10.12688/f1000research.29032.2
- **Bowtie2**: Langmead & Salzberg (2012). _Fast gapped-read alignment with Bowtie 2_. Nature Methods. https://doi.org/10.1038/nmeth.1923
- **MACS2**: Zhang et al. (2008). _Model-based Analysis of ChIP-Seq (MACS)_. Genome Biology. https://doi.org/10.1186/gb-2008-9-9-r137
- **deepTools**: Ramírez et al. (2016). _deepTools2: a next generation web server for deep-sequencing data analysis_. Nucleic Acids Research. https://doi.org/10.1093/nar/gkw257
- **Cutadapt**: Martin (2011). _Cutadapt removes adapter sequences from high-throughput sequencing reads_. EMBnet.journal. https://doi.org/10.14806/ej.17.1.200
- **HTSeq**: Putri et al. (2022). _Analysing high-throughput sequencing data in Python with HTSeq 2.0_. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac166
- **DESeq2**: Love et al. (2014). _Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2_. Genome Biology. https://doi.org/10.1186/s13059-014-0550-8
- **ChIPseeker**: Yu et al. (2015). _ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization_. Bioinformatics. https://doi.org/10.1093/bioinformatics/btv145
