FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="39a510ae3eae7658986d8f0a2f384fe10f50c50bb418a4043247905ab5859137"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/bedtools/intersect/environment.yaml
#   prefix: /conda-envs/f64f273a2db056abfb1bb96490c40bbe
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bedtools =2.31.1
#     - htslib =1.22
RUN mkdir -p /conda-envs/f64f273a2db056abfb1bb96490c40bbe
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/bedtools/intersect/environment.yaml /conda-envs/f64f273a2db056abfb1bb96490c40bbe/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/4f56e195d46da93fbd8fc814617891ce
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/4f56e195d46da93fbd8fc814617891ce
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/bowtie2/build/environment.yaml /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/cutadapt/pe/environment.yaml
#   prefix: /conda-envs/fe64d7cb6796b53876f4d2b0026849e0
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =5.1
RUN mkdir -p /conda-envs/fe64d7cb6796b53876f4d2b0026849e0
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/cutadapt/pe/environment.yaml /conda-envs/fe64d7cb6796b53876f4d2b0026849e0/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/7cdb9e3180904a671be6b170bf7651d2
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.7.2
RUN mkdir -p /conda-envs/7cdb9e3180904a671be6b170bf7651d2
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/fastqc/environment.yaml /conda-envs/7cdb9e3180904a671be6b170bf7651d2/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/84188cb48cfebe3cadc377f2f5001295
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.22
#     - snakemake-wrapper-utils =0.7.2
RUN mkdir -p /conda-envs/84188cb48cfebe3cadc377f2f5001295
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/faidx/environment.yaml /conda-envs/84188cb48cfebe3cadc377f2f5001295/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/93bab9ed6097f300f8667bc95cf9fd61
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.22
RUN mkdir -p /conda-envs/93bab9ed6097f300f8667bc95cf9fd61
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.0.0/bio/samtools/index/environment.yaml /conda-envs/93bab9ed6097f300f8667bc95cf9fd61/environment.yaml

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/c62fabd74158877ceb04451c2e72ec98
#   name: R
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.2
#     - r-reshape2=1.4.4
#     - r-rcolorbrewer=1.1_3
#     - r-ggrepel=0.9.4
#     - r-scales=1.3.0
#     - bioconductor-chipseeker=1.36.0
#     - bioconductor-genomicfeatures=1.52.1
#     - bioconductor-rtracklayer=1.60.0
#     - r-viridislite=0.4.2
#     - bioconductor-deseq2=1.40.2
#     - bioconductor-ihw=1.28.0
#     - r-openxlsx=4.2.5
#   #  - bioconductor-reactomepa=1.46
RUN mkdir -p /conda-envs/c62fabd74158877ceb04451c2e72ec98
COPY workflow/envs/R.yaml /conda-envs/c62fabd74158877ceb04451c2e72ec98/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/6bab4889727c9a54b84b4a193a411453
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - deeptools=3.5.4
#     - python=3.10
#     - pandas=2.1.4
#     - pysam=0.22.0
#     - wiggletools=1.2.11
#     - ucsc-wigtobigwig=447
#     - pybedtools=0.10.0
RUN mkdir -p /conda-envs/6bab4889727c9a54b84b4a193a411453
COPY workflow/envs/deeptools.yaml /conda-envs/6bab4889727c9a54b84b4a193a411453/environment.yaml

# Conda environment:
#   source: workflow/envs/diffbind.yaml
#   prefix: /conda-envs/fcc3fce9754fdad4cad3ccd18fc615fe
#   name: diffbind
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.2
#     - bioconductor-diffbind=3.12.0
#     - bioconductor-chipseeker=1.38.0
#     - bioconductor-genomicfeatures=1.54.1
#     - bioconductor-rtracklayer=1.62.0
#     - bioconductor-profileplyr=1.18.0
RUN mkdir -p /conda-envs/fcc3fce9754fdad4cad3ccd18fc615fe
COPY workflow/envs/diffbind.yaml /conda-envs/fcc3fce9754fdad4cad3ccd18fc615fe/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yaml
#   prefix: /conda-envs/2835662405b3a6fa87f9b0640552f9b7
#   name: mapping
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bowtie2=2.5.2
#     - samtools=1.19
#     - pigz=2.8
#     - multiqc=1.19
#     - python=3.10
#     - pyfaidx=0.8.0
RUN mkdir -p /conda-envs/2835662405b3a6fa87f9b0640552f9b7
COPY workflow/envs/mapping.yaml /conda-envs/2835662405b3a6fa87f9b0640552f9b7/environment.yaml

# Conda environment:
#   source: workflow/envs/peaks.yaml
#   prefix: /conda-envs/29a42af52cf5ad3813c83dad0ee6823e
#   name: peaks
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.10
#     - macs2=2.2.7.1
#     - htseq=2.0.4
#     - pandas=2.0.3
RUN mkdir -p /conda-envs/29a42af52cf5ad3813c83dad0ee6823e
COPY workflow/envs/peaks.yaml /conda-envs/29a42af52cf5ad3813c83dad0ee6823e/environment.yaml

# Conda environment:
#   source: workflow/envs/seqtk.yaml
#   prefix: /conda-envs/7e1fcde41edf6003683de16db2e301a2
#   name: seqtk
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - seqtk=1.4
#     - pigz=2.8
RUN mkdir -p /conda-envs/7e1fcde41edf6003683de16db2e301a2
COPY workflow/envs/seqtk.yaml /conda-envs/7e1fcde41edf6003683de16db2e301a2/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/f64f273a2db056abfb1bb96490c40bbe --file /conda-envs/f64f273a2db056abfb1bb96490c40bbe/environment.yaml && \
    conda env create --prefix /conda-envs/4f56e195d46da93fbd8fc814617891ce --file /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml && \
    conda env create --prefix /conda-envs/fe64d7cb6796b53876f4d2b0026849e0 --file /conda-envs/fe64d7cb6796b53876f4d2b0026849e0/environment.yaml && \
    conda env create --prefix /conda-envs/7cdb9e3180904a671be6b170bf7651d2 --file /conda-envs/7cdb9e3180904a671be6b170bf7651d2/environment.yaml && \
    conda env create --prefix /conda-envs/84188cb48cfebe3cadc377f2f5001295 --file /conda-envs/84188cb48cfebe3cadc377f2f5001295/environment.yaml && \
    conda env create --prefix /conda-envs/93bab9ed6097f300f8667bc95cf9fd61 --file /conda-envs/93bab9ed6097f300f8667bc95cf9fd61/environment.yaml && \
    conda env create --prefix /conda-envs/c62fabd74158877ceb04451c2e72ec98 --file /conda-envs/c62fabd74158877ceb04451c2e72ec98/environment.yaml && \
    conda env create --prefix /conda-envs/6bab4889727c9a54b84b4a193a411453 --file /conda-envs/6bab4889727c9a54b84b4a193a411453/environment.yaml && \
    conda env create --prefix /conda-envs/fcc3fce9754fdad4cad3ccd18fc615fe --file /conda-envs/fcc3fce9754fdad4cad3ccd18fc615fe/environment.yaml && \
    conda env create --prefix /conda-envs/2835662405b3a6fa87f9b0640552f9b7 --file /conda-envs/2835662405b3a6fa87f9b0640552f9b7/environment.yaml && \
    conda env create --prefix /conda-envs/29a42af52cf5ad3813c83dad0ee6823e --file /conda-envs/29a42af52cf5ad3813c83dad0ee6823e/environment.yaml && \
    conda env create --prefix /conda-envs/7e1fcde41edf6003683de16db2e301a2 --file /conda-envs/7e1fcde41edf6003683de16db2e301a2/environment.yaml && \
    conda clean --all -y