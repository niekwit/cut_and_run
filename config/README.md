# Samples.csv

Use this file to describe sample meta data:

`sample`: sample name that matches read file name without extension

`genotype`: genotype of the sample

`treatment`: treatment condition

`factor`: antibody used or input/IgG

`batch`: numerical value for batch 

`control`: (optional) name of the matched IgG/input sample

`reference`: `yes` for the reference condition (reserved for DiffBind, not currently run)

# downsample.csv (optional)

Use this file to downsample reads before trimming (with seqtk, seed set by `downsample_seed` in config.yaml):

`sample`: IP or control sample name

`percentage`: percentage of reads to keep (0-100; samples not listed or set to 100 keep all reads)
