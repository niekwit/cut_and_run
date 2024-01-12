import pysam
import os
import pandas as pd

# Get number of threads for samtools
threads = str(snakemake.threads)

# Get all non-spike-in BAM files
bams = snakemake.input["bam"]

# Get all spike-in BAM files
spike_in_bams = snakemake.input["si_bam"]

# Empty data frame to store spike-in read count and scale factors
df = pd.DataFrame(columns=["sample", "sample_reads", "spike_in_reads", "spike_in_reads_per_sample_reads", "scale_factor"])

# Calculate ratio of non-spike-in to spike-in reads
for bam, spike_in_bam in zip(bams, spike_in_bams):
    # Get sample name and add to data frame
    sample = os.path.basename(bam).split(".")[0]
    df.loc[sample, "sample"] = sample
        
    # Get number of spike-in reads and add to data frame
    spike_in_reads = int(pysam.view("-c", "-@", threads, spike_in_bam))
    df.loc[sample, "spike_in_reads"] = spike_in_reads
    
    # Get number of non-spike-in reads and add to data frame
    non_spike_in_reads = int(pysam.view("-c","-@", threads, bam))
    df.loc[sample, "sample_reads"] = non_spike_in_reads
    
    # Calculate ratio and add to data frame
    ratio =  spike_in_reads / non_spike_in_reads
    df.loc[sample, "spike_in_reads_per_sample_reads"] = ratio
    
# set scale factor with lowest spike-in read count to 1 and correct other accordingly
df["scale_factor"] = df["spike_in_reads_per_sample_reads"].min() / df["spike_in_reads_per_sample_reads"] 

# Write scale factors to file
df.to_csv(snakemake.output[0], index=False)

