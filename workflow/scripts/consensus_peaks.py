import os
import pybedtools
import pandas as pd

bed_files = snakemake.input["beds"]
conditions = snakemake.params["conditions"]
out_bed_files = snakemake.output["bed_out"]
out_bed_files_intermediate = snakemake.output["bed_out_intermediate"]
k = snakemake.params["keep"]
extend_by = snakemake.params["extend_by"]
max_size = snakemake.params["max_size"]
chrom_sizes_file = snakemake.input["chrom_sizes"]

# Load chrom sizes into dictionary
chrom_sizes = {}
with open(chrom_sizes_file) as f:
    for line in f:
        (chrom, size) = line.split()
        chrom_sizes[chrom] = int(size)

for condition in conditions:   
    k = snakemake.params["keep"]
    
    log = [x for x in snakemake.log if condition in x][0]
    with open(log, "a") as log:
        # Empty data frame to store consensus peaks
        df = pd.DataFrame(columns=["chrom", "start", "end"])
        
        # Get all bed files that contain the condition
        beds = [bed for bed in bed_files if condition in bed]
        
        # Check if k is larger than the number of bed files
        # otherwise set it to the number of bed files
        if k > len(beds):
            k = len(beds)
            log.write(f"Warning: k is larger than the number of bed files for {condition}. Setting k to {k}.\n")
        
        # Load each individual bed file into a pandas data frame
        ind_peaks = []
        for bed in beds:
            tmp = pd.read_csv(bed, sep="\t", header=None, low_memory=False)
            ind_peaks.append(tmp)
        
        # Get overlapping regions between all bed files
        x = pybedtools.BedTool()
        input_list = [pybedtools.BedTool(bed).fn for bed in beds]
        consensus_peaks = x.multi_intersect(i=input_list)
            
        # Define column names
        # 5 base columns and one for each bed file
        names = ["chrom", "start", "end", "num", "list"]
        base_file_names = [os.path.basename(x) for x in input_list]
        [names.append(base_file_names[i]) for i in range(len(input_list))]
        
        # Convert to pandas data frame
        consensus_peaks = consensus_peaks.to_dataframe(names=names)

        # Save intermediate bed file
        bed_out_intermediate = [x for x in out_bed_files_intermediate if condition in x][0]
        consensus_peaks.to_csv(bed_out_intermediate, sep="\t", header=False, index=False)
        
        # Get total number of overlapping regions before filtering
        total = len(consensus_peaks)
        
        # Filter out regions that are not overlapping in at least k bed files
        # name column contains the number of overlapping regions
        consensus_peaks = consensus_peaks[consensus_peaks["num"] >= k]
        
        # Number of peaks not in k bed files
        skipped_peaks = total - len(consensus_peaks)
        
        # Instead of only the overlapping region, get the region that contains 
        # the whole area of the overlapping regions
        extended_peaks = 0
        for row in consensus_peaks.itertuples():
            chrom = row.chrom
            start = row.start
            end = row.end
        
            starts = []
            ends = []
            for x in ind_peaks:
                    x = x[(x[0] == chrom) & (x[1] <= int(end)) & (x[2] >= int(start))]
                    if not x.empty:
                        starts.append(x[1].min())
                        ends.append(x[2].max())
            start = min(starts)
            end = max(ends)
        
            # Check if regions are short enough to extend
            if end - start < max_size:
                start = start - extend_by
                if start < 1:
                    start = 1
                # Do not extend end beyond chromosome boundary
                end = end + extend_by
                if end > chrom_sizes[chrom]:
                    end = chrom_sizes[chrom]
                extended_peaks += 1
            
            # Add peak to df
            tmp = pd.DataFrame([[chrom, start, end]], columns=["chrom", "start", "end"])
            df = pd.concat([df,tmp])
        
        # Remove duplicate lines from df
        df = df.drop_duplicates()
        
        # Name peaks
        df["name"] = [f"peak_{i}" for i in range(1, len(df) + 1)]
        
        # Save consensus peaks to bed file
        bed_out = [x for x in out_bed_files if condition in x][0]
        df.to_csv(bed_out, sep="\t", header=False, index=False)
        
        log.write(f"Total number of peaks analysed: {total}\n")
        log.write(f"Skipped peaks: {skipped_peaks}\n")
        log.write(f"Total overlapping peaks: {len(consensus_peaks)}\n")
        log.write(f"Overlapping and extended peaks: {extended_peaks}\n")