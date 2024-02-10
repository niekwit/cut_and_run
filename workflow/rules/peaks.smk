if config["peak_calling"]["macs2"]["use_macs2"]:
    if not config["peak_calling"]["macs2"]["broad"]:
        sys.stderr.write("MAC2S narrow peak calling selected...\n")
        rule call_peaks_macs2_narrow:
            input: 
                unpack(macs2_input)
            output:
                xls="results/peaks/narrow/{bw_input_dir}/{ip_sample}/{ip_sample}_peaks.xls",
                peak="results/peaks/narrow/{bw_input_dir}/{ip_sample}/{ip_sample}_peaks.narrowPeak",
                bed="results/peaks/narrow/{bw_input_dir}/{ip_sample}/{ip_sample}_summits.bed",
                #unpack(macs2_output) don't think you can use function in output
            params:
                genome=resources.genome,
                outdir= lambda wc, output: os.path.dirname(output[0]),
                mode="narrow",
                qvalue=config["peak_calling"]["macs2"]["qvalue"],
                control=config["peak_calling"]["control_available"],
                extra="",
            threads: config["resources"]["macs2"]["cpu"]
            resources:
                runtime=config["resources"]["macs2"]["time"]
            log:
                "logs/macs2/narrow/{bw_input_dir}/{ip_sample}.log"
            script:
                "../scripts/macs2.py"
    else:
        sys.stderr.write("MAC2S broad peak calling selected...\n")
        rule call_peaks_macs2_broad:
            input: 
                unpack(macs2_input)
            output:
                xls="results/peaks/broad/{bw_input_dir}/{ip_sample}/{ip_sample}_peaks.xls",
                peak="results/peaks/broad/{bw_input_dir}/{ip_sample}/{ip_sample}_peaks.broadPeak",
                bed="results/peaks/broad/{bw_input_dir}/{ip_sample}/{ip_sample}_summits.bed",
                gapped="results/peaks/broad/{bw_input_dir}/{ip_sample}/{ip_sample}_peaks.gappedPeak",
                #unpack(macs2_output) don't think you can use a function in output
            params:
                genome=resources.genome,
                outdir= lambda wc, output: os.dirname(output[0]),
                qvalue=config["peak_calling"]["macs2"]["qvalue"],
                bc=config["peak_calling"]["macs2"]["broad_cutoff"],
                control=config["peak_calling"]["control_available"],
                mode="broad",
                extra="",
            threads: config["resources"]["macs2"]["cpu"]
            resources:
                runtime=config["resources"]["macs2"]["time"]
            log:
                "logs/macs2/broad/{bw_input_dir}/{ip_sample}.log"
            script:
                "../scripts/macs2.py"

    rule diffbind: # Compute differential peaks using DiffBind
        input:
            unpack(diffbind_input)
        output:
            dba="dba.RData",
        params:
            control=config["peak_calling"]["control_available"],
        script:
            "scripts/diffbind.R"

elif config["peak_calling"]["htseq_count"]["use_htseq_count"]:
    rule call_peaks_htseq_count:
        input:
            bam="results/mapped/{bw_input_dir}/{sample}.bam",
            gtf=resources.gtf,
        output:
            counts="results/peaks/htseq_count/{bw_input_dir}/{sample}.tsv",
            #anno_counts="results/peaks/htseq_count/{bw_input_dir}/{sample}_annotated.tsv"
        params:
            mode=config["peak_calling"]["htseq_count"]["mode"],
            f=config["peak_calling"]["htseq_count"]["feature"],
            mapq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["peak_calling"]["htseq_count"]["extra"],
        threads: config["resources"]["deeptools"]["cpu"] * 2
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/htseq_count/{bw_input_dir}/{sample}.log"
        conda:
            "../envs/peaks.yaml"
        shell:
            "htseq-count "
            "-m union "
            "-f bam " # data is bam format
            "-r pos " # bam is sorted on position
            "-s no " # Cut & Run data is not stranded
            "-t {params.f} " # get signal over whole gene, instead of just exons
            "-i gene_id " # use gene_id as identifier
            "-a {params.mapq} " # MAPQ cutoff
            "--additional-attr=gene_name " # use for annotation
            "--additional-attr=gene_biotype " # use for annotation
            "-n {threads} "
            "{params.extra} "
            "{input.bam} "
            "{input.gtf} "
            "2> {log} | "
            r"sed 's/\t\t/\tNA\t/g' > {output.counts}" # replace empty fields with NA (genes with no gene_name attributes)


    rule differential_peaks_DESeq2:
        input:
            counts=expand("results/peaks/htseq_count/{bw_input_dir}/{sample}.tsv", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
        output:
            xlsx="results/peaks/DESeq2/{bw_input_dir}/differential_peaks.xlsx",
            rdata="results/peaks/DESeq2/{bw_input_dir}/dds.RData",
        params:
            alpha=config["peak_calling"]["htseq_count"]["DESeq2"]["alpha"],
            control=config["peak_calling"]["htseq_count"]["DESeq2"]["deseq2_apply_control"],
            cfo=config["peak_calling"]["htseq_count"]["DESeq2"]["cumulative_filter_out"],
            sg=config["peak_calling"]["htseq_count"]["DESeq2"]["smallest_group"],
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/DESeq2/differential_peaks_{bw_input_dir}.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/differential_peaks_DESeq2.R"

else:
    sys.stderr.write("WARNING: Skipping peak calling (no peak calling method selected)...\n")
    