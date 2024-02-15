if config["peak_calling"]["macs2"]["use_macs2"]:
    if not config["peak_calling"]["macs2"]["broad"]:
        sys.stderr.write("MAC2S narrow peak calling selected...\n")
        if control_available():
            rule call_peaks_macs2_narrow:
                input: 
                    bam="results/mapped/bl_removed/{ip_sample}.bam",
                    bai="results/mapped/bl_removed/{ip_sample}.bam.bai",
                    cbam="results/mapped/bl_removed/{control_sample}.bam",
                    cbais="results/mapped/bl_removed/{control_sample}.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_peaks.xls",
                    peak="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_peaks.narrowPeak",
                    bed="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_summits.bed",
                    #flag=touch("results/peaks/narrow/macs2.done"),
                params:
                    genome=resources.genome,
                    outdir= lambda wc, output: os.path.dirname(output[0]),
                    mode="narrow",
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"]
                conda:
                    "../envs/peaks.yaml",
                log:
                    "logs/macs2/narrow/bl_removed/{ip_sample}_vs_{control_sample}.log"
                script:
                    "../scripts/macs2.py"
            

            rule call_peaks_macs2_narrow_replicates:
                input:
                    bams=expand("results/mapped/bl_removed/{ip_sample}.bam", bw_input_dir=BW_INPUT_DIR, ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/bl_removed/{ip_sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, ip_sample=IP_SAMPLES),
                    cbams=expand("results/mapped/bl_removed/{control_sample}.bam", bw_input_dir=BW_INPUT_DIR, control_sample=CONTROL_SAMPLES), 
                    cbais=expand("results/mapped/bl_removed/{control_sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_peaks.xls", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                    peak=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_peaks.narrowPeak", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                    bed=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_summits.bed", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                params:
                    genome=resources.genome,
                    mode="narrow",
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"] * 4
                conda:
                    "../envs/peaks.yaml",
                log:
                    expand("logs/macs2/narrow/bl_removed/{condition}.log", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS)
                script:
                    "../scripts/macs2_replicates.py"
        else:
            rule call_peaks_macs2_narrow:
                input: 
                    bam="results/mapped/bl_removed/{ip_sample}.bam",
                    bai="results/mapped/bl_removed/{ip_sample}.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_peaks.xls",
                    peak="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_peaks.narrowPeak",
                    bed="results/peaks/narrow/bl_removed/{ip_sample}/{ip_sample}_summits.bed",
                    #flag=touch("results/peaks/narrow/macs2.done"),
                params:
                    genome=resources.genome,
                    outdir= lambda wc, output: os.path.dirname(output[0]),
                    mode="narrow",
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"]
                conda:
                    "../envs/peaks.yaml",
                log:
                    "logs/macs2/narrow/bl_removed/{ip_sample}.log"
                script:
                    "../scripts/macs2.py"

            rule call_peaks_macs2_narrow_replicates:
                input:
                    bams=expand("results/mapped/bl_removed/{ip_sample}.bam", bw_input_dir=BW_INPUT_DIR, ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/bl_removed/{ip_sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, ip_sample=IP_SAMPLES),
                    #cbams=expand("results/mapped/bl_removed/{control_sample}.bam", bw_input_dir=BW_INPUT_DIR, control_sample=CONTROL_SAMPLES), 
                    #cbais=expand("results/mapped/bl_removed/{control_sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_peaks.xls", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                    peak=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_peaks.narrowPeak", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                    bed=expand("results/peaks/narrow/bl_removed/{condition}/{condition}_summits.bed", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS),
                params:
                    genome=resources.genome,
                    mode="narrow",
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"] * 4
                conda:
                    "../envs/peaks.yaml",
                log:
                    expand("logs/macs2/narrow/bl_removed/{condition}.log", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS)
                script:
                    "../scripts/macs2_replicates.py"

    else:
        sys.stderr.write("MAC2S broad peak calling selected...\n")
        if control_available():
            rule call_peaks_macs2_broad:
                input: 
                    bam="results/mapped/bl_removed/{ip_sample}.bam",
                    bai="results/mapped/bl_removed/{ip_sample}.bam.bai",
                    cbam="results/mapped/bl_removed/{control_sample}.bam",
                    cbais="results/mapped/bl_removed/{control_sample}.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_peaks.xls",
                    peak="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_peaks.broadPeak",
                    gapped="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_vs_{control_sample}_peaks.gappedPeak",
                    #flag=touch("results/peaks/broad/macs2.done"),
                params:
                    genome=resources.genome,
                    outdir= lambda wc, output: os.path.dirname(output[0]),
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    bc=config["peak_calling"]["macs2"]["broad_cutoff"],
                    control=control_available(),
                    mode="broad",
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"]
                conda:
                    "../envs/peaks.yaml",
                log:
                    "logs/macs2/broad/bl_removed/{ip_sample}_vs_{control_sample}.log"
                script:
                    "../scripts/macs2.py"

            rule call_peaks_macs2_broad_replicates:
                input:
                    bams=expand("results/mapped/bl_removed/{ip_sample}.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/bl_removed/{ip_sample}.bam.bai", ip_sample=IP_SAMPLES),
                    cbams=expand("results/mapped/bl_removed/{control_sample}.bam", control_sample=CONTROL_SAMPLES), 
                    cbais=expand("results/mapped/bl_removed/{control_sample}.bam.bai", control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.xls", condition=CONDITIONS),
                    peak=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.broadPeak", condition=CONDITIONS),
                    gapped=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.gappedPeak", condition=CONDITIONS),
                params:
                    genome=resources.genome,
                    mode="broad",
                    bc=config["peak_calling"]["macs2"]["broad_cutoff"],
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"] * 4
                conda:
                    "../envs/peaks.yaml",
                log:
                    expand("logs/macs2/narrow/bl_removed/{condition}.log", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS)
                script:
                    "../scripts/macs2_replicates.py"
        else:
            rule call_peaks_macs2_broad:
                input: 
                    bam="results/mapped/bl_removed/{ip_sample}.bam",
                    bai="results/mapped/bl_removed/{ip_sample}.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_peaks.xls",
                    peak="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_peaks.broadPeak",
                    gapped="results/peaks/broad/bl_removed/{ip_sample}/{ip_sample}_peaks.gappedPeak",
                    #flag=touch("results/peaks/broad/macs2.done"),
                params:
                    genome=resources.genome,
                    outdir= lambda wc, output: os.path.dirname(output[0]),
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    bc=config["peak_calling"]["macs2"]["broad_cutoff"],
                    control=control_available(),
                    mode="broad",
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"]
                conda:
                    "../envs/peaks.yaml",
                log:
                    "logs/macs2/broad/bl_removed/{ip_sample}.log"
                script:
                    "../scripts/macs2.py"

            rule call_peaks_macs2_broad_replicates:
                input:
                    bams=expand("results/mapped/bl_removed/{ip_sample}.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/bl_removed/{ip_sample}.bam.bai", ip_sample=IP_SAMPLES),
                    #cbams=expand("results/mapped/bl_removed/{control_sample}.bam", control_sample=CONTROL_SAMPLES), 
                    #cbais=expand("results/mapped/bl_removed/{control_sample}.bam.bai", control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.xls", condition=CONDITIONS),
                    peak=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.broadPeak", condition=CONDITIONS),
                    gapped=expand("results/peaks/broad/bl_removed/{condition}/{condition}_peaks.gappedPeak", condition=CONDITIONS),
                params:
                    genome=resources.genome,
                    mode="broad",
                    bc=config["peak_calling"]["macs2"]["broad_cutoff"],
                    qvalue=config["peak_calling"]["macs2"]["qvalue"],
                    control=control_available(),
                    csv=csv,
                    extra=config["peak_calling"]["macs2"]["extra"],
                threads: config["resources"]["macs2"]["cpu"]
                resources:
                    runtime=config["resources"]["macs2"]["time"] * 4
                conda:
                    "../envs/peaks.yaml",
                log:
                    expand("logs/macs2/narrow/bl_removed/{condition}.log", bw_input_dir=BW_INPUT_DIR, condition=CONDITIONS)
                script:
                    "../scripts/macs2_replicates.py"
    
    rule annotate_peaks:
        input:
            xls="results/peaks/{peak_mode}/bl_removed/{condition}/{condition}_peaks.xls",
            adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        output:
            bed="results/peaks/{peak_mode}/bl_removed/{condition}/{condition}_peaks.bed",
            txt="results/peaks/{peak_mode}/bl_removed/{condition}/{condition}_annotated.peaks.txt",
            bar="results/plots/peaks/{peak_mode}/bl_removed/{condition}/{condition}_annotated_peaks_bar.pdf",
            tss="results/plots/peaks/{peak_mode}/bl_removed/{condition}/{condition}_annotated_peaks_tss.pdf",
            #enrichment="results/plots/peaks/{peak_mode}/bl_removed/{condition}/{condition}_annotated_peaks_enrichment.pdf",
        params:
            genome=genome,
            pm=PEAK_MODE,
            extra=""
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/annotate_peaks/{peak_mode}/bl_removed/{condition}/{condition}.log"
        conda:
            "../envs/diffbind.yaml"
        script:
            "../scripts/annotate_peaks.R"

    '''
    if run_diffbind():
        rule diffbind:
            input:
                bam=expand("results/mapped/bl_removed/{sample}.bam", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
                bai=expand("results/mapped/bl_removed/{sample}.bam.bai", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
                flag=touch(f"results/peaks/{PEAK_MODE}/macs2.done"),
            output:
                dba=f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/dba.RData",
                pca=f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/PCA.pdf",
                sc=f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/SampleCorrelation.pdf",
                pp=f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/ProfilePLot.pdf",
                _dir=directory(f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/"),
                sh=f"results/peaks/diffbind/{PEAK_MODE}/bl_removed/samplesheet.csv",
            params:
                control=control_available(),
                genome=genome,
            threads: config["resources"]["diffbind"]["cpu"]
            resources:
                runtime=config["resources"]["diffbind"]["time"],
                mem_mb=120000,
            conda:
                "../envs/diffbind.yaml"
            log:
                "logs/diffbind/{peak_mode}/bl_removed/diffbind.log"
            script:
                "../scripts/diffbind.R"
    '''
    
  

elif config["peak_calling"]["htseq_count"]["use_htseq_count"]:
    rule call_peaks_htseq_count:
        input:
            bam="results/mapped/bl_removed/{sample}.bam",
            gtf=resources.gtf,
        output:
            counts="results/peaks/htseq_count/bl_removed/{sample}.tsv",
            #anno_counts="results/peaks/htseq_count/bl_removed/{sample}_annotated.tsv"
        params:
            mode=config["peak_calling"]["htseq_count"]["mode"],
            f=config["peak_calling"]["htseq_count"]["feature"],
            mapq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["peak_calling"]["htseq_count"]["extra"],
        threads: config["resources"]["deeptools"]["cpu"] * 2
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/htseq_count/bl_removed/{sample}.log"
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
            counts=expand("results/peaks/htseq_count/bl_removed/{sample}.tsv", bw_input_dir=BW_INPUT_DIR, sample=SAMPLES),
        output:
            xlsx="results/peaks/DESeq2/bl_removed/differential_peaks.xlsx",
            rdata="results/peaks/DESeq2/bl_removed/dds.RData",
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
            "logs/peaks/DESeq2/differential_peaks_bl_removed.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/differential_peaks_DESeq2.R"

#rule pathway_analysis:
#    input: