if config["peak_calling"]["macs2"]["use_macs2"]:
    if not config["peak_calling"]["macs2"]["broad"]:
        logger.info("MAC2S narrow peak calling selected...")
        if control_available():
            rule macs2_narrow:
                input: 
                    bam="results/mapped/{ip_sample}.bl.bam",
                    bai="results/mapped/{ip_sample}.bl.bam.bai",
                    cbam="results/mapped/{control_sample}.bl.bam",
                    cbais="results/mapped/{control_sample}.bl.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls",
                    peak=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.narrowPeak",
                    bed=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_summits.bed",
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
                    f"logs/macs2/narrow/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
                script:
                    "../scripts/macs2.py"
            

            rule macs2_narrow_replicates:
                input:
                    bams=expand("results/mapped/{ip_sample}.bl.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/{ip_sample}.bl.bam.bai", ip_sample=IP_SAMPLES),
                    cbams=expand("results/mapped/{control_sample}.bl.bam", control_sample=CONTROL_SAMPLES), 
                    cbais=expand("results/mapped/{control_sample}.bl.bam.bai", control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS_NO_CONTROL),
                    peak=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS_NO_CONTROL),
                    bed=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_summits.bed", conditions=CONDITIONS_NO_CONTROL),
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
                    expand(f"logs/macs2/narrow/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS_NO_CONTROL)
                script:
                    "../scripts/macs2_replicates.py"
        else:
            rule macs2_narrow:
                input: 
                    bam="results/mapped/{ip_sample}.bl.bam",
                    bai="results/mapped/{ip_sample}.bl.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.xls",
                    peak=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.narrowPeak",
                    bed=f"results/macs2_narrow/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_summits.bed",
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
                    f"logs/macs2/narrow/fdr{fdr}/{{ip_sample}}.log"
                script:
                    "../scripts/macs2.py"

            rule macs2_narrow_replicates:
                input:
                    bams=expand("results/mapped/{ip_sample}.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/{ip_sample}.bam.bai", ip_sample=IP_SAMPLES),
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS_NO_CONTROL),
                    peak=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS_NO_CONTROL),
                    bed=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_summits.bed", conditions=CONDITIONS_NO_CONTROL),
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
                    expand(f"logs/macs2/narrow/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS_NO_CONTROL)
                script:
                    "../scripts/macs2_replicates.py"
                

        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_narrow/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.narrowPeak", conditions=CONDITIONS_NO_CONTROL),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_narrow/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_narrow/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

    else:
        logger.info("MAC2S broad peak calling selected...")
        if control_available():
            logger.info("Control samples found...")
            rule macs2_broad:
                input: 
                    bam="results/mapped/{ip_sample}.bl.bam",
                    bai="results/mapped/{ip_sample}.bl.bam.bai",
                    cbam="results/mapped/{control_sample}.bl.bam",
                    cbais="results/mapped/{control_sample}.bl.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.xls",
                    peak=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.broadPeak",
                    gapped=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.gappedPeak",
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
                    f"logs/macs2/broad/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
                script:
                    "../scripts/macs2.py"

            
            rule macs2_broad_replicates:
                input:
                    bams=expand("results/mapped/{ip_sample}.bl.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/{ip_sample}.bl.bam.bai", ip_sample=IP_SAMPLES),
                    cbams=expand("results/mapped/{control_sample}.bl.bam", control_sample=CONTROL_SAMPLES), 
                    cbais=expand("results/mapped/{control_sample}.bl.bam.bai", control_sample=CONTROL_SAMPLES), 
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS_NO_CONTROL),
                    peak=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.broadPeak", conditions=CONDITIONS_NO_CONTROL),
                    gapped=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.gappedPeak", conditions=CONDITIONS_NO_CONTROL),
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
                    expand(f"logs/macs2/narrow/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS_NO_CONTROL)
                script:
                    "../scripts/macs2_replicates.py"

        else:
            rule macs2_broad:
                input: 
                    bam="results/mapped/{ip_sample}.bl.bam",
                    bai="results/mapped/{ip_sample}.bl.bam.bai",
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv",
                output:
                    xls=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.xls",
                    peak=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.broadPeak",
                    gapped=f"results/macs2_broad/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_peaks.gappedPeak",
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
                    f"logs/macs2/broad/fdr{fdr}/{{ip_sample}}.log"
                script:
                    "../scripts/macs2.py"

            rule macs2_broad_replicates:
                input:
                    bams=expand("results/mapped/{ip_sample}.bl.bam", ip_sample=IP_SAMPLES),
                    bais=expand("results/mapped/{ip_sample}.bl.bam.bai", ip_sample=IP_SAMPLES),
                    egs="results/effective_genome_sizes/effective_genome_sizes.csv"
                output:
                    xls=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls", conditions=CONDITIONS_NO_CONTROL),
                    peak=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.broadPeak", conditions=CONDITIONS_NO_CONTROL),
                    gapped=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.gappedPeak", conditions=CONDITIONS_NO_CONTROL),
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
                    expand(f"logs/macs2/narrow/fdr{fdr}/{{conditions}}.log", conditions=CONDITIONS_NO_CONTROL)
                script:
                    "../scripts/macs2_replicates.py"
            

        rule peak_annotation_plots:
            input:
                bed=expand(f"results/macs2_broad/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.broadPeak", conditions=CONDITIONS_NO_CONTROL),
                gtf=resources.gtf,
            output:
                dt=f"results/plots/macs2_broad/fdr{fdr}/peaks_distance_to_TSS.pdf",
                fd=f"results/plots/macs2_broad/fdr{fdr}/peak_distributions.pdf",
            log: f"logs/plots/fdr{fdr}/peak_annotation_plots.log"
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"
    
    rule annotate_peaks:
        input:
            xls=f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.xls",
            adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        output:
            bed=f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_peaks.bed",
            txt=f"results/{PEAK_MODE}/fdr{fdr}/{{conditions}}/{{conditions}}_annotated.peaks.txt",
        params:
            pm=PEAK_MODE,
            gtf=resources.gtf,
            extra=""
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            f"logs/annotate_peaks/{PEAK_MODE}/fdr{fdr}/{{conditions}}.log"
        conda:
            "../envs/diffbind.yaml"
        script:
            "../scripts/annotate_peaks.R"

    '''
    if run_diffbind():
        rule diffbind:
            input:
                **diffbind_input(wildcards)
            output:
                dba=f"results/diffbind/{PEAK_MODE}/fdr{fdr}/dba.RData",
                pca=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/PCA.pdf",
                sc=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/SampleCorrelation.pdf",
                pp=f"results/plots/diffbind/{PEAK_MODE}/fdr{fdr}/ProfilePLot.pdf",
                _dir=directory(f"results/diffbind/{PEAK_MODE}/fdr{fdr}/"),
                sh=f"results/diffbind/{PEAK_MODE}/fdr{fdr}/samplesheet.csv",
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
                "logs/diffbind/{peak_mode}/diffbind.log"
            script:
                "../scripts/diffbind.R"
    '''
    
  

elif config["peak_calling"]["htseq_count"]["use_htseq_count"]:
    rule call_peaks_htseq_count:
        input:
            bam="results/mapped/{sample}.bl.bam",
            bai="results/mapped/{sample}.bl.bam.bai",
            gtf=resources.gtf,
        output:
            counts="results/htseq_count/{sample}.tsv",
        params:
            mode=config["peak_calling"]["htseq_count"]["mode"],
            f=config["peak_calling"]["htseq_count"]["feature"],
            mapq=config["bowtie2"]["MAPQ_cutoff"],
            extra=config["peak_calling"]["htseq_count"]["extra"],
        threads: config["resources"]["deeptools"]["cpu"] * 2
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/htseq_count/{sample}.log"
        conda:
            "../envs/peaks.yaml"
        shell:
            "htseq-count "
            "-m {params.mode} "
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
            counts=expand("results/htseq_count/{sample}.tsv", sample=SAMPLES),
        output:
            xlsx="results/htseq_count/DESeq2/differential_peaks.xlsx",
            rdata="results/htseq_count/DESeq2/dds.RData",
        params:
            alpha=config["peak_calling"]["htseq_count"]["DESeq2"]["alpha"],
            fc=config["peak_calling"]["htseq_count"]["DESeq2"]["fc"],
            control=config["peak_calling"]["htseq_count"]["DESeq2"]["deseq2_apply_control"],
            cfo=config["peak_calling"]["htseq_count"]["DESeq2"]["cumulative_filter_out"],
            sg=config["peak_calling"]["htseq_count"]["DESeq2"]["smallest_group"],
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/peaks/DESeq2/differential_peaks.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/differential_peaks_DESeq2.R"
