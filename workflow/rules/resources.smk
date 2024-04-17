rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    cache: True
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    threads: config["resources"]["samtools"]["cpu"]
    resources: 
        runtime=config["resources"]["samtools"]["time"]
    conda:
        "../envs/mapping.yaml"
    script:
        "../scripts/get_resource.sh"


rule index_fasta:
        input:
            fasta=resources.fasta,
        output:
            index=f"{resources.fasta}.fai",
        log:
            "logs/resources/index_fasta.log"
        threads: config["resources"]["samtools"]["cpu"]
        resources: 
            runtime=config["resources"]["samtools"]["time"]
        params:
            extra="",  # optional params string
        wrapper:
            f"{wrapper_version}/bio/samtools/faidx"


if config["remove_MT_seqs"]:
    # based on https://bioinformatics.stackexchange.com/questions/16378/remove-from-multi-fasta-by-sequence-id?noredirect=1&lq=1
    rule remove_MT_seq_from_fasta:
        input:
            fasta=resources.fasta,
            fai=f"{resources.fasta}.fai",
        output:
            fasta=resources.nomt_fasta,
        params:
            mt=mitochondrial_genome_name(),
        log:
            "logs/resources/remove_MT_seq_from_fasta.log"
        threads: 1
        resources: 
            runtime=15
        conda:
            "../envs/mapping.yaml"
        shell:
            "awk '{{print $1}}' {input.fai} | " # Print field with chromosome name
            "grep -v '{params.mt}' | " # remove MT genome
            "xargs samtools faidx {input.fasta} > {output}"  # create fasta without MT sequence


    rule MT_genome_size:
        input:
            fai=f"{resources.fasta}.fai",
        output:
            size=f"resources/{genome}_mt_genome_size.txt",
        params:
            mt=mitochondrial_genome_name(),
        threads: 1
        resources: 
            runtime=10
        log:
            "logs/resources/MT_genome_size.log"
        conda:
            "../envs/mapping.yaml"
        shell:
            "grep '^{params.mt}' {input} | awk '{{print $2}}' > {output} 2> {log} "

    '''
    rule convert_mt_fasta_to_bed:
        input:
            fasta=resources.mt_fasta,
        output:
            bed=resources.mt_bed,
        log:
            "logs/resources/convert_mt_fasta_to_bed.log"
        conda:
            "../envs/mapping.yaml"
        shell:
            "faidx --transform bed {input} > {output} 2> {log}"
    '''

if config["spike-in"]["apply_spike_in"]:
    use rule get_fasta as get_spike_in_fasta with:
        output:
            resources_spike_in.fasta,
        params:
            url=resources_spike_in.fasta_url,
        log:
            "logs/resources/get_spike_in_fasta.log"
        

use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


use rule get_fasta as get_black_list with:
        output:
            resources.blacklist,
        params:
            url=resources.blacklist_url,
        log:
            "logs/resources/get_black_list.log"


rule convert2ensembl:
    input:
        resources.blacklist,
    output:
        resources.ensembl_blacklist,
    threads: 1
    resources: 
        runtime=10
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/resources/convert_blacklist2ensembl.log"
    shell:
        "sed 's/^chr//' {input} > {output} 2> {log}"


rule bowtie2_build:
    input:
        ref=fasta(config, resources),
    output:
        multiext(
            f"resources/bowtie2_index/{genome}_{resources.build}/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    cache: True
    log:
        "logs/bowtie2_build_genome/build.log",
    params:
        extra="",  # optional parameters
    threads: config["resources"]["index"]["cpu"]
    resources:
        runtime=config["resources"]["index"]["time"],
    wrapper:
        f"{wrapper_version}/bio/bowtie2/build"


rule chrom_sizes:
    input:
        fa=resources.fasta,
        fai=f"{resources.fasta}.fai",
    output:
        f"resources/{resources.genome}_chrom.sizes",
    log:
        "logs/resources/chrom_sizes.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/mapping.yaml"
    shell:
        "awk '{{print $1,$2}}' {input.fai} | "
        r"sed 's/ /\t/' > {output}"


rule create_annotation_file:
    input:
        gtf=resources.gtf,
    output:
        rdata=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
    log:
        "logs/resources/create_annotation_file.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/diffbind.yaml"
    script:
        "../scripts/create_annotation_file.R"