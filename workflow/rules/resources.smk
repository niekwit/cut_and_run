rule get_fasta:
    output:
        resources.fasta,
    retries: 3
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


if config["remove_MT_seqs"]:
    # based on https://bioinformatics.stackexchange.com/questions/16378/remove-from-multi-fasta-by-sequence-id?noredirect=1&lq=1
    rule index_fasta:
        input:
            fasta=resources.fasta,
        output:
            index=f"{resources.fasta}.fai",
        log:
            "logs/resources/index_fasta.log"
        params:
            extra="",  # optional params string
        wrapper:
            "v3.3.6/bio/samtools/faidx"
    
    
    rule remove_MT_seq_from_fasta:
        input:
            fasta=resources.fasta,
            fai=f"{resources.fasta}.fai",
        output:
            fasta=resources.nomt_fasta,
        log:
            "logs/resources/remove_MT_seq_from_fasta.log"
        conda:
            "../envs/mapping.yaml"
        shell:
            "awk '{{print $1}}' {input.fai} | " # Print field with chromosome name
            "grep -v 'MT' | " # remove MT genome
            "xargs samtools faidx {input.fasta} > {output}"  # create fasta without MT sequence


    rule MT_genome_size:
        input:
            fai=f"{resources.fasta}.fai",
        output:
            size=f"resources/{genome}_mt_genome_size.txt",
        log:
            "logs/resources/MT_genome_size.log"
        conda:
            "../envs/mapping.yaml"
        shell:
            "grep '^MT' {input} | awk '{{print $2}}' > {output} 2> {log} "

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


rule compress_resources:
    input:
        f=resources.fasta,
        g=resources.gtf,
        p="results/plots/pca.pdf",# dummy input to make sure this rule is executed at the end
    output:
        f"{resources.fasta}.gz",
        f"{resources.gtf}.gz",
    params:
        pigz_options="",
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    log:
        "logs/resources/compress_resources.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "pigz "
        "-p {threads} "
        "{params.pigz_options} "
        "{input.f} "
        "> {log} 2>&1"

