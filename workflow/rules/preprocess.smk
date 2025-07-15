# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule convert_fastq_to_fasta:
    input:
        lambda wildcards: 
            os.path.join(config["fastq_dir"], f"{wildcards.sample}.fastq.gz")
    output:
        temp("results/fastq/{sample}.fasta"),
    log:
        "logs/convert_fastq_to_fasta/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"


rule split_fasta_file:
    input:
        "results/fastq/{sample}.fasta",
    output:
        temp("results/fastq/split/{sample}.part_{part}.fasta"),
    log:
        "logs/split_fasta/{sample}_{part}.log",
    conda:
        "../envs/python.yaml"
    params:
        outdir=config["base_dir"],
        num_parts=config["num_parts"],
    shell:
        """
        seqkit split2 --by-part {params.num_parts} -w 0 {input} --out-dir {params.outdir}/results/fastq/split/
        """


if config["seq_tech"] == "Illumina":

    rule prepare_fastqs:
        output:
            sentinel=temp("data/fastq/merged_done.txt"),
        log:
            "logs/prepare_fastqs_Illumina/prepare.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/prepare_fastq.sh"

elif config["seq_tech"] == "ONT":

    rule prepare_fastqs:
        output:
            sentinel=temp("data/fastq/merged_done.txt"),
        params:
            run_path=config["ONT"]["fastq_pass_path"],
            sample_name_path=config["ONT"]["sample_name_path"],
            target_fragment_length=config["ONT"]["target_fragment_length"],
            filter_intervall=config["ONT"]["filter_intervall"],
            output_dir="data/fastq/",
        log:
            "logs/prepare_fastqs_ONT/prepare.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/prepare_fastq.py"
