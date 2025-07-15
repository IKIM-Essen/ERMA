# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule get_16S_db:
    output:
        seq="data/silva_db/silva_seq_RNA.fasta.gz",
    params:
        path=get_silva_db_dir(),
        seq=config["silva"]["download_path_seq"],
    conda:
        "../envs/python.yaml"
    log:
        "logs/get_silva_db/get_silva_db.log",
    shell:
        r"""
        mkdir -p {params.path}; 
        if [[ "{params.seq}" == http* ]]; then 
            echo "Downloading from {params.seq}" >> {log}; 
            wget -O {output.seq} {params.seq} 2>> {log}
        else
            echo "Copying from local path {params.seq}" >> {log}; 
            cp {params.seq} {output.seq} 
        fi
        """


rule unzip_silva_db:
    input:
        seq="data/silva_db/silva_seq_RNA.fasta.gz",
    output:
        seq=temp("data/silva_db/silva_seq_RNA.fasta"),
    log:
        "logs/unzip_silva_db/get_silva_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        gzip -dk {input.seq} 2> {log};
        """


rule translate_silva_db:
    input:
        seq="data/silva_db/silva_seq_RNA.fasta",
    output:
        seq=temp("data/silva_db/silva_seq.fasta"),
    log:
        "logs/translate_silva_db/translate_silva_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        "seqtk seq -r {input.seq} > {output.seq}"


rule get_card_db:
    output:
        seq="data/card_db/card_seq.tar.bz2",
    params:
        path=get_card_db_dir(),
        seq=config["card"]["download_path"],
    log:
        "logs/get_card_db/get_silva_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        mkdir -p {params.path}; 
        if [[ "{params.seq}" == http* ]]; then 
            echo "Downloading from {params.seq}" >> {log}; 
            wget -O {output.seq} {params.seq} 2>> {log}
        else
            echo "Copying from local path {params.seq}" >> {log}; 
            cp {params.seq} {output.seq} 
        fi
        """


rule unzip_card_db:
    input:
        seq="data/card_db/card_seq.tar.bz2",
    output:
        seq="data/card_db/protein_fasta_protein_homolog_model.fasta",
        aro_mapping="data/card_db/aro_index.tsv",
    params:
        path=get_card_db_dir(),
    log:
        "logs/unzip_card_db/unzip_card_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo {params.path}
        tar -xvjf {input.seq} -C {params.path} 2> {log};
        """


rule makeblastdb_card:
    input:
        seq="data/card_db/protein_fasta_protein_homolog_model.fasta",
    output:
        db="data/card_db/card_db.dmnd",
    params:
        path=get_card_db_dir(),
    log:
        "logs/diamond_makedb_card/makedb_card.log",
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        diamond makedb --in {input.seq} -d {params.path}/card_db 2> {log};
        """
