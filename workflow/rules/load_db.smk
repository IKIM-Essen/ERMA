# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule get_16S_db:
    output:
        seq=local("data/silva_db/silva_seq_RNA.fasta.gz"),
    params:
        path=get_silva_db_dir(),
        seq=config["download_silva"],
    conda:
        "../envs/python.yaml"
    log:
        local("logs/get_silva_db/get_silva_db.log"),
    shell:
        r"""
        mkdir -p {params.path}; 
        if [[ "{params.seq}" == http* ]]; then 
            echo "Downloading from {params.seq}"; 
            wget -O {output.seq} {params.seq} 2>> {log}
        else
            echo "Copying from local path {params.seq}"; 
            cp {params.seq} {output.seq} 2>> {log} 
        fi
        """


rule get_card_db:
    output:
        seq=local("data/card_db/card_seq.tar.bz2"),
    params:
        path=get_card_db_dir(),
        seq=config["download_card"],
    log:
        local("logs/get_card_db/get_card_db.log"),
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        mkdir -p {params.path}; 
        if [[ "{params.seq}" == http* ]]; then 
            echo "Downloading from {params.seq}"; 
            wget -O {output.seq} {params.seq} 2>> {log}
        else
            echo "Copying from local path {params.seq}"; 
            cp {params.seq} {output.seq} 2>> {log} 
        fi
        """


rule unzip_silva_db:
    input:
        seq=local("data/silva_db/silva_seq_RNA.fasta.gz"),
    output:
        seq=local(temp("data/silva_db/silva_seq_RNA.fasta")),
    log:
        local("logs/unzip_silva_db/get_silva_db.log"),
    conda:
        "../envs/python.yaml"
    shell:
        """
        gzip -dk {input.seq} 2> {log};
        """


rule translate_silva_db:
    input:
        seq=local("data/silva_db/silva_seq_RNA.fasta"),
    output:
        seq=local(temp("data/silva_db/silva_seq.fasta")),
    log:
        local("logs/translate_silva_db/translate_silva_db.log"),
    conda:
        "../envs/python.yaml"
    shell:
        "seqtk seq -r {input.seq} > {output.seq}"


rule unzip_card_db:
    input:
        seq=local("data/card_db/card_seq.tar.bz2"),
    output:
        aro=local(temp("data/card_db/unpacked/aro_index.tsv")),
        seq=local(temp("data/card_db/unpacked/protein_fasta_protein_homolog_model.fasta")),
        unpacked=local(temp(directory("data/card_db/unpacked"))),
    log:
        local("logs/unzip_card_db/unzip_card_db.log"),
    conda:
        "../envs/python.yaml"
    shell:
        """
        mkdir -p {output.unpacked}
        tar -xvjf {input.seq} -C {output.unpacked} 2> {log}
        """


if config["add_uniref_targets"]:
    rule merge_card_with_uniprot:
        input:
            card_fasta=local("data/card_db/unpacked/protein_fasta_protein_homolog_model.fasta")
        output:
            info_tsv=local(temp("data/card_db/protein_fasta_with_uniprot.tsv")),
            merged_fasta=local(temp("data/card_db/protein_fasta_with_uniprot.fasta")),
        params:
            cluster=config["uniprot_cluster"],
            targets=config["uniprot_targets"],
            size=config["max_entry_count"]
        conda:
            "../envs/python.yaml"
        log:
            "logs/uniprot_merge/uniprot_merge.log"
        script:
            "../scripts/fetch_and_merge_uniprot.py"


    rule normalize_merged_db:
        input:
            merged_fasta=local("data/card_db/protein_fasta_with_uniprot.fasta")
        output:
            norm_fasta=local(temp("data/card_db/protein_fasta_with_uniprot_norm.fasta"))
        conda:
            "../envs/python.yaml"
        log:
            "logs/uniprot_merge/uniprot_norm.log"
        shell:
            r"""
            seqkit rmdup -s {input.merged_fasta} 2> {log} \
            | seqtk seq -l 0 - > {output.norm_fasta} 2>> {log}
            """


    rule makeblastdb_card:
        input:
            unpacked=local("data/card_db/unpacked"),
            seq=local("data/card_db/protein_fasta_with_uniprot_norm.fasta"),
        output:
            db=local("data/card_db/card_db.dmnd"),
        params:
            path=get_card_db_dir(),
        log:
            local("logs/diamond_makedb_card/makedb_card.log"),
        conda:
            "../envs/diamond.yaml"
        shell:
            """
            diamond makedb --in {input.seq} -d {params.path}/card_db 2> {log};
            """


elif not config["add_uniref_targets"]:
    rule makeblastdb_card:
        input:
            unpacked=local("data/card_db/unpacked"),
            seq=local("data/card_db/unpacked/protein_fasta_protein_homolog_model.fasta"),
        output:
            db=local("data/card_db/card_db.dmnd"),
        params:
            path=get_card_db_dir(),
        log:
            local("logs/diamond_makedb_card/makedb_card.log"),
        conda:
            "../envs/diamond.yaml"
        shell:
            """
            diamond makedb --in {input.seq} -d {params.path}/card_db 2> {log};
            """
