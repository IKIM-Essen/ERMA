# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule diamond_card:
    input:
        fasta=local("results/fastq/split/{sample}.part_{part}.fasta"),
        card=local("data/card_db/card_db.dmnd"),
    output:
        card_results=local("results/{sample}/{part}/card_results.txt"),
        overview_table=local("results/{sample}/{part}/overview_table.txt"),
    params:
        internal_threads=config["max_threads"],
    log:
        local("logs/diamond_card/{sample}_{part}.log"),
    conda:
        "../envs/diamond.yaml"
    threads: config["max_threads"]
    shell:
        """
        diamond blastx -d {input.card} -q {input.fasta} -o {output.card_results} --outfmt 6 --evalue 1e-5 --quiet --threads {params.internal_threads} 2> {log}
        echo -ne "Number of FastQ input reads,{wildcards.sample},{wildcards.part},$(cat {input.fasta}|grep -c '^>')\n" >> {output.overview_table}
        echo -ne "Diamond output hits,{wildcards.sample},{wildcards.part},$(cat {output.card_results}|wc -l)\n" >> {output.overview_table}
        """


if config["similarity_search_mode"] == "test":

    rule usearch_silva:
        input:
            overview_table=local("results/{sample}/{part}/overview_table.txt"),
            fasta=local("results/fastq/split/{sample}.part_{part}.fasta"),
            silva=local("data/silva_db/silva_seq.fasta"),
        output:
            silva_results=local(temp("results/{sample}/{part}/SILVA_results.txt")),
        log:
            local("logs/blast_silva/{sample}_{part}.log"),
        params:
            internal_threads=config["max_threads"],
            min_id=config["min_similarity"],
        conda:
            "../envs/usearch.yaml"
        threads: config["max_threads"]
        shell:
            """
            usearch -usearch_local {input.fasta} -db {input.silva} -blast6out {output.silva_results} -evalue 1e-5 -threads {params.internal_threads} -strand plus -mincols 200 > {log} 2>&1
            echo -ne "Usearch output hits,{wildcards.sample},{wildcards.part},$(cat {output.silva_results}|wc -l)\n" >> {input.overview_table}
            """


if config["similarity_search_mode"] == "full":

    rule usearch_silva:
        input:
            overview_table=local("results/{sample}/{part}/overview_table.txt"),
            fasta=local("results/fastq/split/{sample}.part_{part}.fasta"),
            silva=local("data/silva_db/silva_seq.fasta"),
        output:
            silva_results=local(temp("results/{sample}/{part}/SILVA_results.txt")),
        log:
            local("logs/blast_silva/{sample}_{part}.log"),
        params:
            internal_threads=config["max_threads"],
            min_id=config["min_similarity"],
        conda:
            "../envs/usearch.yaml"
        threads: config["max_threads"]
        shell:
            """
            usearch -usearch_local {input.fasta} -db {input.silva} -blast6out {output.silva_results} -evalue 1e-5 -threads {params.internal_threads} -strand both -mincols 200 2> {log}
            echo -ne "Usearch output hits,{wildcards.sample},{wildcards.part},$(cat {output.silva_results}|wc -l)\n" >> {input.overview_table}
            """


rule integrate_blast_data:
    input:
        overview_table=local("results/{sample}/{part}/overview_table.txt"),
        card_results=local("results/{sample}/{part}/card_results.txt"),
        silva_results=local("results/{sample}/{part}/SILVA_results.txt"),
        aro_mapping=local("data/card_db/aro_index.tsv"),
    output:
        intermed_card_results=local(
            temp("results/{sample}/{part}/intermed_card_results.csv")
        ),
        intermed_silva_results=local(
            temp("results/{sample}/{part}/intermed_silva_results.csv")
        ),
        integrated_data=local(
            temp("results/{sample}/{part}/integrated_filtered_results.csv")
        ),
    log:
        local("logs/integrate_blast_data/{sample}_{part}.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/integrate_blast_data.py"


rule filter_blast_results:
    input:
        overview_table=local("results/{sample}/{part}/overview_table.txt"),
        integrated_data=local("results/{sample}/{part}/integrated_filtered_results.csv"),
    output:
        filtered_data=local("results/{sample}/{part}/filtered_results.csv"),
    params:
        min_similarity=config["min_similarity"],
    log:
        local("logs/filter_blast_results/{sample}_{part}.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/filter_blast_results.py"


rule gzip_intermediates:
    input:
        silva_results=local("results/{sample}/{part}/SILVA_results.txt"),
        card_results=local("results/{sample}/{part}/card_results.txt"),
        filt_data=local("results/{sample}/{part}/filtered_results.csv"),
    output:
        silva_zip=local("results/{sample}/{part}/SILVA_results.txt.gz"),
        card_zip=local("results/{sample}/{part}/card_results.txt.gz"),
        filt_data_zip=local("results/{sample}/{part}/filtered_results.csv.gz"),
        checkpoint=local(temp("results/{sample}/{part}/checkpoint.txt")),
    log:
        local("logs/gzip_blast_results/{sample}_{part}.log"),
    conda:
        "../envs/python.yaml"
    shell:
        """
        gzip {input.silva_results} 2> {log}
        gzip {input.card_results} 2>> {log}
        gzip {input.filt_data} 2>> {log}
        touch {output.checkpoint}
        """


rule table_combined_genera_abundance:
    input:
        filtered_data=local(
            expand(
                "results/{sample}/{part}/filtered_results.csv.gz",
                sample=samples,
                part=get_numpart_list(),
            )
        ),
    output:
        csv=local("results/abundance/combined_genus_abundance.csv"),
    params:
        sample_name=samples,
    log:
        local("logs/genera_abundance_table.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/table_combined_genera_abundance.py"
