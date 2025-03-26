rule diamond_card:
    input:
        fasta = local("results/fastq/{sample}.part_{part}.fasta"),
        card = local("data/card_db/card_db.dmnd")
    output:
        card_results = local("results/{sample}/{part}/card_results.txt")
    params:
        internal_threads = config["max_threads"],        
    log:
        local("logs/diamond_card/{sample}_{part}.log")
    conda:
        "../envs/diamond.yaml" 
    threads: config["max_threads"]      
    shell:
        """
        diamond blastx -d {input.card} -q {input.fasta} -o {output.card_results} --outfmt 6 --evalue 1e-5 --threads {params.internal_threads} 2> {log}
        """

rule usearch_silva:
    input:
        fasta = local("results/fastq/{sample}.part_{part}.fasta"),
        silva = local("data/silva_db/silva_seq.fasta")
    output:
        silva_results = local("results/{sample}/{part}/SILVA_results.txt")
    log:
        local("logs/blast_silva/{sample}_{part}.log")
    params:
        internal_threads = config["max_threads"],
        min_id = config["min_similarity"]
    conda:
        "../envs/usearch.yaml"
    threads: config["max_threads"]                  
    shell:
        """
        usearch -usearch_local {input.fasta} -db {input.silva} -blast6out {output.silva_results} -evalue 1e-5 -threads {params.internal_threads} -strand plus -mincols 200 2> {log}
        """

rule integrate_blast_data:
    input:
        card_results = local("results/{sample}/{part}/card_results.txt"),
        silva_results = local("results/{sample}/{part}/SILVA_results.txt"),
        aro_mapping = local("data/card_db/aro_index.tsv"),
    output:
        intermed_card_results = temp("results/{sample}/{part}/intermed_card_results.csv"),
        intermed_silva_results = temp("results/{sample}/{part}/intermed_silva_results.csv"),
        integrated_data = local("results/{sample}/{part}/integrated_filtered_results.csv")
    log:
        local("logs/integrate_blast_data/{sample}_{part}.log")    
    conda:            
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/integrate_blast_data.py"

rule filter_blast_results:
    input:
        integrated_data = local("results/{sample}/{part}/integrated_filtered_results.csv")
    output:
        filtered_data = local("results/{sample}/{part}/filtered_results.csv")
    params:
        min_similarity = config["min_similarity"]
    log:
        local("logs/filter_blast_results/{sample}_{part}.log")
    conda:
        "../envs/python.yaml"   
    threads: config["max_threads"]
    script:
        "../scripts/filter_blast_results.py"

rule gzip_intermediates:
    input:
        silva_results = local("results/{sample}/{part}/SILVA_results.txt"),
        card_results = local("results/{sample}/{part}/card_results.txt"),
        int_data = local("results/{sample}/{part}/integrated_filtered_results.csv"),
        filt_data = local("results/{sample}/{part}/filtered_results.csv")
    output:
        silva_zip = local("results/{sample}/{part}/SILVA_results.txt.gz"),
        card_zip = local("results/{sample}/{part}/card_results.txt.gz"),
        integrated_data = local("results/{sample}/{part}/integrated_filtered_results.csv.gz"),     
        filt_data_zip = local("results/{sample}/{part}/filtered_results.csv.gz"),
    log:
        local("logs/gzip_blast_results/{sample}_{part}.log")
    shell:
        """
        gzip {input.silva_results} 2> {log}
        gzip {input.card_results} 2>> {log}
        gzip {input.int_data} 2>> {log}
        gzip {input.filt_data} 2>> {log}
        """