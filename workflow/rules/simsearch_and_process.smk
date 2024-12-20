rule blast_card:
    input:
        fasta = "{base_dir}/data/input/{sample}.part_{part}.fasta",
        card = "{base_dir}/data/blast_db/card_db.pdb"
    output:
        card_results = "{base_dir}/results/{sample}/{part}/card_results.txt"
    params:
        db = "{base_dir}/data/blast_db",
        internal_threads = config["num_parts"],        
    log:
        "{base_dir}/logs/blast_card/{sample}_{part}.log"    
    conda:
        "../envs/blast.yaml" 
    threads: config["max_threads"]      
    shell:
        """
        blastx -query {input.fasta} -db {params.db}/card_db -out {output.card_results} -outfmt 6 -evalue 1e-5 -num_threads {params.internal_threads} 2> {log}
        """

rule usearch_silva:
    input:
        fasta = "{base_dir}/data/input/{sample}.part_{part}.fasta",
        silva = "{base_dir}/data/silva_db/silva_seq.fasta"
    output:
        silva_results = "{base_dir}/results/{sample}/{part}/SILVA_results.txt"
    log:
        "{base_dir}/logs/blast_silva/{sample}_{part}.log" 
    params:
        internal_threads = config["max_threads"],
        min_id = config["min_similarity"]
    conda:
        "../envs/blast.yaml"
    threads: config["max_threads"]                  
    shell:
        """
        usearch -usearch_local {input.fasta} -db {input.silva} -blast6out {output.silva_results} -evalue 1e-5 -threads {params.internal_threads} -strand plus mincols 200 2> {log}
        """

rule gzip_blast_results:
    input:
        silva_results = "{base_dir}/results/{sample}/{part}/SILVA_results.txt",
        card_results = "{base_dir}/results/{sample}/{part}/card_results.txt"
    output:
        silva_zip = "{base_dir}/results/{sample}/{part}/SILVA_results.txt.gz",
        card_zip = "{base_dir}/results/{sample}/{part}/card_results.txt.gz"
    log:
        "{base_dir}/logs/gzip_blast_results/{sample}_{part}.log"        
    shell:
        """
        gzip {input.silva_results} 2> {log}
        gzip {input.card_results} 2>> {log}
        """

rule integrate_blast_data:
    input:
        card_results = "{base_dir}/results/{sample}/{part}/card_results.txt.gz",
        silva_results = "{base_dir}/results/{sample}/{part}/SILVA_results.txt.gz",
        aro_mapping = "{base_dir}/data/card_db/aro_index.tsv",    
    output:
        intermed_card_results = temp("{base_dir}/results/{sample}/{part}/intermed_card_results.csv"),
        intermed_silva_results = temp("{base_dir}/results/{sample}/{part}/intermed_silva_results.csv"),
        integrated_data = "{base_dir}/results/{sample}/{part}/integrated_filtered_results.csv"
    params:
        chunksize = config["chunksize"]
    log:
        "{base_dir}/logs/integrate_blast_data/{sample}_{part}.log"    
    conda:            
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/integrate_blast_data.py"

rule filter_blast_results:
    input:
        integrated_data = "{base_dir}/results/{sample}/{part}/integrated_filtered_results.csv"
    output:
        filtered_data = "{base_dir}/results/{sample}/{part}/filtered_results.csv"
    params:
        min_similarity = config["min_similarity"]
    log:
        "{base_dir}/logs/filter_blast_results/{sample}_{part}.log"         
    conda:
        "../envs/python.yaml"   
    threads: config["max_threads"]
    script:
        "../scripts/filter_blast_results.py"

rule gzip_filtered_blast_data:
    input:
        int_data = "{base_dir}/results/{sample}/{part}/filtered_results.csv"
    log:
        "{base_dir}/logs/gzip_integrated_blast_data/{sample}_{part}.log"    
    output:
        int_data_zip = "{base_dir}/results/{sample}/{part}/filtered_results.csv.gz"
    shell:
        "gzip -9 {input.int_data} 2> {log}"