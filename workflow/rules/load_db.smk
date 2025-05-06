rule get_16S_db:
    output:
        seq = local("data/silva_db/silva_seq_RNA.fasta.gz"),
    params:
        path = get_silva_db_dir(),
        seq = config["silva"]["download-path-seq"],
    conda:
        "../envs/python.yaml"
    log:
        local("logs/get_silva_db/get_silva_db.log")        
    shell:
        """
        mkdir -p {params.path}
        if [ ! -f {output.seq} ]; then \
            wget -O {output.seq} {params.seq}; \
        else \
            echo "File {output.seq} already exists, skipping download."; \
        fi
        """

rule unzip_silva_db:
    input:
        seq = local("data/silva_db/silva_seq_RNA.fasta.gz"),
    output:
        seq = local(temp("data/silva_db/silva_seq_RNA.fasta"))
    log:
        local("logs/unzip_silva_db/get_silva_db.log")
    conda:
        "../envs/python.yaml"        
    shell:
        """
        gzip -dk {input.seq} 2> {log};
        """

rule translate_silva_db:
    input:
        seq = local("data/silva_db/silva_seq_RNA.fasta")
    output:
        seq = local(temp("data/silva_db/silva_seq.fasta"))
    log:
        local("logs/translate_silva_db/translate_silva_db.log")        
    conda:
        "../envs/python.yaml"          
    shell:
        "seqtk seq -r {input.seq} > {output.seq}"

rule get_card_db:
    output:
        seq = local("data/card_db/card_seq.tar.bz2")
    params:
        path = get_card_db_dir(),
        seq = config["card"]["download-path"],
    log:
        local("logs/get_card_db/get_silva_db.log")
    conda:
        "../envs/python.yaml"        
    shell:
        """
        mkdir -p {params.path}
        if [ ! -f {output.seq} ]; then \
            wget -O {output.seq} {params.seq}; \
        else \
            echo "File {output.seq} already exists, skipping download."; \
        fi
        """

rule unzip_card_db:
    input:
        seq = local("data/card_db/card_seq.tar.bz2")
    output:
        seq = local("data/card_db/protein_fasta_protein_homolog_model.fasta"),
        aro_mapping = local("data/card_db/aro_index.tsv")
    params:
        path = get_card_db_dir()
    log:
        local("logs/unzip_card_db/unzip_card_db.log")
    conda:
        "../envs/python.yaml"        
    shell:
        """
        echo {params.path}
        tar -xvjf {input.seq} -C {params.path} 2> {log};
        """

rule makeblastdb_card:
    input:
        seq = local("data/card_db/protein_fasta_protein_homolog_model.fasta")
    output:
        db = local("data/card_db/card_db.dmnd")
    params:
        path = get_card_db_dir()
    log:
        local("logs/diamond_makedb_card/makedb_card.log")
    conda:
        "../envs/diamond.yaml"  
    shell:
        """
        diamond makedb --in {input.seq} -d {params.path}/card_db 2> {log};
        """
