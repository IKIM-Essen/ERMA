rule get_16S_db:
    output:
        seq = local("data/silva_db/silva_seq_RNA.fasta.gz"),
    params:
        seq = config["silva"]["download-path-seq"],
        path = "data/silva_db"
    log:
        "logs/get_16S_db/log.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        mkdir -p {params.path};
        cd {params.path};
        wget -O silva_seq_RNA.fasta.gz {params.seq} 2> {log};
        """

rule unzip_silva_db:
    input:
        seq = local("data/silva_db/silva_seq_RNA.fasta.gz"),
    output:
        seq = local(temp("data/silva_db/silva_seq_RNA.fasta"))
    log:
        "logs/unzip_silva_db/log.log"
    shell:
        """
        gzip -dk {input.seq} 2> {log};
        """

rule translate_silva_db:
    input:
        seq = local("data/silva_db/silva_seq_RNA.fasta")
    output:
        seq = local(temp("data/silva_db/silva_seq.fasta"))
    conda:
        "../envs/python.yaml"          
    shell:
        "seqtk seq -r {input.seq} > {output.seq}"

rule get_card_db:
    output:
        seq = local("data/card_db/card_seq.tar.bz2")
    params:
        seq = config["card"]["download-path"],
        path = "data/card_db"
    log:
        "logs/get_card_db/log.log"
    shell:
        """
        mkdir -p {params.path};
        cd {params.path};
        wget -O card_seq.tar.bz2 {params.seq} 2> {log};
        """

rule unzip_card_db:
    input:
        seq = local("data/card_db/card_seq.tar.bz2")
    output:
        seq = local("data/card_db/protein_fasta_protein_homolog_model.fasta"),
        aro_mapping = local("data/card_db/aro_index.tsv")
    params:
        path = "data/card_db"
    log:
        "logs/unzip_card_db/log.log"
    shell:
        """
        tar -xvjf {input.seq} -C {params.path} 2> {log};
        """

rule makeblastdb_card:
    input:
        seq = local("data/card_db/protein_fasta_protein_homolog_model.fasta")
    output:
        db = local("data/blast_db/card_db.pdb")
    params:
        path = "data/blast_db/card_db"
    log:
        "logs/makeblastdb_card/log.log"
    conda:
        "../envs/blast.yaml"  
    shell:
        """
        makeblastdb -in {input.seq} -dbtype prot -out {params.path} 2> {log};
        """
