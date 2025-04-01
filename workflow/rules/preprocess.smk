rule convert_fastq_to_fasta:
    input:
        local("data/fastq/{sample}.fastq.gz")
    output:
        local(temp("results/fastq/{sample}.fasta"))
    log:
        local("logs/convert_fastq_to_fasta/{sample}.log")
    conda:
        "../envs/python.yaml"        
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"

rule split_fasta_file:
    input:
        local("results/fastq/{sample}.fasta")
    output:
        local(temp("results/fastq/{sample}.part_{part}.fasta"))
    log:
        local("logs/split_fasta/{sample}_{part}.log")
    conda:
        "../envs/python.yaml"        
    params:
        outdir = config["base_dir"],
        num_parts = config["num_parts"],   
    shell:
        """
        seqkit split2 --by-part {params.num_parts} {input} --out-dir {params.outdir}/results/fastq/
        """