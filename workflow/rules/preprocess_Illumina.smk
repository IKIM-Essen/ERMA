rule merge_fastq:
    input:
        r1 = local("data/fastq/{sample}_R1_001.fastq.gz"),
        r2 = local("data/fastq/{sample}_R2_001.fastq.gz")
    output:
        out = local("data/fastq/{sample}.fastq.gz")
    log:
        local("logs/merge_fastq/{sample}.log")
    conda:
        "../envs/python.yaml"                
    shell:
        "NGmerge -1 {input.r1} -2 {input.r2} -o {output.out}"

rule convert_fastq_to_fasta:
    input:
        local("data/fastq/{sample}.fastq.gz")
    output:
        local(temp("data/fastq/{sample}.fasta")),
    log:
        local("logs/convert_fastq_to_fasta/{sample}.log")
    conda:
        "../envs/python.yaml"        
    shell:
        """
        seqtk seq -a {input} > {output} 2> {log}
        """

rule split_fasta_file:
    input:
        local("data/fastq/{sample}.fasta")
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