rule merge_fastq:
    input:
        local(r1="{base_dir}/data/fastq/{sample}_R1_001.fastq.gz"),
        local(r2="{base_dir}/data/fastq/{sample}_R2_001.fastq.gz")
    output:
        local(out = "{base_dir}/data/fastq/temp/{sample}.fastq.gz")
    log:
        local("{base_dir}/logs/merge_fastq/{sample}.log")
    conda:
        "../envs/python.yaml"                
    shell:
        "NGmerge -1 {input.r1} -2 {input.r2} -o {output.out}"

rule convert_fastq_to_fasta:
    input:
        local("{base_dir}/data/fastq/temp/{sample}.fastq.gz")
    output:
        local(temp("{base_dir}/data/fastq/temp/{sample}.fasta")),
    log:
        local("{base_dir}/logs/convert_fastq_to_fasta/{sample}.log")
    conda:
        "../envs/python.yaml"        
    shell:
        """
        seqtk seq -a {input} > {output} 2> {log}
        """

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