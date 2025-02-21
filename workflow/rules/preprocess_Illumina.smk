rule merge_fastq:
    input:
        r1="{base_dir}/data/fastq/{sample}_R1_001.fastq.gz",
        r2="{base_dir}/data/fastq/{sample}_R2_001.fastq.gz"
    output:
        out = temp("{base_dir}/data/fastq/temp/{sample}.fastq.gz")
    log:
        "{base_dir}/logs/merge_fastq/{sample}.log"
    conda:
        "../envs/python.yaml"                
    shell:
        "NGmerge -1 {input.r1} -2 {input.r2} -o {output.out}"

rule convert_fastq_to_fasta:
    input:
        "{base_dir}/data/fastq/temp/{sample}.fastq.gz"
    output:
        temp("{base_dir}/data/fastq/temp/{sample}.fasta"),
    log:
        "{base_dir}/logs/convert_fastq_to_fasta/{sample}.log"    
    conda:
        "../envs/python.yaml"        
    shell:
        """
        seqtk seq -a {input} > {output} 2> {log}
        """

rule split_fasta_file:
    input:
        "{base_dir}/data/fastq/temp/{sample}.fasta",
    output:
        temp("{base_dir}/data/fastq/{sample}.part_{part}.fasta")
    params:
        outdir = config["base_dir"],
        num_parts = config["num_parts"],    
    shell:
        """
        seqkit split2 --by-part {params.num_parts} {input} --out-dir {params.outdir}/data/fastq/
        """