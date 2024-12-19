seq_tech = "".join(config["seq_tech"])

rule merge_fastq:
    input:
        r1="{base_dir}/data/input/{sample}_R1_001.fastq.gz",
        r2="{base_dir}/data/input/{sample}_R2_001.fastq.gz"
    output:
        out = "{base_dir}/data/input/temp/{sample}.fastq.gz"
    log:
        "{base_dir}/logs/merge_fastq/{sample}.log"
    conda:
        "../envs/python.yaml"                
    shell:
        "NGmerge -1 {input.r1} -2 {input.r2} -o {output.out}"

rule decompress_fastq:
    input:
        "{base_dir}/data/input/temp/{sample}.fastq.gz"
    output:
        "{base_dir}/data/input/temp/{sample}.fastq"
    log:
        "{base_dir}/logs/decompress_fastq/{sample}.log"    
    shell:
        "gzip -dk {input} 2> {log}"

rule convert_fastq_to_fasta:
    input:
        "{base_dir}/data/input/temp/{sample}.fastq"
    output:
        "{base_dir}/data/input/temp/{sample}.fasta"
    log:
        "{base_dir}/logs/convert_fastq_to_fasta/{sample}.log"    
    conda:
        "../envs/python.yaml"        
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"

rule split_fasta_file:
    input:
        "{base_dir}/data/input/temp/{sample}.fasta"
    output:
        "{base_dir}/data/input/{sample}.part_{part}.fasta"
    params:
        outdir = config["base_dir"],
        num_parts = config["num_parts"],    
    shell:
        """
        seqkit split2 --by-part {params.num_parts} {input} --out-dir {params.outdir}/data/input/
        """