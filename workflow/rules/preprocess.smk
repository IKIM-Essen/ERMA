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
        seqkit split2 --by-part {params.num_parts} -w 0 {input} --out-dir {params.outdir}/results/fastq/
        """

if config["seq_tech"] == "Illumina":
    rule prepare_fastqs:
        output:
            sentinel=local(temp("data/fastq/merged_done.txt")),
        conda:
            "../envs/python.yaml"    
        script:
            "../scripts/prepare_fastq.sh"

elif config["seq_tech"] == "ONT":
    rule prepare_fastqs:
        output:
            sentinel=local(temp("data/fastq/merged_done.txt")),
        params:
            run_path = config["ONT"]["fastq_pass_path"],
            sample_name_path = config["ONT"]["sample_name_path"],
            target_fragment_length = config["ONT"]["target_fragment_length"],
            filter_intervall = config["ONT"]["filter_intervall"],
            output_dir = "data/fastq/",   
        conda:
            "../envs/python.yaml"    
        script:
            "../scripts/prepare_fastq.py"