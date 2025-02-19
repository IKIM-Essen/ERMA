samples = [os.path.basename(f).replace(".fastq.gz", "") for f in glob.glob(os.path.join(config["base_dir"], "data", "fastq", "*.fastq.gz"))]

rule fastqc:
    input:
        "{base_dir}/data/fastq/{sample}.fastq.gz"
    output:
        html="{base_dir}/results/fastqc/{sample}.html",
        zip="{base_dir}/results/fastqc/{sample}_fastqc.zip",
    log:
        "{base_dir}/logs/input/{sample}.log",
    threads: 8
    wrapper:
        "v1.3.1/bio/fastqc"


rule multiqc_report:
    input:
        expand(
            "{base_dir}/results/fastqc/{sample}_fastqc.zip",
            sample=samples,base_dir=config["base_dir"]
        ),
    output:
        report(
            "{base_dir}/results/qc/multiqc.html",
            caption="../../report/genus_top_hits.rst",
            htmlindex="multiqc.html",
            category="4. QC",
        ),
    log:
        "{base_dir}/logs/multiqc/multiqc.log",
    wrapper:
        "v3.3.3/bio/multiqc"