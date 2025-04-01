samples = [os.path.basename(f).replace(".fastq.gz", "") for f in glob.glob(os.path.join(config["base_dir"], "data", "fastq", "*.fastq.gz"))]

rule fastqc:
    input:
        local("data/fastq/{sample}.fastq.gz")
    output:
        html=local("results/fastqc/{sample}.html"),
        zip=local("results/fastqc/{sample}_fastqc.zip"),
    log:
        local("logs/input/{sample}.log"),
    threads: 8
    resources: mem_mb = 1024
    wrapper:
        "v5.8.3/bio/fastqc"

rule multiqc_report:
    input:
        local(expand(
            "results/fastqc/{sample}_fastqc.zip",
            sample=samples
        )),
    output:
        report(
            local("results/qc/multiqc.html"),
            caption="../../report/genus_top_hits.rst",
            htmlindex="multiqc.html",
            category="4. QC",
        ),
    log:
        local("logs/multiqc/multiqc.log"),
    wrapper:
        "v5.8.3/bio/multiqc"