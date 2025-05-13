samples = [os.path.basename(f).replace(".fastq.gz", "") for f in glob.glob(os.path.join(config["base_dir"], config["fastq_dir"], "*.fastq.gz"))]

rule fastqc:
    input:
        lambda wildcards: local(os.path.join(config["fastq_dir"], f"{wildcards.sample}.fastq.gz"))
    output:
        html=local("results/fastqc/{sample}.html"),
        zip=local("results/fastqc/{sample}_fastqc.zip"),
    log:
        local("logs/input/{sample}.log"),
    threads: 8
    resources: mem_mb = 1024
    wrapper:
        "v6.0.1/bio/fastqc"

rule multiqc_report:
    input:
        local(expand(
            "results/fastqc/{sample}_fastqc.zip",
            sample=samples
        )),
    output:
        report(
            local("results/qc/multiqc.html"),
            caption="../../report/multiqc.rst",
            category="4. QC",
            labels={
                "File": "MultiQC Report"
            }     
        ),
    log:
        local("logs/multiqc/multiqc.log"),
    wrapper:
        "v6.0.1/bio/multiqc"  

rule merge_overview_per_sample:
    input:
        checkpoint = local(expand("results/{sample}/{part}/checkpoint.txt",sample=samples,part=get_numpart_list())),
        overview_tables = local(expand("results/{{sample}}/{part}/overview_table.txt",part=get_numpart_list()))
    output:
        report(
            local("results/{sample}/overview_table.html"),
            caption = "../../report/count_overview_per_sample.rst",
            category="2. Single Sample Abundance Data",
            subcategory="{sample}",
            labels={
                "sample": "{sample}",
                "table": "Count Overview"
            }
        )             
    params:
        sample_name = "{sample}",
    log:
        local("logs/merge_overview/{sample}/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/merge_overview_smpl.py"

rule merge_overview_to_one:
    input:
        checkpoint = local(expand("results/{sample}/{part}/checkpoint.txt",sample=samples,part=get_numpart_list())),
        overview_tables = local(expand("results/{sample}/{part}/overview_table.txt",sample=samples,part=get_numpart_list()))
    output:
        local("results/qc/overview_table.txt"),
    params:
        sample_name = samples,
    log:
        local("logs/merge_overview/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/merge_overview_all.py"

rule plot_overview:
    input:
        overview_table = local("results/qc/overview_table.txt")
    output:
        report(
            local("results/qc/overview_plot.png"),
            caption = "../../report/count_overview.rst",
            category="4. QC",
            labels={
                "File": "Count Overview"
            }       
        )
    params:
        sample_name = samples,
    log:
        local("logs/plot_overview/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/plot_overview.py"