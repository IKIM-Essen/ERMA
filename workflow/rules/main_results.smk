rule genera_abundance_table:
    input:
        filtered_data = local("results/{sample}/{part}/filtered_results.csv.gz"),
    output:
        report(
            local(temp("results/{sample}/{part}/genus_abundance.html")),
            caption = "../../report/genus_abundance_table.rst",
            category="2. Single Sample Abundance Data",
            subcategory="{sample}",
            labels={
                "sample": "{sample}",
                "table":"Genera Abundance"}
        ),
    params:
        sample_name = "{sample}"
    log:
        local("logs/{sample}/{part}/genera_abundance_table.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_table.py"

rule combined_genera_abundance_table:
    input:
        filtered_data = local(expand("results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list())),
    output:
        csv = local("results/abundance/combined_genus_abundance.csv"),
    params:
        sample_name = samples,
    log:
        local("logs/genera_abundance_table.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/combined_genera_abundance_table.py"

rule abundance_bubble_plot:
    input:
        abundance_data = local("results/abundance/combined_genus_abundance.csv"),
    output:
        report(
            local("results/abundance/combined_genus_abundance_bubbleplot.html"),
            caption = "../../report/abundance_bubble_plot.rst",
            category="1. Combined Abundance Data",
            labels={"figure":"Abundance Bubble Plot"}
        ),    
    log:
        local("logs/genera_abundance_plot.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_plot.py"

rule reads_per_AMR:
    input:
        abundance_data = local("results/abundance/combined_genus_abundance.csv"),
    output:
        report(
            local("results/abundance/reads_per_found_AMR.html"),
            caption = "../../report/reads_per_AMR.rst",
            category="1. Combined Abundance Data",
            labels={"table":"Reads per AMR"}
        ),    
    log:
        local("logs/genera_abundance_plot.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/reads_per_amr.py"

rule plot_overview:
    input:
        overview_table = local("results/qc/overview_table.txt")
    output:
        report(
            local("results/qc/overview_plot.png"),
            caption = "../../report/boxplot.rst",
            category="4. QC",
            labels={
                "File": "Overview"
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