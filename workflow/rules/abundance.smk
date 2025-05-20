# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule single_genera_abundance_table:
    input:
        filtered_data=local(
            expand(
                "results/{{sample}}/{part}/filtered_results.csv.gz",
                part=get_numpart_list(),
            )
        ),
    output:
        report(
            local("results/{sample}/genus_abundance.html"),
            caption="../../report/genus_abundance_table.rst",
            category="2. Single Sample Abundance Data",
            subcategory="{sample}",
            labels={"sample": "{sample}", "table": "Genera Abundance"},
        ),
    params:
        sample_name="{sample}",
        parts=get_numpart_list(),
    log:
        local("logs/{sample}/genera_abundance_table.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/single_genera_abundance_table.py"


rule combined_genera_abundance_table:
    input:
        filtered_data=local(
            expand(
                "results/{sample}/{part}/filtered_results.csv.gz",
                sample=samples,
                part=get_numpart_list(),
            )
        ),
    output:
        csv=local("results/abundance/combined_genus_abundance.csv"),
    params:
        sample_name=samples,
    log:
        local("logs/genera_abundance_table.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/combined_genera_abundance_table.py"


rule abundance_bubble_plot:
    input:
        abundance_data=local("results/abundance/combined_genus_abundance.csv"),
    output:
        report(
            local("results/abundance/combined_genus_abundance_bubbleplot.html"),
            caption="../../report/abundance_bubble_plot.rst",
            category="1. Combined Abundance Data",
            labels={"figure": "Abundance Bubble Plot"},
        ),
    log:
        local("logs/genera_abundance_plot.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/combined_genera_abundance_plot.py"


rule reads_per_AMR:
    input:
        abundance_data=local("results/abundance/combined_genus_abundance.csv"),
    output:
        report(
            local("results/abundance/reads_per_found_AMR.html"),
            caption="../../report/reads_per_AMR.rst",
            category="1. Combined Abundance Data",
            labels={"table": "Reads per AMR"},
        ),
    log:
        local("logs/genera_abundance_plot.log"),
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/combined_reads_per_amr.py"
