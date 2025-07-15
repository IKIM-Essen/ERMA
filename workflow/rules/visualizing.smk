# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.

rule plot_abundance_data:
    input:
        abundance_data="results/abundance/combined_genus_abundance.csv",
    output:
        report(
            "results/abundance/abundance_data.html",
            caption="../../report/abundance_data.rst",
            category="1. Combined Abundance Data",
            labels={"HTML": "Abundance data"}),
    log:
        "logs/genera_abundance_plot.log",
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/plot_abundance_data.py"


rule plot_abundance_bubble:
    input:
        abundance_data="results/abundance/combined_genus_abundance.csv",
    output:
        report(
            "results/abundance/combined_genus_abundance_bubbleplot.html",
            caption="../../report/abundance_bubble_plot.rst",
            category="1. Combined Abundance Data",
            labels={"figure": "Abundance Bubble Plot"}),
    log:
        "logs/genera_abundance_plot.log",
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/plot_abundance_bubble.py"

        
rule plot_stacked_bar_abundance:
    input:
        abundance_data="results/abundance/combined_genus_abundance.csv",
    output:
        report(
            "results/abundance/stacked_bar_abundance_plot.html",
            caption="../../report/stacked_bar_abundance_plot.rst",
            category="1. Combined Abundance Data",
            labels={"figure": "Stacked Bar Abundance Plot"}),
    log:
        "logs/stacked_bar_abundance_plot.log",
    params:
        min_abundance = config["min_abundance"]
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/plot_stacked_bar_abundance.py"


rule plot_overview_table:
    input:
        overview_table="results/qc/overview_table.txt",
    output:
        report(
            "results/qc/overview_table.html",
            caption="../../report/count_overview.rst",
            category="4. QC",
            labels={"HTML": "Count Overview Table"}),
    log:
        "logs/plot_overview_table/combined.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_overview_table.py"

rule plot_attrition:
    input:
        overview_table="results/qc/overview_table.txt",
    output:
        report(
            "results/qc/attrition_plot.png",
            caption="../../report/attrition_plot.rst",
            category="4. QC",
            labels={"Plot": "Attrition"}),
    params:
        sample_name=samples,
    log:
        "logs/plot_attrition/combined.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_attrition.py"
