# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule generate_percidt_genus:
    input:
        filtered_data=local(
            expand(
                "results/{{sample}}/{part}/filtered_results.csv.gz",
                part=get_numpart_list(),
            )
        ),
    output:
        report(
            local("results/boxplots/{sample}_idt_genus_plot.png"),
            caption="../../report/identity_read_count_per_genus.rst",
            category="2. Samplewise Genus per Percentage Identity",
            subcategory="{sample}",
            labels={"sample": "{sample}", "Plot": "Identity/Read Count per Genus"},
        ),
    log:
        local("logs/generate_percidt_genus/{sample}.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/boxplot_percidt_per_genus.py"


rule plot_alignment_length_boxplot:
    input:
        filtered_data=local(
            expand(
                "results/{sample}/{part}/filtered_results.csv.gz",
                sample=samples,
                part=get_numpart_list(),
            )
        ),
    output:
        report(
            local("results/boxplots/combined_allength_boxplot.png"),
            caption="../../report/boxplot.rst",
            category="3. Boxplots",
            labels={"boxplot": "Alignment Length"},
        ),
    params:
        sample_name=samples,
    log:
        local("logs/plot_alignment_length_boxplot/combined.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/boxplot_align_lengths.py"


rule plot_percentage_identity_boxplot:
    input:
        filtered_data=local(
            expand(
                "results/{sample}/{part}/filtered_results.csv.gz",
                sample=samples,
                part=get_numpart_list(),
            )
        ),
    output:
        report(
            local("results/boxplots/combined_percidt_boxplot.png"),
            caption="../../report/boxplot.rst",
            category="3. Boxplots",
            labels={"boxplot": "Percentage Identity"},
        ),
    params:
        sample_name=samples,
    log:
        local("logs/plot_percentage_identity_boxplot/combined.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/boxplot_percidt.py"


rule plot_evalue_boxplot:
    input:
        filtered_data=local(
            expand(
                "results/{sample}/{part}/filtered_results.csv.gz",
                sample=samples,
                part=get_numpart_list(),
            )
        ),
    output:
        report(
            local("results/boxplots/combined_evalue_boxplot.png"),
            caption="../../report/boxplot.rst",
            category="3. Boxplots",
            labels={"boxplot": "E-Value"},
        ),
    params:
        sample_name=samples,
    log:
        local("logs/plot_evalue_boxplot/combined.log"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/boxplot_evalue.py"
