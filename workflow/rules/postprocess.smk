# Copyright 2024 Adrian Dörr.
# Licensed under the MIT License (https://opensource.org/license/mit)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule tar_single_sample_dirs:
    input:    
        local("results/abundance/stacked_bar_abundance_plot.html"),
        local("results/abundance/combined_genus_abundance_bubbleplot.html"),
        local("results/abundance/abundance_data.html"),
        local("results/boxplots/combined_allength_boxplot.png"),
        local("results/boxplots/combined_evalue_boxplot.png"),
        local("results/boxplots/combined_percidt_boxplot.png"),
        local("results/qc/multiqc.html"),
        local("results/qc/attrition_plot.png"),
        local("results/qc/overview_table.html"),
    output:
        local("results/single_sample_similarity_search_data.tar.gz")
    log:
        local("logs/tar/tar.log")
    script:
        "../scripts/postprocess_tar_intermediates.py"        


rule concatenate_logs:
    input:
        local("results/single_sample_similarity_search_data.tar.gz")
    output:
        local("logs/logs.json")
    script:
        "../scripts/postprocess_concat_logs.py"
