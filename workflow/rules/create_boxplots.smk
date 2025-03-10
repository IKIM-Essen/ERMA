rule generate_percidt_genus:
    input:
        filtered_data = local(expand("results/{{sample}}/{part}/filtered_results.csv.gz",part=get_numpart_list()))
    output:
        report(
            local("results/{sample}/genus_idt_per_genus_plot.png"),
            caption = "../../report/genus_top_hits.rst",
            category="2. Genus percentage Identity",
        )
    log:
        local("logs/generate_percidt_genus/{sample}.log")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/percidt_per_genus.py"

rule plot_alignment_length_boxplot:
    input:
        csv_files = local(expand("results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list()))
    output:
        report(
            local("results/boxplots/combined_allength_boxplot.png"),
            caption = "../../report/genus_top_hits.rst",
            category="3. Statistics",        
        )        
    params:
        sample_name = samples,
    log:
        local("logs/plot_alignment_length_boxplot/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/align_lengths_boxplots.py"

rule plot_percentage_identity_boxplot:
    input:
        csv_files = local(expand("results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list()))
    output:
        report(
            local("results/boxplots/combined_percidt_boxplot.png"),
            caption = "../../report/genus_top_hits.rst",
            category="3. Statistics",        
        )
    params:
        sample_name = samples,
    log:
        local("logs/plot_percentage_identity_boxplot/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/percidt_boxplots.py"

rule plot_evalue_boxplot:
    input:
        csv_files = local(expand("results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list()))
    output:
        report(
            local("results/boxplots/combined_evalue_boxplot.png"),
            caption = "../../report/genus_top_hits.rst",
            category="3. Statistics"
        )
    params:
        sample_name = samples,
    log:
        local("logs/plot_evalue_boxplot/combined.log")
    conda:
        "../envs/python.yaml"     
    script:
        "../scripts/evalue_boxplots.py"