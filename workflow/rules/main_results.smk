runname = "".join(config["runname"])
seq_tech = "".join(config["seq_tech"])
        
rule genera_abundance_table:
    input:
        filtered_data = local(expand("results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list())),
    output:
        report(
            local("results/abundance/combined_genus_abundance.html"),
            caption = "../../report/genus_top_hits.rst",
            category="1. Abundance"
        ),
        csv = local("results/abundance/combined_genus_abundance.csv"),
    params:
        sample_name = samples,
    log:
        local("logs/genera_abundance_table.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_table.py"

rule abundance_bubble_plot:
    input:
        abundance_data = local("results/abundance/combined_genus_abundance.csv"),
    output:
        report(
            local("results/abundance/combined_genus_abundance_bubbleplot.html"),
            caption = "../../report/genus_top_hits.rst",
            category="1. Abundance"
        ),
        report(
            local("results/abundance/reads_per_found_AMR.html"),
            caption = "../../report/genus_top_hits.rst",
            category="1. Abundance"
        ),
    params:
        abundance_filter = 0.001
    log:
        local("logs/genera_abundance_plot.log")
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_plot.py"