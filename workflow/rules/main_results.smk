runname = "".join(config["runname"])
seq_tech = "".join(config["seq_tech"])
        
rule genera_abundance_table:
    input:
        filtered_data = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv.gz",sample=samples,part=get_numpart_list()),
    output:
        report(
            "{base_dir}/results/abundance/combined_genus_abundance.csv",
            caption = "../../report/genus_top_hits.rst",
            category="3. General data"
        )
    params:
        sample_name = samples,
    log:
        "{base_dir}/logs/genera_abundance_table.log"                  
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_table.py"

rule abundance_bubble_plot:
    input:
        abundance_data = "{base_dir}/results/abundance/combined_genus_abundance.csv",
    output:
        report(
            "{base_dir}/results/abundance/combined_genus_abundance_bubbleplot.html",
            caption = "../../report/genus_top_hits.rst",
            category="1. Main result"
        ),
        report(
            "{base_dir}/results/abundance/reads_per_found_AMR.csv",
            caption = "../../report/genus_top_hits.rst",
            category="1. Main result"
        ),
        "{base_dir}/results/abundance/reads_per_found_AMR_per_sample.csv",
    params:
        abundance_filter = 0.01
    log:
        "{base_dir}/logs/genera_abundance_plot.log"                  
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_plot.py"

rule found_AMR_per_sample:
    input:
        abundance_data = "{base_dir}/results/abundance/reads_per_found_AMR_per_sample.csv",
    output:
        report(
            "{base_dir}/results/abundance/AMR_per_sample.png",
            caption = "../../report/genus_top_hits.rst",
            category="1. Main result"
        ),
    log:
        "{base_dir}/logs/genera_abundance_plot.log"                  
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/AMR_per_sample.py"