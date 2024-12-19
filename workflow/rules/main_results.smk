runname = "".join(config["runname"])
seq_tech = "".join(config["seq_tech"])

if seq_tech == "Illumina":
    samples = list(set([re.split(r'_R\d_',r)[0] for r in reads]))
elif seq_tech == "ONT":
    samples = reads
else:
    raise ValueError("Invalid sequencing technology specified. Check config file and README.")

rule generate_genus_distribution_plot:
    input:
        filtered_data = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv", 
                                sample=samples,part=get_numpart_list())
    output:
        report(
            "{base_dir}/results/{sample}/genus_distribution_plot.html",
            caption = "../../report/genus_top_hits.rst",
            htmlindex="genus_distribution_plot.html",
            category="1. Genus distribution/{sample}"
        ),
        "{base_dir}/results/{sample}/combined_df.csv",
    params:
        sample_name = samples,
    log:
        "{base_dir}/logs/generate_genus_distribution_plot/{sample}.log"         
    conda:
        "../envs/python.yaml"    
    threads: config["max_threads"]              
    script:
        "../scripts/generate_genus_distribution_plot.py"
        
rule genera_abundance_table:
    input:
        filtered_data = expand("{{base_dir}}/results/{sample}/{part}/filtered_results.csv",sample=samples,part=get_numpart_list())
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
            category="0. Main result"
        ),
        report(
            "{base_dir}/results/abundance/reads_per_found_AMR.csv",
            caption = "../../report/genus_top_hits.rst",
            category="0. Main result"
        ),
    params:
        abundance_filter = 0.001
    log:
        "{base_dir}/logs/genera_abundance_plot.log"                  
    conda:
        "../envs/python.yaml"
    threads: config["max_threads"]
    script:
        "../scripts/genera_abundance_plot.py"