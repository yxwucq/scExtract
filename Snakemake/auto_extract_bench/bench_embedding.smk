# Author: Yuxuan Wu

import os

configfile: "config.yaml"
print(f"configfile: config.yaml")

project_dir = config["project_dir"]
# name_list = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
name_list = glob_wildcards(project_dir + "/{sample}" + "/raw_data/{sample}_raw.h5ad").sample
print(f"project_dir: {project_dir}")

if config['applied_files'] != 'all':
    name_list = [x for x in name_list if config['applied_files'] in x]
print(f"Pipeline applied to {name_list}")

DEBUG = bool(config['debug'])
def tempd(file):
    return temp(file) if not DEBUG else file

rule all:
    input:
        expand(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_metrics.txt"), sample=name_list)

rule Benchmark:
    input:
        output_benchmark=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_benchmark.h5ad"),
        output_config_pkl=os.path.join(project_dir, "{sample}", config["config_pkl"]),
        config_file=os.path.join(project_dir, config["init_config_ini"]),
    params:
        true_key=config["true_key"],
        method=config["method"],
        predict_group_key=config["predict_group_key"],
        similarity_key=config["similarity_key"],
    output:
        output_benchmark=temp(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_benchmark_embedding.h5ad")),
        output_metrics=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_metrics.txt"),
    shell: """
        scExtract benchmark \
            --adata_path {input.output_benchmark} \
            --output_path {output.output_benchmark} \
            --config_path {input.config_file} \
            --result_metrics_path {output.output_metrics} \
            --method embedding \
            --true_group_key {params.true_key} \
            --predict_group_key {params.predict_group_key} \
            --similarity_key {params.similarity_key} \
            --output_config_pkl {input.output_config_pkl}
    """