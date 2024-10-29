# Author: Yuxuan Wu

import os

configfile: "config.yaml"
print(f"configfile: config.yaml")

project_dir = config["project_dir"]
name_list = glob_wildcards(project_dir + "/{sample}" + "/raw_data/{sample}_raw.h5ad").sample
print(f"project_dir: {project_dir}")

if config['applied_files'] != 'all':
    name_list = [x for x in name_list if config['applied_files'] in x]
print(f"Pipeline applied to {name_list}")

if config["AddEmbedding.user_dataset"] is None:
    config["AddEmbedding.user_dataset"] = ''
print(f"User_dataset: {config['AddEmbedding.user_dataset']}")

DEBUG = bool(config['debug'])
def tempd(file):
    return temp(file) if not DEBUG else file

rule all:
    input:
        # finished=os.path.join(project_dir, f"integrate_input_{config['output_suffix']}", "finished.stamp"),
        merged_adata=os.path.join(project_dir, "merged.h5ad"),

rule AutoExtract:
    input:
        input_adata=os.path.join(project_dir, "{sample}", "raw_data", "{sample}_raw.h5ad"),
        pdf_file=os.path.join(project_dir, "{sample}", "raw_data", "{sample}.pdf"),
        config_file=os.path.join(project_dir, config["init_config_ini"]),
    params:
        output_dir=os.path.join(project_dir, "{sample}"),
        output_name="{sample}_" + config["output_suffix"] + "_extracted.h5ad",
    output:
        output_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_extracted.h5ad"),
        output_config_pkl=os.path.join(project_dir, "{sample}", config["config_pkl"]),
        output_log=os.path.join(project_dir, "{sample}", config["log_file"]),
    shell: """
        scExtract auto_extract \
            --adata_path {input.input_adata} \
            --pdf_path {input.pdf_file} \
            --config_path {input.config_file} \
            --output_dir {params.output_dir} \
            --output_name {params.output_name} \
            --output_config_pkl {output.output_config_pkl} \
            --output_log {output.output_log} \
            --benchmark_no_context_key no_context_annotation
    """

# Directly integrate the extracted data
rule Integrate:
    input:
        merge_output_adata=expand(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_extracted.h5ad"), sample=name_list),
        merged_embedding_dict=os.path.join(project_dir, "embedding_dict.pkl"),
        config_file=os.path.join(project_dir, config["init_config_ini"]),
    output:
        merged_adata=os.path.join(project_dir, "merged.h5ad"),
    params:
        method=config["method"],
    shell: """
        scExtract integrate \
            --file_list {input.merge_output_adata} \
            --embedding_dict_path {input.merged_embedding_dict} \
            --config_path {input.config_file} \
            --method {params.method} \
            --downsample \
            --output_path {output.merged_adata} \ 
    """

# Step-wise integration
# rule AddEmbedding:
#     input:
#         merge_output_adata=expand(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_extracted.h5ad"), sample=name_list),
#         config_file=os.path.join(project_dir, config["init_config_ini"]),
#     output:
#         merged_embedding_dict=os.path.join(project_dir, f"{config['output_suffix']}_embedding_dict.pkl"),
#     params:
#         user_dataset=config["AddEmbedding.user_dataset"],
#     shell: """
#         scExtract extract_celltype_embedding \
#             --file_list {params.user_dataset} {input.merge_output_adata} \
#             --config_path {input.config_file} \
#             --output_embedding_pkl {output.merged_embedding_dict}
#     """
#
# rule Integrate_Input:
#     input:
#         merge_output_adata=expand(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_extracted.h5ad"), sample=name_list),
#         merged_embedding_dict=os.path.join(project_dir, f"{config['output_suffix']}_embedding_dict.pkl"),
#     output:
#         merged_adata_input=directory(os.path.join(project_dir, f"integrate_input_{config['output_suffix']}")),
#         finished=os.path.join(project_dir, f"integrate_input_{config['output_suffix']}", "finished.stamp"),
#     params:
#         user_dataset=config["AddEmbedding.user_dataset"],
#     shell: """
#         mkdir -p {output.merged_adata_input}
#         mv {input.merge_output_adata} {output.merged_adata_input}
#         if [ -n "{params.user_dataset}" ]; then
#             mv {params.user_dataset} {output.merged_adata_input}
#         fi
#         mv {input.merged_embedding_dict} {output.merged_adata_input}
#         touch {output.finished}
#     """ 
