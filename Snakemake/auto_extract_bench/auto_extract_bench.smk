# Author: Yuxuan Wu

import os

configfile: "config_extract.yaml"
print(f"configfile: config_extract.yaml")

project_dir = config["project_dir"]
# name_list = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
name_list = glob_wildcards(project_dir + "/{sample}" + "/raw_data/{sample}_raw.h5ad").sample
print(f"project_dir: {project_dir}")
# name_list = [x for x in name_list if 'sample16' in x]

if config['applied_files'] != 'all':
    name_list = [x for x in name_list if config['applied_files'] in x]
print(f"Pipeline applied to {name_list}")

DEBUG = bool(config['debug'])
def tempd(file):
    return temp(file) if not DEBUG else file

rule all:
    input:
        expand(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "benchmark.h5ad"), sample=name_list)

rule AutoExtract:
    input:
        input_adata=os.path.join(project_dir, "{sample}", "raw_data", "{sample}_raw.h5ad"),
        pdf_file=os.path.join(project_dir, "{sample}", "raw_data", "{sample}.pdf"),
        config_file=os.path.join(project_dir, config["init_config_ini"]),
    params:
        output_dir=os.path.join(project_dir, "{sample}"),
        output_name="{sample}_" + config["output_suffix"] + "extracted.h5ad",
    output:
        output_adata=tempd(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "extracted.h5ad")),
        output_config_pkl=os.path.join(project_dir, "{sample}", config["config_pkl"]),
        output_log=os.path.join(project_dir, "{sample}", config["log_file"]),
    shell: """
        scExtract auto_extract \
            --adata_path {input.input_adata} \
            --pdf_path {input.pdf_file} \
            --output_dir {params.output_dir} \
            --config_path {input.config_file} \
            --output_name {params.output_name} \
            --output_config_pkl {output.output_config_pkl} \
            --output_log {output.output_log} \
            --benchmark_no_context_key no_context_annotation
    """

# rule AddSingleR:
#     input:
#         output_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_extracted.h5ad"),
#     params:
#         ref_data=config["ref_data"],
#         ref_features=config["ref_features"],
#         ref_labels=config["ref_labels"],
#         singler_key=config["singler_key"],
#     output:
#         with_singler_adata=tempd(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_with_singler.h5ad")),
#     shell: """
#         scExtract add_singler_annotation \
#             --adata_path {input.output_adata} \
#             --output_path {output.with_singler_adata} \
#             --ref_data {params.ref_data} \
#             --ref_features {params.ref_features} \
#             --ref_labels {params.ref_labels} \
#             --key_added {params.singler_key} \
#             --cache_dir {wildcards.sample}_singler_cache
#     """

# rule AddCellTypist:
#     input:
#         with_singler_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_with_singler.h5ad"),
#         pdf_file=os.path.join(project_dir, "{sample}", "raw_data", "{sample}.pdf"),
#         config_file=os.path.join(project_dir, config["init_config_ini"]),
#     output:
#         with_celltypist_adata=tempd(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_with_celltypist.h5ad")),
#     shell: """
#         scExtract add_celltypist_annotation \
#             --pdf_path {input.pdf_file} \
#             --output_path {output.with_celltypist_adata} \
#             --adata_path {input.with_singler_adata} \
#             --config_path {input.config_file} \
#     """

# rule AddscType:
#     input:
#         with_celltypist_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_with_celltypist.h5ad"),
#         pdf_file=os.path.join(project_dir, "{sample}", "raw_data", "{sample}.pdf"),
#         config_file=os.path.join(project_dir, config["init_config_ini"]),
#     output:
#         with_sctype_adata=tempd(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "_with_sctype.h5ad")),
#     shell: """
#         scExtract add_sctype_annotation \
#             --pdf_path {input.pdf_file} \
#             --output_path {output.with_sctype_adata} \
#             --adata_path {input.with_celltypist_adata} \
#             --config_path {input.config_file} \
#     """

rule IntersectTrue:
    input:
        with_sctype_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "extracted.h5ad"),
        true_adata=os.path.join(project_dir, "{sample}", "raw_data" ,"{sample}_true.h5ad"),
    params:
        true_key=config["true_key"],
    output:
        with_true_adata=tempd(os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "with_true.h5ad")),
    run:
        import scanpy as sc
        adata = sc.read_h5ad(input.with_sctype_adata)
        adata_true = sc.read_h5ad(input.true_adata)
        adata = adata[adata.obs.index.isin(adata_true.obs.index)].copy()
        adata.obs[config["true_key"]] = adata_true.obs[config["true_key"]]
        adata.write(output.with_true_adata)

rule Benchmark:
    input:
        with_true_adata=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "with_true.h5ad"),
        output_config_pkl=os.path.join(project_dir, "{sample}", config["config_pkl"]),
        config_file=os.path.join(project_dir, config["init_config_ini"]),
    params:
        true_key=config["true_key"],
        method=config["method"],
        predict_group_key=config["predict_group_key"],
        similarity_key=config["similarity_key"],
    output:
        output_benchmark=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "benchmark.h5ad"),
        output_metrics=os.path.join(project_dir, "{sample}", "{sample}_" + config["output_suffix"] + "metrics.txt"),
    shell: """
        scExtract benchmark \
            --adata_path {input.with_true_adata} \
            --output_path {output.output_benchmark} \
            --config_path {input.config_file} \
            --result_metrics_path {output.output_metrics} \
            --method {params.method} \
            --true_group_key {params.true_key} \
            --predict_group_key {params.predict_group_key} \
            --similarity_key {params.similarity_key} \
            --output_config_pkl {input.output_config_pkl}
    """