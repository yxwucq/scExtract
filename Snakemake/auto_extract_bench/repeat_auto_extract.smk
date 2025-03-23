# Author: Yuxuan Wu
import os

repeat_times = 6
repeats = list(map(str, range(1, repeat_times+1)))
sample_name = glob_wildcards("raw_data/{sample}_raw.h5ad").sample

rule all:
    input:
        expand("{sample}_repeat_{repeat}.h5ad", repeat=repeats, sample=sample_name)

rule AutoExtract:
    input:
        input_adata=os.path.join("raw_data", "{sample}_raw.h5ad"),
        pdf_file=os.path.join("raw_data", "{sample}.pdf"),
        config_file='../repeat_auto_extract_config.ini',
    params:
        output_name="{sample}_repeat_{repeat}.h5ad",
    output:
        output_adata="{sample}_repeat_{repeat}.h5ad",
        output_config_pkl="{sample}_repeat_{repeat}_config.pkl",
        output_log="{sample}_repeat_{repeat}.log",
    shell: """
        scExtract auto_extract \
            --adata_path {input.input_adata} \
            --pdf_path {input.pdf_file} \
            --output_dir . \
            --config_path {input.config_file} \
            --output_name {params.output_name} \
            --output_config_pkl {output.output_config_pkl} \
            --output_log {output.output_log} \
            --benchmark_no_context_key no_context_annotation
    """