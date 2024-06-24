import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from auto_extract.auto_extract import auto_extract
from benchmark.benchmark import benchmark_annotation
from methods_comparison.singler_anno import add_singler_annotation

import scanpy as sc
import os

sample = sys.argv[1]
print(f"Processing {sample}...")

raw_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/benchmark_datasets/{sample}/raw_data/"
processed_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/benchmark_datasets/{sample}/processed_data/"

## 01.annotation
print("Processing auto_extract...")
auto_extract(
    adata_path=os.path.join(raw_data_dir, f"{sample}_raw.h5ad"),
    pdf_path=os.path.join(raw_data_dir, f"{sample}.pdf"),
    output_dir=processed_data_dir,
    output_name=f"{sample}_processed.h5ad",
    benchmark_no_context_key='no_context_annotation'
)

## 02.add singler annotation
print("Processing singler...")
add_singler_annotation(adata_path=os.path.join(processed_data_dir, f"{sample}_processed.h5ad"))

## 03.add true & benchmark
print("Processing benchmark...")
adata = sc.read_h5ad(os.path.join(processed_data_dir, f"{sample}_processed.h5ad"))
adata_true = sc.read_h5ad(os.path.join(processed_data_dir, f"{sample}_true.h5ad"))
adata = adata[adata.obs.index.isin(adata_true.obs.index)].copy()
adata.obs['cell_type'] = adata_true.obs['cell_type']
adata.write(os.path.join(processed_data_dir, f"{sample}_processed.h5ad"))

benchmark_annotation(
    adata_path=os.path.join(processed_data_dir, f"{sample}_processed.h5ad"),
    output_path=os.path.join(processed_data_dir, f"{sample}_processed_benchmark.h5ad"),
    true_group_key='cell_type',
    config_path=os.path.join(processed_data_dir, f"config.pkl"),
    result_metrics_path=os.path.join(processed_data_dir, f"{sample}_benchmark_metrics.txt"),
    method = "embedding",
    predict_group_key='scextract,no_context_annotation,singler_annotation',
    similarity_key='similarity_scextract,similarity_no_context_annotation,similarity_singler'
)