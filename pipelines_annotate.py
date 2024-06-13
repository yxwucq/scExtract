import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from auto_extract.auto_extract import auto_extract
from benchmark.benchmark import benchmark_annotation
from methods_comparison.singler_anno import add_singler_annotation

import scanpy as sc
import os

sample = sys.argv[1]
print(f"Processing {sample}...")

raw_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/raw_data/"
processed_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/processed_data/"

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
