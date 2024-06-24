import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from auto_extract.auto_extract import auto_extract
from benchmark.benchmark import benchmark_annotation
from methods_comparison.singler_anno import add_singler_annotation

import scanpy as sc
import os

failed_samples = []

sample_list = os.listdir('/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/')
print(f"Total samples: {len(sample_list)}")

sample_list = [sample for sample in sample_list if os.path.exists(os.path.join(f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/raw_data/", f"{sample}_raw.h5ad"))]
print(f"Samples with raw data: {len(sample_list)}")

sample_list = [sample for sample in sample_list if not os.path.exists(os.path.join(f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/processed_data/", 'config.pkl'))]
print(f"Samples without config.pkl: {len(sample_list)}")

for sample in sample_list:
    raw_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/raw_data/"
    processed_data_dir = f"/home/wu/datb1/AutoExtractSingleCell/skin_real_world_datasets/{sample}/processed_data/"
    print(f"Processing {sample}...")
    if os.path.exists(os.path.join(raw_data_dir, 'config.pkl')):
        continue
    if os.path.exists(os.path.join(raw_data_dir, f"{sample}_raw.h5ad")):
        try:
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
        except:
            failed_samples.append(sample)
            continue
        
print(f"Failed samples: {failed_samples}")