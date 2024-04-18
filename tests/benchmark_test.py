import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from benchmark.benchmark import benchmark_annotation

import scanpy as sc

sample = 'sample1'
print(f'Running {sample} test')

if sample == 'sample1':
    adata = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed.h5ad')
    adata_true = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/adata_all_true.h5ad')
    adata = adata[adata.obs.index.isin(adata_true.obs.index)].copy()
    adata.obs['cell_type'] = adata_true.obs['leiden']
    adata.write('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed.h5ad')

    benchmark_annotation(
        adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed.h5ad',
        true_group_key='cell_type',
        predict_group_key='leiden,no_context_annotation,singler_annotation',
        ontology='cl',
        method='ols_api',
        similarity_key='similarity_scextract,similarity_no_context_annotation,similarity_singler'
    )

elif sample == 'sample2':
    adata = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/processed_data/processed.h5ad')
    adata_true = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/processed_data/adata_all_true_extended.h5ad')
    adata = adata[adata.obs.index.isin(adata_true.obs.index)].copy()
    adata.obs['cellxgene_cell_type'] = adata_true.obs['cell_type']
    adata.obs['cell_type'] = adata_true.obs['author_cell_type_expand']
    
    benchmark_annotation(
        adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/processed_data/processed.h5ad',
        true_group_key='cell_type',
        predict_group_key='louvain',
        ontology='cl',
        method='ols_api',
        similarity_key='similarity_scextract'
    )