import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from benchmark.benchmark import benchmark_annotation

import scanpy as sc

sample = 'sample1'

if sample == 'sample1':
    adata = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed.h5ad')
    adata_true = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/adata_all_true.h5ad')
    adata = adata[adata.obs.index.isin(adata_true.obs.index)].copy()
    adata.obs['cell_type'] = adata_true.obs['leiden']

    adata_benchmark = benchmark_annotation(
        adata=adata,
        true_group_key='cell_type',
        predict_group_key='louvain',
        ontology='cl',
        method='ols_api'
    )

    print(adata_benchmark.obs.similarity.value_counts())
    adata_benchmark.write_h5ad('/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed_benchmark.h5ad')
    

