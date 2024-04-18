import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from methods_comparison.singler_anno import add_singler_annotation

sample = 'sample1'
if sample == 'sample1':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data/processed.h5ad'

add_singler_annotation(adata_path=adata_path)