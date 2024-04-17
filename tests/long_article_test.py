import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from auto_extract.auto_extract import auto_extract

sample = 'sample1'

if sample == 'sample1':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/raw_data/raw_data_batch.h5ad',
    pdf_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/raw_data/raw_data_batch.h5ad',
    output_dir='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data',
    output_name='processed.h5ad'

auto_extract(
    adata_path=adata_path,
    pdf_path=pdf_path,
    output_dir=output_dir,
    output_name=output_name
)