import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from auto_extract.auto_extract import auto_extract

sample = 'sample1'
print(f'Running {sample} test')

if sample == 'sample1':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/raw_data/raw_data_batch.h5ad'
    pdf_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/raw_data/extracted_text.txt'
    output_dir='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data'
    output_name='processed.h5ad'
    benchmark_no_context_key='no_context_annotation'

elif sample == 'sample2':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/raw_data/BMO_raw.h5ad'
    pdf_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/raw_data/s41592-024-02172-2.pdf'
    output_dir='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/processed_data'
    output_name='processed.h5ad'
    benchmark_no_context_key='no_context_annotation'

auto_extract(
    adata_path=adata_path,
    pdf_path=pdf_path,
    output_dir=output_dir,
    output_name=output_name,
    benchmark_no_context_key=benchmark_no_context_key
)