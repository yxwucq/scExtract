import sys
sys.path.append('/home/wu/datb1/AutoExtractSingleCell/scExtract')
from scextract.auto_extract.auto_extract import auto_extract

config_path = '/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/config.ini'
sample = 'sample1'
print(f'Running {sample} test')

if sample == 'sample1':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/sample1_revision/raw_data/sample1_raw.h5ad'
    pdf_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/sample1_revision/raw_data/sample1.pdf'
    output_dir='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample1/processed_data'
    output_name='processed.h5ad'
    benchmark_no_context_key='no_context_annotation'

elif sample == 'sample2':
    adata_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/raw_data/BMO_raw.h5ad'
    pdf_path='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/raw_data/s41592-024-02172-2.pdf'
    output_dir='/home/wu/datb1/AutoExtractSingleCell/scExtract/tests/test_data/sample2/processed_data'
    output_name='processed.h5ad'
    benchmark_no_context_key='no_context_annotation'

if __name__ == '__main__':
    auto_extract(
        adata_path=adata_path,
        pdf_path=pdf_path,
        output_dir=output_dir,
        output_name=output_name,
        benchmark_no_context_key=benchmark_no_context_key,
        config_path=config_path
    )