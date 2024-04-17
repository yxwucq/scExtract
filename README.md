# scExtract: Automatic annotation of single-cell RNA-seq data from the literature

scExtract is a tool for automating the extraction, processing, and annotation of single-cell data from literature. The tool is designed to use LLMs to extract relevant information from scientific articles and process the data.

## Usage
The input of the program are an anndata object and a PDF/txt file containing the article from which the single-cell data is to be extracted. 

First copy `auto_extract/config_sample.py` to `auto_extract/config.py`, then fill your api provider(Claude3 and OpenAI models are both supported) in `auto_extract/config.py`:
```
class Config:
    API_KEY = 'YOUR_API_KEY'
    API_BASE_URL = "YOUR_API_BASE_URL" # for third party openai api
    TYPE = "openai" # claude or openai
    MODEL = "claude-3-sonnet-20240229" # model processing the article
    TOOL_MODEL = "claude-3-opus-20240229" # model for short messages
```
Then directly excute through `python main.py`:
```
usage: main.py [-h] --adata_path ADATA_PATH --pdf_path PDF_PATH [--output_dir OUTPUT_DIR] [--output_name OUTPUT_NAME]

Automatically extract, process, and annotate single-cell data from literature.

options:
  -h, --help            show this help message and exit
  --adata_path ADATA_PATH, -i ADATA_PATH
                        Path to the raw data in AnnData format.
  --pdf_path PDF_PATH, -p PDF_PATH
                        Path to the PDF file containing the article. should in pdf or txt format.
  --output_dir OUTPUT_DIR, -d OUTPUT_DIR
                        Directory to save the processed data.
  --output_name OUTPUT_NAME, -o OUTPUT_NAME
                        Name of the output file.
```
The extraction follows these steps, the processing decisions/parameters are all article-based.
1. Filter: Including Min_genes, Min_cells, Mitochondria_counts_percentage, etc.
2. Preprocess: Including Normalization, Log1p_transform, Highly_variable_genes_selection, etc.
3. Unsupervised clustering: Leiden clustering or Louvain clustering to similar groups as defined in the text
4. Marker gene identification: Find marker genes based on differential expression
5. Annotation: Cell type annotation
6. Reannotating clusters(Optional): Query the low-confidence annotations associated gene expression and reannotate them

For detailed configuration, refer to `config.py`.

## Benchmark

run `benchmark_annotation` function in a style like `tests/benchmark_test.py`
```
adata_benchmark = benchmark_annotation(
    adata=adata, # annotated dataset
    true_group_key='cell_type', # manually labeled true value
    predict_group_key='louvain', # cluster method in auto_annotation step
    ontology='cl', # use cell ontology to benchmark
    method='ols_api' # directly use ols_api from ebi
)
```
The returned `adata_benchmark` contains `similarity` in the obs field.

## Example
### sample1
Wang, S., Drummond, M.L., Guerrero-Juarez, C.F. et al. Single cell transcriptomics of human epidermis identifies basal stem cell transition states. Nat Commun 11, 4239 (2020). https://doi.org/10.1038/s41467-020-18075-7

![sample1](src/sample1_benchmark.png)
