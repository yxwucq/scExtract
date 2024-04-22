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
options:
usage: main.py [-h] {auto_extract,benchmark,add_singler_annotation} ...

Automatically extract, process, and annotate single-cell data from literature.

positional arguments:
  {auto_extract,benchmark,add_singler_annotation}
                        sub-command help
    auto_extract        Automatically extract, process, and annotate single-cell data from literature.
    benchmark           Benchmark annotation results using true labels.
    add_singler_annotation
                        Annotate single-cell data using py&c++ implementation of singler.

options:
  -h, --help            show this help message and exit
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

For a whole process, including extract annotation, add other method, compare with ground truth for benchmarking, you can directly run
`python pipelines.py sample{i}` in the following folder structure:

```
.
├── processed_data
│   └── sample{i}_true.h5ad # Contains `cell_type` col in obs for benchmarking
└── raw_data
    ├── sample{i}.pdf
    └── sample{i}_raw.h5ad # Contains `Batch` col in obs for possible batch correction
```

### step-by-step

run ` python main.py benchmark` function
```
usage: main.py benchmark [-h] --adata_path ADATA_PATH [--output_path OUTPUT_PATH] --true_group_key TRUE_GROUP_KEY [--predict_group_key PREDICT_GROUP_KEY]
                         [--ontology ONTOLOGY] [--method METHOD] [--similarity_key SIMILARITY_KEY]

options:
  -h, --help            show this help message and exit
  --adata_path ADATA_PATH, -i ADATA_PATH
                        Path to the processed data in AnnData format.
  --output_path OUTPUT_PATH, -o OUTPUT_PATH
                        Path to save the output file. If not specified, the input file will be overwritten.
  --true_group_key TRUE_GROUP_KEY, -t TRUE_GROUP_KEY
                        Key of the true group in adata.obs.
  --predict_group_key PREDICT_GROUP_KEY, -p PREDICT_GROUP_KEY
                        Key of the predicted group in adata.obs. Support multiple keys separated by comma.
  --ontology ONTOLOGY, -l ONTOLOGY
                        Ontology to use for annotation.
  --method METHOD, -m METHOD
                        Method to use for annotation.
  --similarity_key SIMILARITY_KEY, -s SIMILARITY_KEY
                        Key to save the similarity results. Support multiple keys separated by comma. Order should be the same as predict_group_key.
```

## Example
### sample1
Wang, S., Drummond, M.L., Guerrero-Juarez, C.F. et al. Single cell transcriptomics of human epidermis identifies basal stem cell transition states. Nat Commun 11, 4239 (2020). https://doi.org/10.1038/s41467-020-18075-7

![sample1](src/sample1_benchmark.png)

## Other methods
### singleR
see `python main.py add_singler_annotation -h`
```
usage: main.py add_singler_annotation [-h] --adata_path ADATA_PATH [--output_path OUTPUT_PATH] [--ref_data REF_DATA] [--ref_features REF_FEATURES]
                                      [--ref_labels REF_LABELS] [--cache_dir CACHE_DIR]

options:
  -h, --help            show this help message and exit
  --adata_path ADATA_PATH, -i ADATA_PATH
                        Path to the processed data in AnnData format.
  --output_path OUTPUT_PATH, -o OUTPUT_PATH
                        Path to save the output file. If not specified, the input file will be overwritten.
  --ref_data REF_DATA, -d REF_DATA
                        Reference data to use for annotation.
  --ref_features REF_FEATURES, -f REF_FEATURES
                        Reference features to use for annotation.
  --ref_labels REF_LABELS, -l REF_LABELS
                        Reference labels to use for annotation.
  --cache_dir CACHE_DIR, -c CACHE_DIR
                        Directory to save the cache files.
```

### LLM without article context
In `python main.py auto_extract`, add `--benchmark_no_context_key`
```
--benchmark_no_context_key BENCHMARK_NO_CONTEXT_KEY, -b BENCHMARK_NO_CONTEXT_KEY
                        If specified, Directly get annotation from marker genes without article context for benchmarking, the result will be saved in
                        adata.obs[benchmark_no_context_key].
```