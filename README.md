# scExtract: Automatic annotation and integration of single-cell RNA-seq datasets from the literature

scExtract is a tool for automating the extraction, processing, and annotation of single-cell data from literature. The tool is designed to use LLMs to extract relevant information from scientific articles and process the data.

The input of the program are an anndata object and a PDF/txt file containing the article from which the single-cell data is to be extracted. 

## Usage

- Step1: Clone repo from github
```
git clone https://github.com/yxwucq/scExtract
```

- Step2: Copy `auto_extract/config_sample.py` to `auto_extract/config.py`, then fill your api provider(Claude3 and OpenAI models are both supported) in `auto_extract/config.py`:

```
class Config:
    API_KEY = 'YOUR_API_KEY'
    API_BASE_URL = "YOUR_API_BASE_URL" # for third party openai api
    TYPE = "openai" # claude or openai
    MODEL = "claude-3-sonnet-20240229" # model processing the article
    TOOL_MODEL = "claude-3-opus-20240229" # model for short messages
```

If you want to benchmark auto-annotation using similarity of annotated-text or later integrate your datasets, you need convert annotation to embedding using LLM, by setting `CONVERT_EMBEDDING = True` in `config.py`:

```
class Config:
    ...
    CONVERT_EMBEDDING = True
    EMBEDDING_MODEL = 'text-embedding-3-large'
    API_STYLES = '' # values from ['azure', 'same', 'openai']
    EMBEDDING_API_KEY = ''
    EMBEDDING_ENDPOINT = ''
```

If you using GPTs for text extraction, you can set `API_STYLES = same` to use same setting. But there are no official text-to-embedding model in Claudes, you can set `API_STYLES = azure|openai` and according api key and endpoint to use Openai t2e model while using Claudes for text extraction.

- Step3: Add requirements and install locally
```
pip install -f requirements.txt
pip install -e .
```

Finally, directly excute by typing `scExtract`:

```
scExtract -h   
```

### Annotate

Using `auto_extract` subcommand.

```
scExtract auto_extract \
    -i ADATA_PATH \
    -p PDF_PATH \
    -d OUTPUT_DIR \
    -o OUTPUT_NAME
```

The extraction follows these steps, the processing decisions/parameters are all article-based.
1. Filter: Including Min_genes, Min_cells, Mitochondria_counts_percentage, etc.
2. Preprocess: Including Normalization, Log1p_transform, Highly_variable_genes_selection, etc.
3. Unsupervised clustering: Leiden clustering or Louvain clustering to similar groups as defined in the text
4. Marker gene identification: Find marker genes based on differential expression
5. Annotation: Cell type annotation
6. Reannotating clusters(Optional): Query the low-confidence annotations associated gene expression and reannotate them

For detailed configuration, refer to `config.py`.

### Integration

Using `integrate` subcommand.

```
scExtract integrate \
    -f FILE_LIST \
    -m scExtract Method to use for integration. Support scExtract and cellhint. \
    --prior_method llm Method to use for creating the prior similarity matrix. Support ontology, llm and local.
```

For large dataset computed on HPC without internet access, you can first generate text embedding by individual dataset using `python integration/extract_celltype_embedding.py`. Then integrate using local provided dict object:

```
scExtract integrate \
    -f FILE_LIST \
    -m scExtract Method to use for integration. Support scExtract and cellhint. \
    --prior_method local \
    --embedding_dict_path EMBEDDING_DICT_PATH Path to the cell type embedding dictionary.
```

### Benchmark

For benchmark, using `benchmark` subcommand:
```
scExtract benchmark \
    -i INPUT_ADATA \
    -o OUTPUT_ADATA \
    -r METRICS_FILE \
    --true_group_key TRUE_KEY \
    --predict_group_key KEY1,KEY2,... \
    --similarity_key SIMI_KEY1,SIMI_KEY2,...
```

## Other methods
### add singleR annotation
see `scExtract add_singler_annotation -h`

## Example

![Similarity to author-defined cell type](src/similarity.png)

### sample
Muto, Y., Wilson, P.C., Ledru, N. et al. Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. Nat Commun 12, 2190 (2021). https://doi.org/10.1038/s41467-021-22368-w

`Cell_type` is author-defined cell type, `scExtract` is cell type extracted from scExtract, `no_context_anno` is cell type extracted without context information, `singler` is cell type from singleR. `Tissue`, `Certainty` are from scExtract.

![sample8](src/sample8_benchmark.png)
