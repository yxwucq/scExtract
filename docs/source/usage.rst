===============
Usage
===============

LLM API Configuration
-----------------------

* First, initialize scExtract configuration by running:

.. code-block:: bash

    scExtract init

This will generate a config file (default :code:`config.ini`) with the following structure:

.. code-block:: bash

    [API]
    API_KEY = Fill in your API key here
    API_BASE_URL = Fill in your API base URL here, e.g. https://api.openai.com/v1
    # Supported API styles: openai, claude
    TYPE = openai
    MODEL = Fill in your API model here, e.g. gpt-4o-mini
    TOOL_MODEL = Fill in your API model here, e.g. gpt-4o-mini
    ...

For annotation-to-embedding conversion:

.. code-block:: bash
    
    # Whether to convert the embedding for later use
    CONVERT_EMBEDDING = false
    EMBEDDING_MODEL = text-embedding-3-large
    API_STYLES = azure
    # Note: If you are using the openai API and set the API_STYLES to same,
    # there is no need to specify the EMBEDDING_API_KEY and EMBEDDING_ENDPOINT
    EMBEDDING_API_KEY = Fill in your API key here
    EMBEDDING_ENDPOINT = Fill in your API endpoint here, e.g. https://api.openai.com/v1

- Set `CONVERT_EMBEDDING = true` if you want to:
    - Benchmark auto-annotation using similarity of annotated text
    - Later integrate your datasets in a prior-aware manner

- If using GPTs for text extraction:
    - Set `API_STYLES = same` to use same settings

- If using Claude:
    - Set `API_STYLES = azure|openai` since Claude lacks official text-to-embedding models
    - Provide additional API key and endpoint to use OpenAI's text-to-embedding model

Raw Data Preparation
-----------------------

* :code:`scExtract` accepts raw data in .h5ad format. Raw counts / Normalized counts / Log-transformed counts are all supported and will be automatically detected.
* :code:`Batch` column must be specified in :code:`adata.obs` for possible batch effect correction.
* For most single-cell data forms deposited in public repositories, we provide various notebooks for preprocessing, see :doc:`notebooks`.

Data Processing
-----------------------

* To process raw data, run:

.. code-block:: bash

    scExtract auto_extract -i adata.h5ad -p paper.pdf -o processed.h5ad

This will automatically extract processing parameters from the paper and process the data accordingly.

Benchmark Annotation
-----------------------

* To benchmark cell type annotation, run (suppose the true cell type is stored in :code:`cell_type` column):

.. code-block:: bash

    scExtract benchmark -i processed.h5ad -o benchmark.h5ad -r benchmark_results.csv \
        --true_group_key cell_type --predict_group_key scExtract --similarity_key scExtract_similarity