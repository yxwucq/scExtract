============
Installation
============

System Requirements
=====================

- Python 3.10 or higher
- 64-bit operating system (Linux, macOS, or Windows)
- Minimum 8GB RAM (16GB or more recommended)

Basic Installation
=====================

Install the latest version from Source:

.. code-block:: bash

    # Create a new conda environment
    conda create -n scExtract python=3.10
    conda activate scExtract

    # Install scExtract
    git clone https://github.com/yxwucq/scExtract
    cd scExtract
    pip install -e .

Dependencies
=============

Core Dependencies
-----------------

scExtract depends on the following Python packages:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Package
     - Version
     - Description
   * - python
     - >=3.10 <3.13
     - Programming language
   * - anndata
     - >=0.10.7
     - Data structure for single-cell data
   * - scanpy
     - >=1.10.1
     - Single-cell analysis tools
   * - harmonypy
     - >=0.0.9
     - Batch effect correction
   * - leidenalg
     - >=0.10.2
     - Community detection
   * - louvain
     - >=0.8.2
     - Clustering algorithm
   * - mygene
     - >=3.2.2
     - Gene annotation
   * - oaklib
     - >=0.6.3
     - Ontology analysis
   * - openai
     - >=1.17.0
     - OpenAI API integration
   * - pypdf
     - >=4.2.0
     - PDF processing
   * - singler
     - >=0.1.2
     - Cell type annotation
   * - singlecellexperiment
     - >=0.4.4
     - Single cell data structure
   * - termcolor
     - >=2.4.0
     - Terminal color output
   * - tqdm
     - >=4.66.2
     - Progress bars

Optional Dependencies
---------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Package
     - Version
     - Purpose
   * - anthropic
     - >=0.25.1
     - Anthropic API integration

