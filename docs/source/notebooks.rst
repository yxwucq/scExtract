=======================================
Notebooks for Preprocessing GSE Data
=======================================

We provide a series of notebooks to help you preprocess your data. These notebooks are designed to be easy to use on public GSE datasets.
The notebooks are designed to locate in :code:`raw_data` folder. You can choose to run the notebooks according to the types of data you have:

Batched data with 10X format: `batch_10X Example Notebook <https://github.com/yxwucq/scExtract/blob/master/docs/source/_notebooks/batch_10X.ipynb>`_

Batched data with csv.gz format: `batch_csv Example Notebook <https://github.com/yxwucq/scExtract/blob/master/docs/source/_notebooks/batch_csv.ipynb>`_

if txt.gz format, change delimiter to :code:`\\t`

Single data with 10X format: `single_10X Example Notebook <https://github.com/yxwucq/scExtract/blob/master/docs/source/_notebooks/single_10X.ipynb>`_

Single data with txt.gz format: `single_txt Example Notebook <https://github.com/yxwucq/scExtract/blob/master/docs/source/_notebooks/single_txt.ipynb>`_

if csv.gz format, change delimiter to :code:`,`

Note some older 10X data using genes.tsv instead of features.tsv, you should add :code:`Gene Expression` to the column names to make 
it compatible with the latest version of :code:`sc.read_10x_mtx()`.

.. code-block:: python

    for sample in sample_list:
    os.mkdir(os.path.join(dir_path, sample))
    for sample_file in file_list:
        if sample in sample_file:
            os.system(f'mv {os.path.join(dir_path, sample_file)} {os.path.join(dir_path, sample)}')
            if 'barcodes' in sample_file:
                os.system(f'mv {os.path.join(dir_path, sample, sample_file)} {os.path.join(dir_path, sample, "barcodes.tsv.gz")}')
            elif 'genes' in sample_file:
                os.system(f'gunzip {os.path.join(dir_path, sample, sample_file)}')
                os.system(f"sed -i 's/$/\tGene Expression/' {os.path.join(dir_path, sample, sample_file.split('.gz')[0])}")
                os.system(f"gzip {os.path.join(dir_path, sample, sample_file.split('.gz')[0])}")
                os.system(f'mv {os.path.join(dir_path, sample, sample_file)} {os.path.join(dir_path, sample, "features.tsv.gz")}')
            elif 'matrix' in sample_file:
                os.system(f'mv {os.path.join(dir_path, sample, sample_file)} {os.path.join(dir_path, sample, "matrix.mtx.gz")}')
