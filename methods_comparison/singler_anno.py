# can also try r2py version of the code
# https://gist.github.com/gokceneraslan/f953b697ea09b9a8ff890a8707775498

import anndata as ad
import h5py
import numpy as np
from scipy.sparse import csr_matrix
from pathlib import Path
from typing import List
import os
import singler
import singlecellexperiment as sce

def write_10X_h5(adata: ad.AnnData, 
                 file: str,
                 layer: str = 'counts',
                 ) -> None:
    """Writes adata to a 10X-formatted h5 file.
    
    Note that this function is not fully tested and may not work for all cases.
    It will only write the following keys to the h5 file compared to 10X:
    - X
    - obs_names
    - var_names

    Args:
        adata (AnnData object): processed AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Returns:
        None
    """
    if layer == 'X':
        if type(adata.X) != csr_matrix:
            raise ValueError("X layer is not a csr_matrix.")
    elif layer not in adata.layers.keys():
        raise ValueError(f"No {layer} layer found in adata.")
    else:
        adata.X = adata.layers[layer]
    
    if '.h5' not in file: file = f'{file}.h5'
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))

def singler_annotation(h5file_path: str,
                       ref_data: str = "HumanPrimaryCellAtlasData",
                       ref_features: str = "symbol",
                       ref_labels: str = "main",
                       cache_dir: str = "_cache",
                       ) -> List[str]:
    
    if not h5file_path.endswith('.h5'):
        raise ValueError("File must be in h5 format.")
    
    data = sce.read_tenx_h5(h5file_path)
    mat = data.assay("counts")
    features = [str(x) for x in data.row_data["name"]]

    results = singler.annotate_single(
        mat,
        features,
        ref_data = ref_data,
        ref_features = ref_features,
        ref_labels = ref_labels,
        cache_dir = cache_dir
    )

    return results.column("best")

def add_singler_annotation(adata_path: str,
                           h5file_path: str = '_tmpfile.h5',
                           ref_data: str = "HumanPrimaryCellAtlas",
                           ref_features: str = "symbol",
                           ref_labels: str = "main",
                           cache_dir: str = "_cache",
                           key_added: str = "singler_annotation",
                           output_path: str = None,
                           ) -> None:
    """"
    Add singler annotation to AnnData object.
    """
    
    adata = ad.read_h5ad(adata_path)
    
    adata_tmp = adata.copy()
    write_10X_h5(adata_tmp, h5file_path)
    
    del adata_tmp
    
    singler_anno_results = singler_annotation(h5file_path, 
                                              ref_data,
                                              ref_features, 
                                              ref_labels, 
                                              cache_dir)
    
    os.remove(h5file_path)
    adata.obs[key_added] = singler_anno_results
    
    if output_path is not None:
        adata.write(output_path)
    else:
        adata.write(adata_path)