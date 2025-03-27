from typing import Dict, Any
import logging

import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad
import warnings

from .parse_params import Params

def is_normalized_and_log1p(expression_matrix: np.ndarray) -> Dict[str, bool]:
    """
    Check if the expression matrix is normalized to a specific total counts per cell and log1p transformed.

    Parameters:
    expression_matrix (np.ndarray): The expression matrix to check.

    Returns:
    dict: A dictionary with keys 'normalized' and 'log1p' indicating whether the matrix is normalized and log1p transformed.
    """
    result = {'normalized': False, 'log1p': False}

    # Check if the matrix is normalized to a specific total counts per cell
    total_counts = expression_matrix.sum(axis=1)
    total_counts_per_cell = expression_matrix.sum(axis=1).mean()
    if np.allclose(total_counts, total_counts_per_cell, atol=100):  # Allowing some tolerance
        result['normalized'] = True

    # Check if the matrix is log1p transformed
    # Log1p transformed data should have most values between 0 and a small positive number
    if np.all(expression_matrix >= 0) and np.all(expression_matrix < 10):
        result['log1p'] = True

    return result

def filter(adata: ad.AnnData, 
           params: Params,
           ) -> ad.AnnData:

    if type(adata.X) is np.ndarray: # possible nan value
        adata.X = csr_matrix(np.nan_to_num(adata.X).copy())

    check_result = is_normalized_and_log1p(adata.X.toarray())
    if check_result['normalized'] or check_result['log1p']: 
        logging.warning("The expression matrix is already normalized and/or log1p transformed. \
                       The filtering steps will be skipped.")
        return adata

    # Filter cells
    if params['filter_cells_min_genes'] is not None:
        sc.pp.filter_cells(adata, min_genes=min(params['filter_cells_min_genes'], 500))
    if params['filter_cells_max_genes'] is not None:
        sc.pp.filter_cells(adata, max_genes=max(params['filter_cells_max_genes'], 5000))
    if params['filter_cells_min_counts'] is not None:
        sc.pp.filter_cells(adata, min_counts=min(params['filter_cells_min_counts'], 2000))
    if params['filter_cells_max_counts'] is not None:
        sc.pp.filter_cells(adata, max_counts=max(params['filter_cells_max_counts'], 40000))
    if params['filter_genes_min_cells'] is not None:
        sc.pp.filter_genes(adata, min_cells=min(params['filter_genes_min_cells'], 50))
    
    # TODO only human datasets
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPL', 'RPS'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    
    if params['filter_mito_percentage_max'] is not None and adata.var['mt'].sum() > 5: # at least 5 mitochondrial genes
        adata = adata[adata.obs['pct_counts_mt'] < params['filter_mito_percentage_max'], :].copy()
    if params['filter_mito_percentage_min'] is not None and adata.var['mt'].sum() > 5:
        adata = adata[adata.obs['pct_counts_mt'] > params['filter_mito_percentage_min'], :].copy()
    if params['filter_ribo_percentage_min'] is not None and adata.var['ribo'].sum() > 5:
        adata = adata[adata.obs['pct_counts_ribo'] > params['filter_ribo_percentage_min'], :].copy()
    if params['filter_ribo_percentage_max'] is not None and adata.var['ribo'].sum() > 5: # at least 5 ribosomal genes
        adata = adata[adata.obs['pct_counts_ribo'] < params['filter_ribo_percentage_max'], :].copy()
        
    return adata

def preprocess(adata: ad.AnnData, 
               params: Params,
               ) -> ad.AnnData:
    
    adata = adata.copy() # avoid inplace modification
    adata.layers['counts'] = adata.X.copy()

    check_result = is_normalized_and_log1p(adata.X.toarray())
    if check_result['log1p']:
        check_result['normalized'] = True # if log1p transformed, it is normalized
        adata.layers['counts'] = csr_matrix(np.expm1(adata.X.toarray())) # revert log1p transformation

    # Normalize total target sum
    if params['normalize_total_target_sum'] is not None and not check_result['normalized']:
        sc.pp.normalize_total(adata, target_sum=params['normalize_total_target_sum'])
    elif check_result['normalized']:
        logging.warning("The expression matrix is already normalized. The normalization step will be skipped.")
    
    # Log1p transform
    if params['log1p_transform']:
        sc.pp.log1p(adata)
    elif check_result['log1p']:
        logging.warning("The expression matrix is already log1p transformed. The log1p transformation step will be skipped.")
        
    # Batch correction
    # if params['batch_correction']:
    #     # key is should be specified TODO
    #     sc.pp.combat(adata, key='Batch')

    # Save log1p transformed data
    adata.raw = adata

    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=params['highly_variable_genes_num'])
    
    # Scale
    if params['scale']:
        sc.pp.scale(adata)
    
    # PCA
    sc.pp.pca(adata, n_comps=max(params['pca_comps'], 100))
    
    # Batch correction scanorama
    if params['batch_correction']:
        sc.external.pp.harmony_integrate(adata, 'Batch')
        
    # Find neighbors
    use_rep = 'X_pca_harmony' if params['batch_correction'] else 'X_pca'
    sc.pp.neighbors(adata, n_neighbors=params['find_neighbors_neighbors_num'], use_rep=use_rep, n_pcs=params['find_neighbors_using_pcs'])
    
    return adata

def clustering(adata: ad.AnnData, 
               params: Params,
               ) -> ad.AnnData:
    
    # Unsupervised clustering
    if params['unsupervised_cluster_method'] == 'leiden':
        clustering_key = 'leiden'
        clustering_method = sc.tl.leiden
    elif params['unsupervised_cluster_method'] == 'louvain':
        clustering_key = 'louvain'
        clustering_method = sc.tl.louvain
    else:
        warnings.warn(f"Invalid clustering method: {params['unsupervised_cluster_method']}, using default: leiden")
        clustering_key = 'leiden'
        clustering_method = sc.tl.leiden
    
    # test resolution to meet the number of groups
    for resolution in [0.1, 0.3, 0.5, 0.7, 0.9]:
        clustering_method(adata, resolution=resolution)
        groups_num = len(adata.obs[clustering_key].unique())
        
        if groups_num >= params['leiden_or_louvain_group_numbers']:
            if abs(groups_num - params['leiden_or_louvain_group_numbers']) > 3:
                warnings.warn(f"Detected {groups_num} groups, which quite differs from the expected {params['leiden_or_louvain_group_numbers']} groups")
            break
        
    # Visualize
    use_rep = 'X_pca_harmony' if params['batch_correction'] else 'X_pca'
    if params['visualize_method'] == 'umap':
        sc.tl.umap(adata)
    elif params['visualize_method'] == 'tsne':
        sc.tl.tsne(adata, use_rep=use_rep)
        
    return adata