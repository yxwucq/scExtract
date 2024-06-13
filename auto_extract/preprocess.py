import scanpy as sc
import numpy as np
from scipy.sparse import csc_matrix
import anndata as ad
import warnings

from .config import Config
from .parse_params import Params

def filter(adata: ad.AnnData, 
           params: Params,
           ) -> ad.AnnData:

    if type(adata.X) is np.ndarray: # possible nan value
        adata.X = csc_matrix(np.nan_to_num(adata.X).copy())

    # Filter cells
    if params['filter_cells_min_genes'] is not None:
        sc.pp.filter_cells(adata, min_genes=params['filter_cells_min_genes'])
    if params['filter_cells_max_genes'] is not None:
        sc.pp.filter_cells(adata, max_genes=params['filter_cells_max_genes'])
    if params['filter_cells_min_counts'] is not None:
        sc.pp.filter_cells(adata, min_counts=params['filter_cells_min_counts'])
    if params['filter_cells_max_counts'] is not None:
        sc.pp.filter_cells(adata, max_counts=params['filter_cells_max_counts'])
    if params['filter_genes_low'] is not None:
        sc.pp.filter_genes(adata, min_cells=params['filter_genes_low'])
    
    # TODO only human datasets
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPL', 'RPS'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    
    if params['filter_mito_percentage_low'] is not None and adata.var['ribo'].sum() > 5: # at least 5 ribosomal genes
        adata = adata[adata.obs['pct_counts_mt'] < params['filter_mito_percentage_low'], :].copy()
    if params['filter_ribo_percentage_low'] is not None and adata.var['mt'].sum() > 5: # at least 5 mitochondrial genes
        adata = adata[adata.obs['pct_counts_ribo'] < params['filter_ribo_percentage_low'], :].copy()
        
    return adata

def preprocess(adata: ad.AnnData, 
               params: Params,
               ) -> ad.AnnData:
    
    adata = adata.copy() # avoid inplace modification
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize total target sum
    if params['normalize_total_target_sum'] is not None:
        sc.pp.normalize_total(adata, target_sum=params['normalize_total_target_sum'])
    
    # Log1p transform
    if params['log1p_transform']:
        sc.pp.log1p(adata)
        
    # Batch correction
    # if params['batch_correction']:
    #     # key is should be specified TODO
    #     sc.pp.combat(adata, key='Batch')
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=params['highly_variable_genes_num'])
    
    # Scale
    if params['scale']:
        sc.pp.scale(adata)
    
    # PCA
    sc.pp.pca(adata, n_comps=max(params['pca_comps'], 60))
    
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
    
    # test resolution to meet the number of groups
    for resolution in [0.1, 0.3, 0.5, 0.7, 0.9]:
        clustering_method(adata, resolution=resolution)
        groups_num = len(adata.obs[clustering_key].unique())
        
        if groups_num >= params['leiden_or_louvain_group_numbers']:
            if abs(groups_num - params['leiden_or_louvain_group_numbers']) > 3:
                warnings.warn(f"Detected {groups_num} groups, which quite differs from the \
                              expected {params['leiden_or_louvain_group_numbers']} groups")
            break
        
    # Visualize
    use_rep = 'X_pca_harmony' if params['batch_correction'] else 'X_pca'
    if params['visualize_method'] == 'umap':
        sc.tl.umap(adata)
    elif params['visualize_method'] == 'tsne':
        sc.tl.tsne(adata, use_rep=use_rep)
        
    return adata