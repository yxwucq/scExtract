import scanpy as sc
import anndata as ad
import warnings
from typing import List, Dict

from .config import Config
from .parse_params import Params

def get_marker_genes(adata: ad.AnnData,
                     params: Params,
                     ) -> tuple[ad.AnnData, dict]:
    """
    Find marker genes for each cluster.
    """
    
    cluster_key = params.get_params['unsupervised_cluster_method']
    
    warnings.filterwarnings('ignore', category=RuntimeWarning)    
    # use wilcoxon rank sum test to find differentially expressed genes
    sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon', n_genes=20, tie_correct=True)
    
    marker_genes = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        marker_genes[cluster] = adata.uns['rank_genes_groups']['names'][cluster][:10]

    return adata, marker_genes
    
def annotate(adata: ad.AnnData,
             annotation_dict: Dict[int, List[str]],
             params: Params,
             ) -> ad.AnnData:
    """
    Annotate clusters based on marker genes.
    """
    
    cluster_key = params.get_params['unsupervised_cluster_method']
    
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(int)
    
    rename_dict = {key: value[0] for key, value in annotation_dict.items()}
    certainty_dict = {key: value[1] for key, value in annotation_dict.items()}
    source_dict = {key: value[2] for key, value in annotation_dict.items()}
    tissue_dict = {key: value[3] for key, value in annotation_dict.items()}
    
    adata.obs['Certainty'] = adata.obs[cluster_key].map(certainty_dict)
    adata.obs['Source'] = adata.obs[cluster_key].map(source_dict)
    adata.obs['Tissue'] = adata.obs[cluster_key].map(tissue_dict)

    adata.obs[cluster_key] = adata.obs[cluster_key].map(rename_dict).astype('category')
    
    return adata