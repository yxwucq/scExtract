import scanpy as sc
import pandas as pd
import anndata as ad
import warnings
from typing import List, Dict

from .parse_params import Params

def get_marker_genes(adata: ad.AnnData,
                     params: Params,
                     fast_mode: bool = False,
                     ) -> tuple[ad.AnnData, dict]:
    """
    Find marker genes for each cluster.
    """
    
    cluster_key = params.get_params['unsupervised_cluster_method']
    
    warnings.filterwarnings('ignore', category=RuntimeWarning)    
    warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
    # use wilcoxon rank sum test to find differentially expressed genes
    if not fast_mode:
        sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon', n_genes=20, tie_correct=True, use_raw=False)
    else:
        sc.tl.rank_genes_groups(adata, cluster_key, method='t-test', n_genes=20, use_raw=False)
    
    marker_genes = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        marker_genes[int(cluster)] = list(adata.uns['rank_genes_groups']['names'][cluster][:10].copy())

    return adata, marker_genes
    
def annotate(adata: ad.AnnData,
             annotation_dict: Dict[int, List[str]],
             params: Params,
             final: bool = False,
             verbose_annot: bool = False,
             ) -> ad.AnnData:
    """
    Annotate clusters based on marker genes.
    """
    
    cluster_key = params.get_params['unsupervised_cluster_method']
    
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(int)
    
    rename_dict = {key: value[0] for key, value in annotation_dict.items()}
    certainty_dict = {key: value[1] for key, value in annotation_dict.items()}
    source_dict = {key: value[2] for key, value in annotation_dict.items()}
    if verbose_annot:
        tissue_dict = {key: value[3] for key, value in annotation_dict.items()}
        disease_dict = {key: value[4] for key, value in annotation_dict.items()}
        developmental_stage_dict = {key: value[5] for key, value in annotation_dict.items()}
    
    adata.obs['Certainty'] = adata.obs[cluster_key].map(certainty_dict)
    adata.obs['Source'] = adata.obs[cluster_key].map(source_dict)
    if verbose_annot:
        adata.obs['Tissue'] = adata.obs[cluster_key].map(tissue_dict)
        adata.obs['Disease'] = adata.obs[cluster_key].map(disease_dict)
        adata.obs['Developmental_stage'] = adata.obs[cluster_key].map(developmental_stage_dict)

    if final:
        adata.obs[cluster_key] = adata.obs[cluster_key].map(rename_dict).str.replace('/', '|').astype('category')
    else:
        adata.obs[f"{cluster_key}_rough"] = adata.obs[cluster_key].map(rename_dict).str.replace('/', '|').astype('category')
    
    return adata

def query_datasets(adata: ad.AnnData,
                   params: Params,
                   percision: int = 3, # Precision of the mean expression values
                     ) -> Dict[str, List[float]]:
    
    query_df = adata.obs.copy()
    # query mean expression of clusters
    mean_expression_dict = {}
    for gene in params['genes_to_query']:
        if gene not in adata.var_names:
            continue
        query_df[gene] = adata.raw[:, gene].X.toarray().flatten()
        mean_expression_dict[gene] = query_df.groupby(params.get_params['unsupervised_cluster_method'])[gene].mean().tolist()
        mean_expression_dict[gene] = [round(x, percision) for x in mean_expression_dict[gene]]
    
    return mean_expression_dict

def simple_annotate(adata: ad.AnnData,
                    annotation_dict: Dict[int, str],
                    params: Params,
                    key_added: str,
                    cluster_key: str = None,
                    ) -> ad.AnnData:
    
    if cluster_key is None:
        cluster_key = params.get_params['unsupervised_cluster_method']
    
    if key_added in adata.obs.keys():
        raise ValueError(f"{key_added} already exists in adata.obs.")
    
    try:
        adata.obs[cluster_key] = adata.obs[cluster_key].astype(int)
    except:
        pass
    
    adata.obs[key_added] = adata.obs[cluster_key].map(annotation_dict)
    adata.obs[key_added] = adata.obs[key_added].astype('category')
    
    return adata