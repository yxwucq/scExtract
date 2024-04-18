from typing import Optional, Tuple
from oaklib import get_adapter
from oaklib.datamodels.search import SearchConfiguration, SearchTermSyntax, SearchProperty

import os
import numpy as np
import logging
import anndata as ad
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

def request_ols(query_cell_type: str,
                ontology: str = "cl",
                query_fields: str = "label",
                exact: bool = False,
                rows: int = 10,
                retries: int = 5,
                backoff_factor: float = 0.1,
                ) -> Optional[Tuple[str, str]]:
    
    """
    Request OLS API for cell type information, return the first obo_id and label
    """
    
    query_cell_type = '+'.join(query_cell_type.split())
    url = f"http://www.ebi.ac.uk/ols4/api/search?q={query_cell_type}&ontology={ontology}&isDefiningOntology=true"
    
    # default url performs better, which is strange
    if query_fields != "label":
        url += f"&queryFields={query_fields}"
    if exact:
        url += f"&exact={exact}"
    if rows != 10:
        url += f"&rows={rows}"
    
    s = requests.Session()
    retries = Retry(total=retries,
                    backoff_factor=backoff_factor,
                    status_forcelist=[ 500, 502, 503, 504 ])
    
    s.mount('http://', HTTPAdapter(max_retries=retries))
    
    resp = s.get(url)
    results = resp.json()['response']['docs']
    
    if len(results) > 0:
        for i in range(len(results)):
            if results[i]['obo_id'].startswith(ontology.upper()):
                return results[i]['obo_id'], results[i]['label']
    else:
        return None, None

def request_ols_ubergraph(query_cell_type: str,
                          ontology: str = "cl") -> Optional[Tuple[str, str]]:
                          
    adapter = get_adapter(f"ubergraph:{ontology}")
    for curie in adapter.basic_search(query_cell_type):
        return curie, adapter.label(curie)
    return None, None

def ontology_pairwise_similarity(obo_id1: Optional[str],
                                 obo_id2: Optional[str],
                                 ontology: str = "cl",
                                 ) -> float:
    
    if obo_id1 is None or obo_id2 is None:
        return None
    
    adapter = get_adapter(f"sqlite:obo:{ontology}")
    results = adapter.pairwise_similarity(obo_id1, obo_id2)
    
    return round(results['jaccard_similarity'], 3)  

def benchmark_annotation(adata_path : str,
                         true_group_key: str,
                         predict_group_key: str = None,
                         ontology: str = "cl",
                         method: str = "ols_api",
                         similarity_key: str = "similarity",
                         output_path: str = None
                        ) -> None:
    
    adata = ad.read_h5ad(adata_path)
        
    if predict_group_key is None:
        if 'leiden' in adata.obs.keys() and 'louvain' in adata.obs.keys():
            raise ValueError("Multiple clustering results found. \
                             Please specify the predict_group_key.")
        elif 'leiden' in adata.obs.keys():
            predict_group_key = 'leiden'
        elif 'louvain' in adata.obs.keys():
            predict_group_key = 'louvain'
        else:
            raise ValueError("No clustering results found. \
                             Please specify the predict_group_key.")
    
    if method == "ols_api":
        request_func = request_ols
    elif method == "ols_ubergraph":
        request_func = request_ols_ubergraph
    else:
        raise ValueError("Unknown method.")
    
    predict_group_key_list = predict_group_key.split(',')
    similarity_key_list = similarity_key.split(',')

    for predict_group_key, similarity_key in zip(predict_group_key_list, similarity_key_list):
        logging.info(f"Processing Predict group key: {predict_group_key}, Add to Similarity key: {similarity_key}")
        adata.obs[predict_group_key + '_cl_label'] = None
        adata.obs[predict_group_key + '_cl_obo_id'] = None
        adata.obs[true_group_key + '_cl_label'] = None
        adata.obs[true_group_key + '_cl_obo_id'] = None
        
        predict_label_dict, predict_obo_dict = {}, {}
        true_label_dict, true_obo_dict = {}, {}
        
        for cell_type in adata.obs[predict_group_key].unique():
            obo_id, label = request_func(cell_type, ontology=ontology)
            logging.info(f"Predicted: {cell_type}, Label: {label}, OBO ID: {obo_id}")
            predict_label_dict[cell_type] = label
            predict_obo_dict[cell_type] = obo_id
            
        for cell_type in adata.obs[true_group_key].unique():
            obo_id, label = request_func(cell_type, ontology=ontology)
            logging.info(f"True: {cell_type}, Label: {label}, OBO ID: {obo_id}")
            true_label_dict[cell_type] = label
            true_obo_dict[cell_type] = obo_id
        
        adata.obs[predict_group_key + '_cl_label'] = adata.obs[predict_group_key].map(predict_label_dict)
        adata.obs[predict_group_key + '_cl_obo_id'] = adata.obs[predict_group_key].map(predict_obo_dict)
        adata.obs[true_group_key + '_cl_label'] = adata.obs[true_group_key].map(true_label_dict)
        adata.obs[true_group_key + '_cl_obo_id'] = adata.obs[true_group_key].map(true_obo_dict)
        
        similarity_dict = {}
        similarity_dict[(None, None)] = None
        
        for predict_cl in predict_obo_dict.values():
            for true_cl in true_obo_dict.values():
                similarity_dict[(predict_cl, true_cl)] = ontology_pairwise_similarity(predict_cl, true_cl, ontology)
                
        adata.obs[similarity_key] = [similarity_dict[(x, y)] for x, y in zip(adata.obs[predict_group_key + '_cl_obo_id'], 
                                                                        adata.obs[true_group_key + '_cl_obo_id'])]
    
    for similarity_key in similarity_key_list:
        logging.info(f"Similarity key: {similarity_key}, Mean: {np.mean(adata.obs[similarity_key])}")
    
    if output_path is not None:
        adata.write(output_path)
    else:
        adata.write(adata_path)