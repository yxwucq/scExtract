from typing import Optional, Tuple
from oaklib import get_adapter
from oaklib.datamodels.search import SearchConfiguration, SearchTermSyntax, SearchProperty
from sklearn.metrics import adjusted_rand_score

import os
import numpy as np
import logging
import anndata as ad
import pandas as pd
import pickle
import requests
from requests.adapters import HTTPAdapter, Retry

from ..auto_extract.parse_params import Params 
from ..auto_extract.agent import get_cell_type_embedding_by_llm

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
    
    # remove plural to make the query more accurate
    if query_cell_type.endswith('s'):
        query_cell_type = query_cell_type[:-1]
    
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
        return None, None
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
                         config_path: str = 'config.ini',
                         predict_group_key: Optional[str] = None,
                         method: str = "embedding",
                         output_config_pkl: Optional[str] = None,
                         similarity_key: str = "similarity",
                         ontology: str = "cl",
                         result_metrics_path: str = None,
                         output_path: str = None) -> None:
    
    if 'ols' in method:
        benchmark_annotation_ols(adata_path = adata_path,
                                 true_group_key = true_group_key,
                                 predict_group_key = predict_group_key,
                                 method = method,
                                 similarity_key = similarity_key,
                                 ontology = ontology,
                                 result_metrics_path = result_metrics_path,
                                 output_path = output_path)
        
    elif method == 'embedding':
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

        predict_group_key_list = predict_group_key.split(',')
        similarity_key_list = similarity_key.split(',')

        if output_config_pkl is not None:
            with open(output_config_pkl, 'rb') as f:
                config = pickle.load(f)
        else: # avoid loading config file
            config = Params(config_path)
            config.embedding_dict = {}
        
        ari_list = []
        predict_group_key_list_cp = predict_group_key_list.copy()
        similarity_key_list_cp = similarity_key_list.copy()
        
        for predict_group_key, similarity_key in zip(predict_group_key_list, similarity_key_list):
            if predict_group_key == 'scextract':
                if 'leiden' in adata.obs.keys():
                    predict_group_key = 'leiden'
                elif 'louvain' in adata.obs.keys():
                    predict_group_key = 'louvain'
                else:
                    raise ValueError("No clustering results found.")
            
            logging.info(f"Processing Predict group key: {predict_group_key}, Add to Similarity key: {similarity_key}")
            predict_cell_type_list = adata.obs[predict_group_key].unique()
            true_cell_type_list = adata.obs[true_group_key].unique()
            
            predict_cell_type_list = [x for x in predict_cell_type_list if x != '']
            if len(predict_cell_type_list) == 0:
                logging.warning(f"No cell type found in {predict_group_key}")
                predict_group_key_list_cp.remove(predict_group_key)
                similarity_key_list_cp.remove(similarity_key)
                continue
            
            ari_score = adjusted_rand_score(adata.obs[true_group_key], adata.obs[predict_group_key])
            ari_list.append(ari_score)
            
            similarity_dict = {}
            similarity_dict[(None, None)] = None
            
            # get cell type embedding by llm
            if not all([x in config.embedding_dict.keys() for x in predict_cell_type_list]):
                embedding_list = get_cell_type_embedding_by_llm(predict_cell_type_list, config_path=config_path)
                for i in range(len(predict_cell_type_list)):
                    config.embedding_dict[predict_cell_type_list[i]] = embedding_list[i]
            
            if not all([x in config.embedding_dict.keys() for x in true_cell_type_list]):
                embedding_list = get_cell_type_embedding_by_llm(true_cell_type_list, config_path=config_path)
                for i in range(len(true_cell_type_list)):
                    config.embedding_dict[true_cell_type_list[i]] = embedding_list[i]       

            for x in predict_cell_type_list:
                for y in true_cell_type_list:
                    similarity_dict[(x, y)] = np.dot(config.embedding_dict[x], config.embedding_dict[y]) / \
                                                (np.linalg.norm(config.embedding_dict[x]) * np.linalg.norm(config.embedding_dict[y]))

            adata.obs[similarity_key] = adata.obs.apply(lambda x: similarity_dict[(x[predict_group_key], x[true_group_key])], axis=1)
        
        similarity_key_list = similarity_key_list_cp
        predict_group_key_list = predict_group_key_list_cp
        
        for similarity_key in similarity_key_list:
            logging.info(f"Similarity key: {similarity_key}, Mean: {np.mean(adata.obs[similarity_key])}")
        
        for i, predict_group_key in enumerate(predict_group_key_list):
            logging.info(f"ARI for {predict_group_key}: {ari_list[i]}")
        
        with open(result_metrics_path, 'w') as f:
            for i, predict_group_key in enumerate(predict_group_key_list):
                f.write(f"ARI for {predict_group_key}: {ari_list[i]}\n")
            for similarity_key in similarity_key_list:
                f.write(f"Similarity key: {similarity_key}, Mean: {np.mean(adata.obs[similarity_key])}\n")
                f.write(f"Group-level Similarity key: {similarity_key}, Mean: {adata.obs.groupby(true_group_key)[similarity_key].mean().mean()}\n")
        
        if output_path is not None:
            adata.write(output_path)
        else:
            adata.write(adata_path)    
            
    else:
        raise ValueError("Unknown method.")

def benchmark_annotation_ols(adata_path : str,
                         true_group_key: str,
                         predict_group_key: Optional[str] = None,
                         ontology: str = "cl",
                         method: str = "ols_api",
                         similarity_key: str = "similarity",
                         result_metrics_path: str = None,
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

    ari_list = []
    predict_group_key_list_cp = predict_group_key_list.copy()
    similarity_key_list_cp = similarity_key_list.copy()
        
    for predict_group_key, similarity_key in zip(predict_group_key_list, similarity_key_list):
        if predict_group_key == 'scextract':
            if 'leiden' in adata.obs.keys():
                predict_group_key = 'leiden'
            elif 'louvain' in adata.obs.keys():
                predict_group_key = 'louvain'
            else:
                raise ValueError("No clustering results found.")
        
        ari_score = adjusted_rand_score(adata.obs[true_group_key], adata.obs[predict_group_key])
        ari_list.append(ari_score)
        
        logging.info(f"Processing Predict group key: {predict_group_key}, Add to Similarity key: {similarity_key}")
        predict_cell_type_list = adata.obs[predict_group_key].unique()
        
        predict_cell_type_list = [x for x in predict_cell_type_list if x != '']
        if len(predict_cell_type_list) == 0:
            logging.warning(f"No cell type found in {predict_group_key}")
            predict_group_key_list_cp.remove(predict_group_key)
            similarity_key_list_cp.remove(similarity_key)
            continue
        
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
    
    similarity_key_list = similarity_key_list_cp
    predict_group_key_list = predict_group_key_list_cp
    
    for similarity_key in similarity_key_list:
        logging.info(f"Similarity key: {similarity_key}, Mean: {np.mean(adata.obs[similarity_key])}")
    
    for i, predict_group_key in enumerate(predict_group_key_list):
        logging.info(f"ARI for {predict_group_key}: {ari_list[i]}")

    with open(result_metrics_path, 'w') as f:
        for i, predict_group_key in enumerate(predict_group_key_list):
            f.write(f"ARI for {predict_group_key}: {ari_list[i]}\n")
        for similarity_key in similarity_key_list:
            f.write(f"Similarity key: {similarity_key}, Mean: {np.mean(adata.obs[similarity_key])}\n")
            f.write(f"Group-level Similarity key: {similarity_key}, Mean: {adata.obs.groupby(true_group_key)[similarity_key].mean().mean()}\n")
 
    if output_path is not None:
        adata.write(output_path)
    else:
        adata.write(adata_path)