from typing import List, Optional
from tqdm import tqdm
import anndata as ad
import os
import cellhint_prior
from scipy.sparse import csr_matrix
from openai import AzureOpenAI, OpenAI
from scipy.sparse.csgraph import minimum_spanning_tree
import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import argparse
import logging
import sys

from ..auto_extract.agent import get_cell_type_embedding_by_llm

def extract_celltype_embedding(file_list: List[str],
                               output_embedding_pkl: str,
                               output_individual_config_pkls: Optional[str] = None,
                               cell_type_column: Optional[str] = None,
                               config_path: str = 'config.ini',
                               ):
    """
    Extract cell type embeddings from the processed datasets.
    """
    
    embeddings_dict = {}
    
    logging.info(f"Extracting cell type embeddings from {len(file_list)} processed datasets.")
    for file_path in tqdm(file_list):
        adata = ad.read_h5ad(file_path)
        if cell_type_column is not None:
            cell_types = adata.obs[cell_type_column].unique().tolist()
        else:
            if 'cell_type' in adata.obs.columns:
                cell_types = adata.obs['cell_type'].astype(str).str.replace('_', ' ').unique().tolist()
            elif 'leiden' in adata.obs.columns:
                cell_types = adata.obs['leiden'].astype(str).str.replace('_', ' ').unique().tolist()
            elif 'louvain' in adata.obs.columns:
                cell_types = adata.obs['louvain'].astype(str).str.replace('_', ' ').unique().tolist()
            else:
                raise ValueError("Unknown cell type column.")
        logging.info(f"Writing {len(cell_types)} cell types to the embedding dictionary.")
        
        cell_types = [x.replace('/', '|') for x in cell_types]
        
        if output_individual_config_pkls is not None:
            file_config_path_list = output_individual_config_pkls.split()
            for config_path in file_config_path_list:
                if os.path.exists(output_individual_config_pkls):
                    with open(config_path, 'rb') as f:
                        config = pickle.load(f)
                    if not all([x in config.embedding_dict.keys() for x in cell_types]):
                        emb = get_cell_type_embedding_by_llm(cell_types, config_path=config_path)
                        for i in range(len(cell_types)):
                            embeddings_dict[cell_types[i]] = emb[i]
        else:
            emb = np.array(get_cell_type_embedding_by_llm(cell_types, config_path=config_path))
        for key, value in zip(cell_types, emb):
            embeddings_dict[key] = value
    
        del adata
        
    logging.info(f"Toal {len(embeddings_dict)} cell type embeddings are extracted.")
    
    with open(output_embedding_pkl, 'wb') as f:
        pickle.dump(embeddings_dict, f)