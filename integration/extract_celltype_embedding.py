from typing import List, Optional
from tqdm import tqdm
import anndata as ad
import os
import cellhint
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

from auto_extract.agent import get_cell_type_embedding_by_llm

def extract_celltype_embedding(file_list: str,
                               output_pkl: str,
                               config_path: str = 'config.pkl',
                               ):
    """
    Extract cell type embeddings from the processed datasets.
    """
    
    embeddings_dict = {}
    
    logging.info(f"Extracting cell type embeddings from {len(file_list)} processed datasets.")
    for file_path in tqdm(file_list):
        adata = ad.read_h5ad(file_path)
        if 'leiden' in adata.obs.columns:
            cell_types = adata.obs['leiden'].astype(str).str.replace('_', ' ').unique().tolist()
        elif 'louvain' in adata.obs.columns:
            cell_types = adata.obs['louvain'].astype(str).str.replace('_', ' ').unique().tolist()
        else:
            raise ValueError("Unknown cell type column.")
        logging.info(f"Writing {len(cell_types)} cell types to the embedding dictionary.")
        
        cell_types = [x.replace('/', '|') for x in cell_types]
        
        file_config_path = file_path.replace('/*_processed.h5ad', '/') + config_path 
        if os.path.exists(file_config_path):
            with open(config_path, 'rb') as f:
                config = pickle.load(f)
            if not all([x in config.embedding_dict.keys() for x in cell_types]):
                emb = get_cell_type_embedding_by_llm(cell_types)
                for i in range(len(cell_types)):
                    embeddings_dict[cell_types[i]] = emb[i]
        else:
            emb = np.array(get_cell_type_embedding_by_llm(cell_types))
        for key, value in zip(cell_types, emb):
            embeddings_dict[key] = value
    
        del adata
        
    logging.info(f"Toal {len(embeddings_dict)} cell type embeddings are extracted.")
    
    with open(output_pkl, 'wb') as f:
        pickle.dump(embeddings_dict, f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_list', nargs='+', type=str, required=True)
    parser.add_argument('--output_pkl', type=str, required=True)
    args = parser.parse_args()
    
    extract_celltype_embedding(args.file_list, args.output_pkl)