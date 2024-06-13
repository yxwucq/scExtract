from typing import List, Optional
from tqdm import tqdm
import anndata as ad
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

def get_cell_type_embedding_by_llm(cell_types: List[str]) -> List[np.ndarray]:
    """
    Get cell type embeddings by using the OpenAI API.
    """
    api_key = ""
    azure_endpoint = ""
    client = AzureOpenAI(
    api_key=api_key, api_version="2024-02-01", azure_endpoint=azure_endpoint
    )
    response = client.embeddings.create(
        model="text-embedding-3-large", input=cell_types
    )
    emb = [x.embedding for x in response.data]
    return emb

def extract_celltype_embedding(file_list: str,
                               output_pkl: str):
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
        emb = np.array(get_cell_type_embedding_by_llm(cell_types))
        for key, value in zip(cell_types, emb):
            embeddings_dict[key] = value
    
    logging.info(f"Toal {len(embeddings_dict)} cell type embeddings are extracted.")
    
    with open(output_pkl, 'wb') as f:
        pickle.dump(embeddings_dict, f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_list', nargs='+', type=str, required=True)
    parser.add_argument('--output_pkl', type=str, required=True)
    args = parser.parse_args()
    
    extract_celltype_embedding(args.file_list, args.output_pkl)