import mygene
import pandas as pd
import scanpy as sc
import anndata as ad
import logging
from tqdm import tqdm

def convert_ensembl_to_symbol(adata: ad.AnnData) -> ad.AnnData:
    """
    Convert Ensembl IDs to gene symbols.
    """
    
    mg = mygene.MyGeneInfo()
    shape_before = adata.shape
    
    logging.info('Finding gene symbols for Ensembl IDs...')
    query_dict_list = mg.querymany(adata.var.index, scopes='ensembl.gene', fields='symbol')
    
    logging.info('Converting Ensembl IDs to gene symbols...')
    convert_dict = {}
    for x in query_dict_list:
        if x['query'] in convert_dict.keys():
            print('duplicated Ensembl IDs: '+x['query'])
        if 'symbol' in x:
            convert_dict[x['query']] = x['symbol']
        else:
            convert_dict[x['query']] = x['query']
    
    for ind, rows in tqdm(adata.var.iterrows()):
        if ind in convert_dict.keys():
            adata.var.loc[ind, 'gene_symbol'] = convert_dict[ind]
        else:
            adata.var.loc[ind, 'gene_symbol'] = ind
    
    adata_extracted = adata[:, adata.var['gene_symbol'] != adata.var.index].copy()
    adata_extracted.var = adata_extracted.var.reset_index(names='Ensembl_ID').set_index('gene_symbol', drop=True).copy()
    
    shape_after = adata_extracted.shape
    
    logging.info(f"Successfully Converted {shape_after[1]} Ensembl IDs to gene symbols.")
    logging.info(f"{shape_before[1] - shape_after[1]} Ensembl IDs were not converted.")
    
    return adata_extracted 

if __name__ == '__main__':
    adata = sc.read_h5ad('/home/wu/datb1/AutoExtractSingleCell/sample2/raw_data/raw_data.h5ad')
    adata_test = adata[:, adata.var.iloc[:10, :].index].copy()
    adata_test

    adata_extracted = convert_ensembl_to_symbol(adata_test)