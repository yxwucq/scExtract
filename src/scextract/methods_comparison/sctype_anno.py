# use python implementation at https://github.com/kris-nader/sc-type-py

import urllib.request
import scanpy as sc
import numpy as np
import pandas as pd

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale, MinMaxScaler
import requests
import xml.etree.ElementTree as ET
import scanpy as sc
from collections import defaultdict
import concurrent.futures
import multiprocessing
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import configparser
import re

from ..auto_extract.agent import Claude3, Openai, get_cell_type_embedding_by_llm
from ..auto_extract.parse_params import Params 

# From https://github.com/kris-nader/sc-type-py
def gene_sets_prepare(path_to_db_file, cell_type):
    # Read data from Excel file
    cell_markers = pd.read_excel(path_to_db_file)
    # Filter by cell_type
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    # Preprocess geneSymbolmore1 and geneSymbolmore2 columns
    #for col in ['geneSymbolmore1', 'geneSymbolmore2']:
    #    cell_markers[col] = cell_markers[col].str.replace(" ", "").str.upper()
    for col in ['geneSymbolmore1', 'geneSymbolmore2']:
        if col in cell_markers.columns:
            cell_markers[col] = cell_markers[col].fillna('').str.replace(" ", "").str.upper()
    # Stack and drop duplicates to get unique gene names
    gene_names = pd.concat([cell_markers['geneSymbolmore1'], cell_markers['geneSymbolmore2']]).str.split(',', expand=True).stack().drop_duplicates().reset_index(drop=True)
    # gene_names = gene_names[gene_names.str.strip() != '']
    gene_names = gene_names[gene_names != 'None'].unique()
    # Get approved symbols for gene names
    res = get_gene_symbols(set(gene_names))
    res = dict(zip(res['Gene'], res['Symbol']))
    # Process gene symbols
    for col in ['geneSymbolmore1', 'geneSymbolmore2']:
        cell_markers[col] = cell_markers[col].apply(lambda row: process_gene_symbols(row, res)).str.replace("///", ",").str.replace(" ", "")
    # Group by cellName and create dictionaries of gene sets
    gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(lambda x: list(set(','.join(x).split(',')))).to_dict()
    gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(lambda x: list(set(','.join(x).split(',')))).to_dict()
    return {'gs_positive': gs_positive, 'gs_negative': gs_negative}

# From https://github.com/kris-nader/sc-type-py
def process_gene_symbols(gene_symbols, res):
    if pd.isnull(gene_symbols):
        return ""
    markers_all = gene_symbols.upper().split(',')
    markers_all = [marker.strip().upper() for marker in markers_all if marker.strip().upper() not in ['NA', '']]
    markers_all = sorted(markers_all)
    if len(markers_all) > 0:
        markers_all = [res.get(marker) for marker in markers_all]
        markers_all = [symbol for symbol in markers_all if symbol is not None]
        markers_all = list(set(markers_all))
        return ','.join(markers_all)
    else:
        return ""
    
# From https://github.com/kris-nader/sc-type-py
def get_gene_symbols(genes):
    data = {"Gene": [], "Symbol": []}
    for gene in genes:
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        url = f"https://rest.genenames.org/fetch/symbol/{gene}"
        response = session.get(url)
        if response.status_code == 200:
            root = ET.fromstring(response.content)
            result_elem = root.find("result")
            if result_elem is not None and result_elem.get("numFound") == "0":
                url = f"https://rest.genenames.org/search/alias_symbol/{gene}"
                response = session.get(url)
                if response.status_code == 200:
                    root = ET.fromstring(response.content)
                    result_elem = root.find("result")
                    if result_elem is not None and result_elem.get("numFound") != "0":
                        symbols = [doc.find('str[@name="symbol"]').text for doc in root.findall('.//doc')]
                        data["Gene"].append(gene)
                        data["Symbol"].append(','.join(symbols))
                    elif result_elem is not None and result_elem.get("numFound") == "0":
                        url = f"https://rest.genenames.org/search/prev_symbol/{gene}"
                        response = session.get(url)
                        if response.status_code == 200:
                            root = ET.fromstring(response.content)
                            result_elem = root.find("result")
                            if result_elem is not None and result_elem.get("numFound") != "0":
                                symbol_element = root.find('.//str[@name="symbol"]').text
                                data["Gene"].append(gene)
                                data["Symbol"].append(symbol_element)
                else:
                    print(f"Failed to retrieve data for gene {gene}. Status code:", response.status_code)
            else:
                symbol_element = root.find('.//str[@name="symbol"]').text
                data["Gene"].append(gene)
                data["Symbol"].append(gene)
        else:
            print(f"Failed to retrieve data for gene {gene}. Status code:", response.status_code)
    df = pd.DataFrame(data)
    return df

# From https://github.com/kris-nader/sc-type-py
def sctype_score(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
    marker_stat = defaultdict(int, {gene: sum(gene in genes for genes in gs.values()) for gene in set(gene for genes in gs.values() for gene in genes)})
    marker_sensitivity = pd.DataFrame({'gene_': list(marker_stat.keys()), 'score_marker_sensitivity': list(marker_stat.values())})
    # Rescaling the score_marker_sensitivity column
    # grab minimum and maximum
    min_value=1
    max_value= len(gs)
    # Apply the formula to the column
    marker_sensitivity['score_marker_sensitivity'] = 1 - (marker_sensitivity['score_marker_sensitivity'] - min_value) / (max_value - min_value)
    # Convert gene names to Uppercase
    if gene_names_to_uppercase:
        scRNAseqData.index = scRNAseqData.index.str.upper()
    # Subselect genes only found in data
    names_gs_cp = list(gs.keys())
    names_gs_2_cp = list(gs2.keys())
    gs_ = {key: [gene for gene in scRNAseqData.index if gene in gs[key]] for key in gs}
    gs2_ = {key: [gene for gene in scRNAseqData.index if gene in gs2[key]] for key in gs2}
    gs__ = dict(zip(names_gs_cp, gs_.values()))
    gs2__ = dict(zip(names_gs_2_cp, gs2_.values()))
    
    all_genes = set().union(*[set(genes) for genes in gs__.values()])
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity['gene_'].isin(all_genes)]
    
    # Z-scale if not
    if not scaled:
        Z = scale(scRNAseqData.T).T
    else:
        Z = scRNAseqData
    # Multiply by marker sensitivity
    for _, row in cell_markers_genes_score.iterrows():
        Z.loc[row['gene_']] *= row['score_marker_sensitivity']
    marker_genes = set().union(*[set(genes) for genes in gs__.values()], 
                              *[set(genes) for genes in gs2__.values()])
    Z = Z.loc[list(marker_genes)]
    # Combine scores
    es = pd.DataFrame(index=gs__.keys(),columns=Z.columns,data=np.zeros((len(gs__), Z.shape[1])))
    for gss_, genes in gs__.items():
        for j in range(Z.shape[1]):
            gs_z = Z.loc[genes, Z.columns[j]]
            gz_2 = Z.loc[gs2__[gss_], Z.columns[j]] * -1 if gs2__ and gss_ in gs2__ else pd.Series(dtype=np.float64)
            sum_t1 = np.sum(gs_z) / np.sqrt(len(gs_z))
            sum_t2 = np.sum(gz_2) / np.sqrt(len(gz_2)) if not gz_2.empty else 0
            if pd.isna(sum_t2):
                sum_t2 = 0
            es.loc[gss_, Z.columns[j]] = sum_t1 + sum_t2
    es = es.dropna(how='all')
    return es

# From https://github.com/kris-nader/sc-type-py
def process_cluster(cluster,adata,es_max,clustering):
    cluster_data = es_max.loc[:, adata.obs.index[adata.obs[clustering] == cluster]]
    es_max_cl = cluster_data.sum(axis=1).sort_values(ascending=False)
    top_scores = es_max_cl.head(10)
    ncells = sum(adata.obs[clustering] == cluster)
    return pd.DataFrame({
        'cluster': cluster,
        'type': top_scores.index,
        'scores': top_scores.values,
        'ncells': ncells
    })

def get_best_model(pdf_path: str, config_path: str = 'config.ini') -> str:
    config = configparser.ConfigParser()
    config.read(config_path)

    if 'openai' in config['API']['TYPE']:
        claude_agent = Openai(pdf_path, config_path)
    elif 'claude' in config['API']['TYPE']:
        claude_agent = Claude3(pdf_path, config_path)
    else:
        raise ValueError(f"Model {config['API']['MODEL']} not supported.")
    
    df_models = pd.read_excel('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx')
    model_list_str = ','.join(list(df_models['tissueType'].unique()))

    params = Params(config_path)
    claude_agent.initiate_propmt()
    
    choose_model_prompt = params.get_prompt('CHOOSE_MODEL_PROMPT').replace('please choose the model that best fits your data:',
                    f"please choose the model that best fits your data: {model_list_str}")
    model_response = claude_agent.chat(choose_model_prompt)
    print(model_response)
    # parse response
    # model: 'Immune_All_Low.pkl'
    model_name = re.search(r"model: (.*)", model_response).group(1)
    model_name = model_name.replace("'", "")
    
    return model_name

def add_sctype_annotation(pdf_path: str,
                              adata_path: str,
                              config_path: str = 'config.ini',
                              key_added: str = 'sctype_annotation',
                              output_path: str = None,
                                ):

    model_name = get_best_model(pdf_path, config_path)
    model_name = model_name.replace('.pkl', '')

    adata = ad.read_h5ad(adata_path)
    scaled_data = pd.DataFrame(adata.X)
    # change column indexes
    scaled_data.columns =adata.var_names
    # Change the row indexes
    scaled_data.index = adata.obs_names
    scaled_data=scaled_data.T

    scRNAseqData=scaled_data
    gs_list=gene_sets_prepare(path_to_db_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
                              cell_type=model_name)
    es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
    
    group_by_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'
    
    unique_clusters = adata.obs[group_by_key].unique()
    # Apply the function to each unique cluster and combine the results into a DataFrame
    cL_results = pd.concat([process_cluster(cluster,adata,es_max,group_by_key) for cluster in unique_clusters])

    # Group by cluster and select the top row based on scores
    sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

    # Set low-confidence clusters to "Unknown"
    sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

    # Iterate over unique clusters
    adata.obs[key_added] = ""
    for cluster in sctype_scores['cluster'].unique():
        # Filter sctype_scores for the current cluster
        cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
        # Get the type for the current cluster
        cl_type_value = cl_type['type'].iloc[0]
        # Update key_added in pbmc.obs for cells belonging to the current cluster
        adata.obs.loc[adata.obs[group_by_key] == cluster, key_added] = cl_type_value
    
    print(f"Annotation added to adata.obs[{key_added}]")
    print(f"Successfully annotated {len(sctype_scores)} clusters")
    
    # save adata
    if output_path is not None:
        adata.write_h5ad(output_path)
    else:
        adata.write_h5ad(adata_path)

# if __name__ == '__main__':
#     add_sctype_annotation(
#         pdf_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/sample14/raw_data/sample14.pdf',
#         adata_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/sample14/sample14_deepseek_v2_5_with_celltypist.h5ad',
#         config_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/config.ini',
#         key_added='sctype_annotation',
#         output_path='/home/wu/datb1/AutoExtractSingleCell/01.benchmark_datasets/sample14/sample14_deepseek_v2_5_with_sctype.h5ad'
#     )