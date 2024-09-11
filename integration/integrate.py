from typing import List, Optional
from tqdm import tqdm
import anndata as ad
import cellhint
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import scanpy as sc
from sklearn.neighbors import KernelDensity
import pandas as pd
import numpy as np
import logging
import pickle
import gc

from benchmark.benchmark import request_ols, ontology_pairwise_similarity
from auto_extract.agent import get_cell_type_embedding_by_llm

def merge_datasets(file_list: List[str],
                   downsample: Optional[bool] = False,
                   target_cells_per_label: Optional[int] = 1000
                   ) -> ad.AnnData:
    """
    Merge multiple datasets into one.
    """
    for file in file_list:
        adata = sc.read(file)
        print(f"Reading {file}, shape: {adata.shape}")
        sample_name = file.split('/')[-1].replace('_processed.h5ad', '')
        adata.obs['Dataset'] = sample_name.replace('_', ' ')
        
        if 'leiden' in adata.obs.columns:
            adata.obs['cell_type'] = adata.obs['leiden'].astype(str).str.replace('_', ' ').copy()
        elif 'louvain' in adata.obs.columns:
            adata.obs['cell_type'] = adata.obs['louvain'].astype(str).str.replace('_', ' ').copy()
        else:
            raise ValueError('No clustering result in the dataset')

        if downsample:
            adata = density_weighted_stratified_sampling(adata, target_cells_per_label)
            print(f"Downsample to: {adata.shape}")
            
        raw_adata = ad.AnnData(X=csr_matrix(adata.raw.X.copy()), obs=adata.obs.copy(), var=adata.var.copy())
        
        if 'adata_all' not in locals():
            adata_all = raw_adata.copy()
        else:
            adata_all = ad.concat([adata_all, raw_adata], join='outer')
            adata_all.obs_names_make_unique()
            
        del adata
        del raw_adata
        gc.collect()

    adata_all.obs = adata_all.obs.drop(columns = ['leiden', 'louvain'])
    del adata_all.var
    return adata_all

def preprocess_merged_dataset(adata_all: ad.AnnData) -> ad.AnnData:
    """
    Preprocess the merged dataset. Should be raw counts data with dataset and cell_type annotations.
    """
    adata_all.raw = adata_all
    sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)
    sc.pp.scale(adata_all)
    sc.tl.pca(adata_all)
    # sc.external.pp.harmony_integrate(adata_all, key = 'dataset')
    sc.pp.neighbors(adata_all, use_rep='X_pca')
    sc.tl.umap(adata_all)
    
    return adata_all

def create_prior_similarity_matrix(df_raw: pd.DataFrame,
                                   prior_method: str,
                                   config_path : str = 'config.ini',
                                   embedding_dict_path: Optional[str] = None) -> pd.DataFrame:
    
    similarity_matrix = np.zeros((len(df_raw), len(df_raw)))
    if prior_method == 'ontology':
        ontology_mapping_dict = {}

        for i in range(len(df_raw)):
            manual_cell_type = df_raw.columns[i].split('_')[1]
            manual_cell_type = manual_cell_type.strip('s')
            if manual_cell_type in ontology_mapping_dict.keys():
                continue
            obd_id, obd_label = request_ols(manual_cell_type)
            ontology_mapping_dict[manual_cell_type] = obd_id
        
        similarity_dict = {}
        for i in tqdm(range(len(df_raw))):
            for j in range(len(df_raw)):
                
                ct1 = df_raw.columns[i].split('_')[1].strip('s')
                ct2 = df_raw.index[j].split('_')[1].strip('s')
                
                if (ct1, ct2) in similarity_dict.keys():
                    similarity_matrix[i, j] = similarity_dict[(ct1, ct2)]
                elif (ct2, ct1) in similarity_dict.keys():
                    similarity_matrix[i, j] = similarity_dict[(ct2, ct1)]
                else:
                    obd_id1 = ontology_mapping_dict[ct1]
                    obd_id2 = ontology_mapping_dict[ct2]
                    if obd_id1 == 'None' or obd_id2 == 'None':
                        similarity_matrix[i, j] = 0
                    else:
                        similarity_matrix[i, j] = ontology_pairwise_similarity(ontology_mapping_dict[ct1], ontology_mapping_dict[ct2])
                    similarity_dict[(ct1, ct2)] = similarity_matrix[i, j]
        similarity_matrix_df = pd.DataFrame(similarity_matrix, index = df_raw.index, columns = df_raw.columns)
        return similarity_matrix_df
    
    elif prior_method == 'llm':
        ct = df_raw.index.str.split('_').str[1].tolist()
        emb = get_cell_type_embedding_by_llm(ct, config_path)
        emb = np.array(emb)
        emb_norm = np.linalg.norm(emb, axis = 1)
        similarity_matrix = np.dot(emb, emb.T) / np.outer(emb_norm, emb_norm)
        similarity_matrix_1 = (similarity_matrix - np.mean(similarity_matrix, axis=1)[:, np.newaxis]) / (1 - np.mean(similarity_matrix, axis=1)[:, np.newaxis])
        similarity_matrix_1[similarity_matrix_1 < 0] = 0
        similarity_matrix_2 = (similarity_matrix - np.mean(similarity_matrix, axis=1)[np.newaxis, :]) / (1 - np.mean(similarity_matrix, axis=1)[np.newaxis, :])
        similarity_matrix_2[similarity_matrix_2 < 0] = 0
        similarity_matrix = (similarity_matrix_1 + similarity_matrix_2) / 2
        similarity_matrix_df = pd.DataFrame(similarity_matrix, index = df_raw.index, columns = df_raw.index)
        return similarity_matrix_df
    
    elif prior_method == 'local':
        if embedding_dict_path is None:
            raise ValueError("Please provide the path to the cell type embedding dictionary.")
        with open(embedding_dict_path, 'rb') as f:
            embedding_dict = pickle.load(f)
            emb_list = [embedding_dict[x] for x in df_raw.index.str.split('_').str[1]]
            emb = np.array(emb_list)
            similarity_matrix = np.dot(emb, emb.T)
            similarity_matrix_1 = (similarity_matrix - np.mean(similarity_matrix, axis=1)[:, np.newaxis]) / (1 - np.mean(similarity_matrix, axis=1)[:, np.newaxis])
            similarity_matrix_1[similarity_matrix_1 < 0] = 0
            similarity_matrix_2 = (similarity_matrix - np.mean(similarity_matrix, axis=0)[np.newaxis, :]) / (1 - np.mean(similarity_matrix, axis=0)[np.newaxis, :])
            similarity_matrix_2[similarity_matrix_2 < 0] = 0
            similarity_matrix = (similarity_matrix_1 + similarity_matrix_2) / 2
            similarity_matrix_df = pd.DataFrame(similarity_matrix, index = df_raw.index, columns = df_raw.index)
            return similarity_matrix_df
        
    else:
        raise ValueError("Unknown prior method.")

def normalize_connectivities(df_raw: pd.DataFrame,
                            prior_similarity_matrix_df: pd.DataFrame,
                            prior_weight: float) -> pd.DataFrame:
        '''
        Normalize the raw connectivities matrix.
        '''
        
        assert df_raw.shape == prior_similarity_matrix_df.shape
        df_raw = df_raw.copy() # avoid changing the original matrix
        
        # combine raw connectivities matrix with prior similarity matrix
        for i in range(len(df_raw)):
            for j in range(len(df_raw)):
                df_raw.iloc[i, j] = (1-prior_weight) * df_raw.iloc[i, j] + prior_weight * prior_similarity_matrix_df.iloc[i, j]
                
                if df_raw.columns[i].split('_')[0] == df_raw.index[j].split('_')[0]:
                    df_raw.iloc[i, j] = 0 # set the same dataset to 0
      
        datasets = df_raw.columns.str.split('_').str[0].unique()

        # normalize the connectivities matrix for each dataset pair
        for dat1 in datasets:
            for dat2 in datasets:
                if dat1 == dat2:
                    continue
                idx1 = df_raw.columns.str.startswith(dat1)
                idx2 = df_raw.columns.str.startswith(dat2)
                # replace 0 with 1
                sum_idx1_idx2 = df_raw.loc[idx1, idx2].sum(axis=1)
                sum_idx1_idx2[sum_idx1_idx2 == 0] = 1
                df_raw.loc[idx1, idx2] = (df_raw.loc[idx1, idx2].values.T / sum_idx1_idx2.values).T.copy()
        
        return df_raw

def density_weighted_stratified_sampling(adata: ad.AnnData,
                                         target_cells_per_label: int = 1000,
                                         bandwidth: float = 0.5) -> ad.AnnData:
    """
    Subsample the dataset by density-weighted stratified
    """
    
    def compute_kde_weights(data, bandwidth):
        kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(data)
        log_dens = kde.score_samples(data)
        # Calculate weights
        weights = 1 / np.exp(log_dens)
        # Normalize weights
        weights /= np.sum(weights)
        return weights

    subsampled_data = []
    for label in adata.obs['cell_type'].unique():
        label_data = adata[adata.obs['cell_type'] == label].copy()
        
        if len(label_data) <= target_cells_per_label:
            subsampled_data.append(label_data)
        else:
            # Compute weights
            weights = compute_kde_weights(label_data.obsm['X_pca'], bandwidth)
            
            # Sample from the label data using the weights
            indices = np.random.choice(
                len(label_data), 
                size=target_cells_per_label, 
                replace=False, 
                p=weights
            )
            subsampled_data.append(label_data[indices])
    
    return sc.concat(subsampled_data)

def integrate_processed_datasets(file_list: List[str],
                                 method: str, # 'cellhint' or 'scExtract'
                                 output_path: str,
                                 config_path : str = 'config.ini',
                                 prior_weight: Optional[float] = 0.5,
                                 prior_method: Optional[str] = 'ontology', # 'ontology' or 'llm'
                                 alignment_path: Optional[str] = None,
                                 embedding_dict_path: Optional[str] = None,
                                 downsample: Optional[bool] = False,
                                 downsample_cells_per_label: Optional[int] = 1000,
                                 **kwargs) -> None:
    
    logging.basicConfig(level=logging.INFO)
    logging.info(f"Integrating {len(file_list)} datasets using {method} method.")
    adata_all = merge_datasets(file_list, downsample, downsample_cells_per_label)
    adata_all.write(output_path) # save the merged dataset
    
    logging.info(f"Merged dataset shape: {adata_all.shape}")
    adata_all = preprocess_merged_dataset(adata_all)
    # sc.pl.umap(adata_all, color = 'dataset')
    
    if method == 'cellhint':
        if alignment_path is None:
            alignment_path = output_path.replace('.h5ad', '.pkl')
        alignment = cellhint.harmonize(adata_all, dataset = 'Dataset', cell_type = 'cell_type', use_rep = 'X_pca', prior_path = embedding_dict_path, prior_weight = prior_weight, **kwargs)
        alignment.write(alignment_path)
    
    elif method == 'scExtract':
        adata_all.obs['dataset_cell_type'] = (adata_all.obs['Dataset'].astype(str) + '_' + adata_all.obs['cell_type'].astype(str)).astype('category')
    
        logging.info("Running PAGA for the merged dataset.")
        sc.pp.neighbors(adata_all, use_rep='X_pca')
        sc.tl.paga(adata_all, groups = 'dataset_cell_type')
        
        # # assign cell type colors
        # sc.pl.umap(adata_all, color = ['dataset', 'cell_type'], ncols = 1)
        
        # # mapping cell_type color to dataset_cell_type
        # adata_all.uns['dataset_cell_type_colors'] = []
        # for i, x in enumerate(adata_all.obs['dataset_cell_type'].cat.categories):
        #     adata_all.uns['dataset_cell_type_colors'].append(adata_all.uns['cell_type_colors'][adata_all.obs['cell_type'].cat.categories.get_loc(x.split('_')[1])])

        # create raw connectivities matrix
        df_raw = pd.DataFrame(adata_all.uns['paga']['connectivities'].toarray()).copy()
        df_raw.columns = adata_all.obs['dataset_cell_type'].cat.categories
        df_raw.index = adata_all.obs['dataset_cell_type'].cat.categories
        
        logging.info("Creating prior similarity matrix from automatic annotation.")
        prior_similarity_matrix_df = create_prior_similarity_matrix(df_raw, config_path, prior_method, embedding_dict_path)
        adata_all.uns['prior_similarity_matrix_df'] = prior_similarity_matrix_df.copy()
        adata_all.uns['raw_similarity_matrix_df'] = df_raw.copy()
        
        logging.info("Normalizing connectivities matrix.")
        df_norm = normalize_connectivities(df_raw, prior_similarity_matrix_df, prior_weight)
        
        # update connectivities in adata_all
        Tcsr = minimum_spanning_tree(csr_matrix(np.exp(-df_norm.values)))
        adata_all.uns['paga']['connectivities'] = Tcsr
        
    # save the integrated dataset
    adata_all.write(output_path)