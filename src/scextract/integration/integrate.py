from typing import List, Optional
from functools import wraps
import time
from tqdm import tqdm
import anndata as ad
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import scanpy as sc
from sklearn.neighbors import KernelDensity
import pandas as pd
import numpy as np
import logging
import pickle
import gc

from ..benchmark.benchmark import request_ols, ontology_pairwise_similarity
from ..auto_extract.agent import get_cell_type_embedding_by_llm

def time_print(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        print(f"Starting {func.__name__}...")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        duration_min = (end_time - start_time) / 60
        print(f"Finished {func.__name__} in {duration_min:.2f} minutes")
        return result
    return wrapper

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
        if 'Dataset' not in adata.obs.columns:
            adata.obs['Dataset'] = sample_name.replace('_', ' ')
        
        if 'cell_type' in adata.obs.columns:
            adata.obs['cell_type'] = adata.obs['cell_type'].astype(str).str.replace('_', ' ').copy()
        elif 'leiden' in adata.obs.columns:
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

    if 'leiden' in adata_all.obs.columns or 'louvain' in adata_all.obs.columns:
        adata_all.obs = adata_all.obs.drop(columns = ['leiden', 'louvain'], errors='ignore')
    del adata_all.var
    return adata_all

def preprocess_merged_dataset(adata_all: ad.AnnData, 
                              dimred=50) -> ad.AnnData:
    """
    Preprocess the merged dataset. Should be raw counts data with dataset and cell_type annotations.
    """
    adata_all.raw = adata_all
    sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)
    sc.pp.scale(adata_all)
    sc.pp.pca(adata_all, n_comps=dimred)
    # sc.pp.neighbors(adata_all, use_rep='X_pca')
    # sc.tl.umap(adata_all)
    
    return adata_all

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
                                 method: str,
                                 output_path: str,
                                 config_path : str = 'config.ini',
                                 alignment_path: Optional[str] = None,
                                 embedding_dict_path: Optional[str] = None,
                                 downsample: Optional[bool] = False,
                                 downsample_cells_per_label: Optional[int] = 1000,
                                 search_factor: int = 5,
                                 approx: bool = False,
                                 use_gpu: bool = False,
                                 batch_size: int = 5000,
                                 dimred: int = 100,
                                 use_pct: bool = False,
                                 **kwargs) -> None:
    
    """
    Integrate multiple processed datasets
    
    Parameters
    ----------
    file_list : List[str]
        List of paths to the processed data in AnnData format.
    method : str
        Method to use for integration. Support 'cellhint', 'cellhint_prior', 'scanorama_prior' and 'scExtract'.
    output_path : str
        Path to save the output file. If not specified, the input file will be overwritten.
    config_path : str
        System config file path.
    prior_weight : float
        Weight of the prior similarity matrix from automatic annotation.
    alignment_path : str
        Path to the output alignment file.
    embedding_dict_path : str
        Path to the cell type embedding dictionary. Only used for local. Can be generated by extract_celltype_embedding
    downsample : bool
        Whether to downsample the cells.
    downsample_cells_per_label : int
        Number of cells to downsample per label.
    search_factor : int
        Search factor for scanorama_prior, only vaiable for approx = True in scanorama_prior.
    approx : bool
        Whether to use approximate search in scanorama_prior and scanorama.
    use_gpu : bool
        Whether to use GPU for scanorama_prior.
    batch_size : int
        Batch size for scanorama_prior.
    dimred : int
        Number of dimensions for PCA.
    use_pct : bool
        Whether to use percent of cells for cellhint.
    **kwargs : dict
        Additional parameters for the scanorama_prior method.
    """
    
    # check if dependencies are installed
    if method == 'cellhint':
        try:
            import cellhint
            cellhint.harmonize = time_print(cellhint.harmonize)
        except ImportError:
            raise ImportError("Please install cellhint to use the cellhint method.")
    elif method == 'scanorama':
        try:
            import scanorama
            scanorama.scanorama.integrate_scanpy = time_print(scanorama.scanorama.integrate_scanpy)
        except ImportError:
            raise ImportError("Please install scanorama to use the scanorama method.")
    elif method == 'cellhint_prior':
        try:
            import cellhint_prior
            cellhint_prior.harmonize = time_print(cellhint_prior.harmonize)
        except ImportError:
            raise ImportError("Please install cellhint_prior to use the cellhint_prior method.")
    elif method == 'scanorama_prior':
        try:
            import scanorama_prior
            scanorama_prior.scanorama.integrate_scanpy = time_print(scanorama_prior.scanorama.integrate_scanpy)
        except ImportError:
            raise ImportError("Please install scanorama_prior to use the scanorama_prior method.")    
    elif method == 'scExtract':
        try:
            import cellhint_prior
            cellhint_prior.harmonize = time_print(cellhint_prior.harmonize)
        except ImportError:
            raise ImportError("Please install cellhint_prior to use the cellhint method.")
        try:
            import scanorama_prior
            scanorama_prior.scanorama.integrate_scanpy = time_print(scanorama_prior.scanorama.integrate_scanpy)
        except ImportError:
            raise ImportError("Please install scanorama_prior to use the scExtract method.")
    
    logging.basicConfig(level=logging.INFO)
    logging.info(f"Integrating {len(file_list)} datasets using {method} method.")
    
    if method == 'scanorama_prior':
        print("Remind: Scanorama_prior requires the harmonized cell type")
        
        assert len(file_list) == 1, "Scanorama_prior only supports merged dataset. \
            First merge the datasets using --method cellhint_prior."
        
        adata_all = sc.read_h5ad(file_list[0])
        # Use log1p normalized data for scanorama
        adata_all = adata_all.raw.to_adata()
        adata_all.raw = adata_all
        sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)

        logging.info(f"Merged dataset shape: {adata_all.shape}")
        
        harmonized_celltype_list = adata_all.obs[f"cell_type"].unique().tolist()
        
        with open(embedding_dict_path, 'rb') as f:
            embedding_dict = pickle.load(f)
            
        emb_list = [embedding_dict[x] for x in harmonized_celltype_list]
        harmonized_celltype_embedding_similarities = np.dot(np.array(emb_list), np.array(emb_list).T)
        harmonized_celltype_embedding_similarities_df = pd.DataFrame(harmonized_celltype_embedding_similarities, index = harmonized_celltype_list, columns = harmonized_celltype_list)

        # split data into batches
        batches = adata_all.obs['Dataset'].cat.categories.tolist()
        adatas = []
        for batch in batches:
            adatas.append(adata_all[adata_all.obs['Dataset'] == batch,].copy())
        
        del adata_all
        gc.collect()
        
        scanorama_prior.scanorama.integrate_scanpy(adatas,
                                         type_similarity_matrix=harmonized_celltype_embedding_similarities_df,
                                         search_factor = search_factor,
                                         approx = approx,
                                         use_gpu = use_gpu,
                                         batch_size = batch_size,
                                         dimred = dimred,
                                         **kwargs
                                         )
        
        adata_all = ad.concat(adatas, join='outer')
        adata_all.obsm['X_scanorama_prior'] = adata_all.obsm['X_scanorama'].copy()
        
        del adatas
        del adata_all.obsm['X_scanorama']
        gc.collect()
        
        sc.pp.neighbors(adata_all, n_neighbors=30, use_rep='X_scanorama_prior')
        sc.tl.umap(adata_all)
        
        adata_all.write(output_path)
    
    elif method == 'scanorama':
        assert len(file_list) == 1, "Scanorama_prior only supports merged dataset. \
            First merge the datasets using --method cellhint_prior."
        
        adata_all = sc.read_h5ad(file_list[0])
        # Use log1p normalized data for scanorama
        adata_all = adata_all.raw.to_adata()
        adata_all.raw = adata_all
        sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)

        logging.info(f"Merged dataset shape: {adata_all.shape}")
        
        # split data into batches
        batches = adata_all.obs['Dataset'].cat.categories.tolist()
        adatas = []
        for batch in batches:
            adatas.append(adata_all[adata_all.obs['Dataset'] == batch,].copy())
        
        scanorama.scanorama.integrate_scanpy(adatas,
                                         approx = approx,
                                         batch_size = batch_size,
                                         dimred = dimred,
                                         **kwargs
                                         )

        adata_all = ad.concat(adatas, join='outer')
        
        del adatas
        gc.collect()

        sc.pp.neighbors(adata_all, n_neighbors=30, use_rep='X_scanorama')
        sc.tl.umap(adata_all)
        
        adata_all.write(output_path)

    else:
        if len(file_list) > 1:
            adata_all = merge_datasets(file_list, downsample, downsample_cells_per_label)
            adata_all.obs['cell_type_raw'] = adata_all.obs['cell_type'].copy()
            adata_all.write(output_path) # save the merged dataset
            adata_all = preprocess_merged_dataset(adata_all, dimred=dimred)

        else:
            adata_all = sc.read_h5ad(file_list[0]) # default to read the merged dataset
            adata_all = adata_all.raw.to_adata()
            adata_all.raw = adata_all
            sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)
            sc.pp.scale(adata_all)
            sc.pp.pca(adata_all, n_comps=dimred)
            
            adata_all.obs['cell_type_raw'] = adata_all.obs['cell_type'].copy()
            adata_all.write(output_path) # save the merged dataset
        
        logging.info(f"Merged dataset shape: {adata_all.shape}")
        
        if method == 'scExtract':
            import scanorama_prior
            import cellhint_prior
            
            if alignment_path is None:
                alignment_path = output_path.replace('.h5ad', '.pkl')
            
            harmonized_celltype_list = adata_all.obs[f"cell_type"].unique().tolist()
            harmonized_celltype_embedding_list = get_cell_type_embedding_by_llm(harmonized_celltype_list, config_path = config_path)
            embedding_dict = dict(zip(harmonized_celltype_list, harmonized_celltype_embedding_list))
            
            alignment = cellhint_prior.harmonize(adata_all, dataset = 'Dataset', cell_type = 'cell_type', use_rep = 'X_pca', embedding_dict=embedding_dict, **kwargs)
            print("Harmonized cell type stored in 'harmonized_cellhint_prior' column.")
            alignment.write(alignment_path)
            adata_all.obs[f"harmonized_cellhint_prior"] = alignment.reannotation.loc[adata_all.obs_names, ['reannotation']].copy()
            adata_all.obs[f"cell_type"] = adata_all.obs[f"harmonized_cellhint_prior"].apply(remove_none_type)
            
            harmonized_celltype_list = adata_all.obs[f"cell_type"].unique().tolist()
            
            harmonized_celltype_embedding_list = get_cell_type_embedding_by_llm(harmonized_celltype_list, config_path = config_path)
            harmonized_celltype_embedding_similarities = np.dot(np.array(harmonized_celltype_embedding_list), np.array(harmonized_celltype_embedding_list).T)
            harmonized_celltype_embedding_similarities_df = pd.DataFrame(harmonized_celltype_embedding_similarities, index = harmonized_celltype_list, columns = harmonized_celltype_list)

            # split data into batches
            # Use log1p normalized data for scanorama
            adata_all = adata_all.raw.to_adata()
            adata_all.raw = adata_all        
            sc.pp.highly_variable_genes(adata_all, batch_key = 'Dataset', subset = True)

            batches = adata_all.obs['Dataset'].cat.categories.tolist()
            adatas = []
            for batch in batches:
                adatas.append(adata_all[adata_all.obs['Dataset'] == batch,].copy())
            
            del adata_all
            gc.collect()
            
            scanorama_prior.scanorama.integrate_scanpy(adatas,
                                            type_similarity_matrix=harmonized_celltype_embedding_similarities_df,
                                            search_factor = search_factor,
                                            use_gpu = use_gpu,
                                            batch_size = batch_size,
                                            dimred = dimred,
                                            **kwargs
                                            )
            
            adata_all = ad.concat(adatas, join='outer')
            adata_all.obsm['X_scanorama_prior'] = adata_all.obsm['X_scanorama'].copy()
            
            del adatas
            del adata_all.obsm['X_scanorama']
            gc.collect()
            
            sc.pp.neighbors(adata_all, n_neighbors=30, use_rep='X_scanorama_prior')
            sc.tl.umap(adata_all)
        
        elif method == 'cellhint':
            if alignment_path is None:
                alignment_path = output_path.replace('.h5ad', '.pkl')
            alignment = cellhint.harmonize(adata_all, dataset = 'Dataset', cell_type = 'cell_type', use_rep = 'X_pca', use_pct=use_pct, **kwargs)
            print("Harmonized cell type stored in 'harmonized_cellhint' column.")
            adata_all.obs[f"harmonized_cellhint"] = alignment.reannotation.loc[adata_all.obs_names, ['reannotation']].copy()
            
            adata_all.obs[f"cell_type"] = adata_all.obs[f"harmonized_cellhint"].apply(remove_none_type)
            cellhint.integrate(adata_all, batch = 'Dataset', cell_type = f"cell_type", use_rep='X_pca')
            
            alignment.write(alignment_path)
        
        elif method == 'cellhint_prior':        
            if alignment_path is None:
                alignment_path = output_path.replace('.h5ad', '.pkl')
            with open(embedding_dict_path, 'rb') as f:
                embedding_dict = pickle.load(f)
            alignment = cellhint_prior.harmonize(adata_all, dataset = 'Dataset', cell_type = 'cell_type', use_rep = 'X_pca', use_pct=use_pct, embedding_dict=embedding_dict, **kwargs)
            print("Harmonized cell type stored in 'harmonized_cellhint_prior' column.")
            adata_all.obs[f"harmonized_cellhint_prior"] = alignment.reannotation.loc[adata_all.obs_names, ['reannotation']].copy()
            adata_all.obs[f"cell_type"] = adata_all.obs[f"harmonized_cellhint_prior"].apply(remove_none_type)
            cellhint_prior.integrate(adata_all, batch = 'Dataset', cell_type = f"cell_type", use_rep='X_pca')

            alignment.write(alignment_path)
            
        else:
            raise ValueError("Unknown method.")        
            
        # save the integrated dataset
        adata_all.write(output_path)

def remove_none_type(x):
    x = x.replace('NONE = ', '')
    x = x.replace(' = NONE', '')
    x = x.replace('UNRESOLVED = ', '')
    x = x.replace(' = UNRESOLVED', '')
    return x

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
        emb = get_cell_type_embedding_by_llm(ct, config_path = config_path)
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