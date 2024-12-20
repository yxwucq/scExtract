import configparser
import re
import mygene
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import logging
from tqdm import tqdm
from ..auto_extract.agent import Claude3, Openai, get_cell_type_embedding_by_llm
from ..auto_extract.parse_params import Params 

def major_vote_top_clusters(df, 
                        leiden_key='leiden', 
                        cell_type_key='cell_type_raw', 
                        dataset_key='Dataset', 
                        top_n=5):
    """
    Analyze cell clusters and return top cell types for each cluster.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame containing cluster and cell type information
    leiden_key : str, default='leiden'
        Column name for cluster identification
    cell_type_key : str, default='cell_type_raw'
        Column name for cell type information
    dataset_key : str, default='Dataset'
        Column name for dataset information
    top_n : int, default=5
        Number of top cell types to return for each cluster
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing analysis results
    """
    # Initialize empty list for storing results
    results = []
    
    # Iterate through each cluster
    for cluster in df[leiden_key].unique():
        # Get data for current cluster
        cluster_df = df[df[leiden_key] == cluster]
        
        # Calculate counts for each cell type
        cell_type_counts = cluster_df[cell_type_key].value_counts()
        
        # Get top N cell types by count
        top_cell_types = cell_type_counts.head(top_n)
        
        # Calculate unique Dataset support for each top cell type
        for cell_type in top_cell_types.index:
            cell_type_df = cluster_df[cluster_df[cell_type_key] == cell_type]
            dataset_support = len(cell_type_df[dataset_key].unique())
            
            # Add results to list
            results.append({
                'cluster': cluster,
                'cell_type': cell_type,
                'cell_count': cell_type_counts[cell_type],
                'dataset_support': dataset_support,
                'percentage': (cell_type_counts[cell_type] / len(cluster_df) * 100).round(2)
            })
    
    # Convert results to DataFrame
    result_df = pd.DataFrame(results)
    
    # Sort by cluster and cell count
    try:
        result_df['cluster'] = result_df['cluster'].astype(int)
    except ValueError:
        pass
    result_df = result_df.sort_values(['cluster', 'cell_count'], ascending=[True, False])
    
    return result_df

def majority_voting_cluster_annotation(cell_type_by_leiden_df: pd.DataFrame,
                                       config_path: str = 'config.ini'):
                                       
    config = configparser.ConfigParser()
    config.read(config_path)

    if 'openai' in config['API']['TYPE']:
        claude_agent = Openai(pdf_path=None, config_path=config_path)
    elif 'claude' in config['API']['TYPE']:
        claude_agent = Claude3(pdf_path=None, config_path=config_path)
    else:
        raise ValueError(f"Model {config['API']['MODEL']} not supported.")
    
    # split df into chunks of 50 rows
    content_list = []
    for i in range(0, len(cell_type_by_leiden_df), 50):
        end_idx = min(i+50, len(cell_type_by_leiden_df))
        content_list.append(cell_type_by_leiden_df.iloc[i:end_idx, :].to_string())
    
    params = Params(config_path)
    
    all_annotation_dict = {}
    all_reasoning = ""
    print('Majority voting annotation in progress...')
    for content in content_list:
        major_voting_prompt = params.get_tool_prompt('MAJOR_VOTE_ANNOTATION_PROMPT') + content
        model_response = claude_agent._tool_retrieve(messages=[{"role": "user", "content": major_voting_prompt}])
        annotation_dict = params.parse_annotation_response(model_response, simple_annotation=True)
        all_annotation_dict.update(annotation_dict)
        all_reasoning += model_response + '\n'
    
    return all_annotation_dict, all_reasoning

def plot_multiple_cluster_composition(df, clusters=None, figsize=(15, 3.5), min_percentage=1, ncols=3):
    """
    Create visualizations of cell type composition for multiple clusters in a grid layout.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing columns: 'cluster', 'cell_type', 'cell_count', 
        'dataset_support', 'percentage'
    clusters : list, optional
        List of cluster IDs to plot. If None, plot all clusters
    figsize : tuple, default=(15, 3.5)
        Figure size per row (width, height)
    min_percentage : float, default=1
        Minimum percentage to display in the plots
    ncols : int, default=3
        Number of clusters to display per row
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    # Determine which clusters to plot
    if clusters is None:
        clusters = sorted(df['cluster'].unique())
    else:
        clusters = [c for c in clusters if c in df['cluster'].unique()]
    
    # Calculate number of rows needed
    nrows = (len(clusters) + ncols - 1) // ncols
    
    # Adjust figure size based on number of rows
    adjusted_height = figsize[1] * nrows
    fig = plt.figure(figsize=(figsize[0], adjusted_height))
    
    # Create grid for subplots with increased column padding
    gs = fig.add_gridspec(nrows, ncols * 3, 
                         width_ratios=np.tile([2, 1, 1], ncols),  # Increased last ratio for padding
                         hspace=0.5,  # Vertical spacing between rows
                         wspace=0.8)  # Increased horizontal spacing between columns
    
    for idx, cluster_id in enumerate(clusters):
        # Calculate row and column position
        row = idx // ncols
        col = (idx % ncols) * 3  # *3 because we have 3 columns per cluster
        
        # Get data for current cluster
        plot_df = df[df['cluster'] == cluster_id].copy()
        plot_df = plot_df[plot_df['percentage'] >= min_percentage]
        
        # Create subplots for this cluster with different widths
        ax1 = fig.add_subplot(gs[row, col])    # Takes 2 width units
        ax2 = fig.add_subplot(gs[row, col+1])  # Takes 1 width unit
        
        # Color palette
        colors = sns.color_palette("husl", n_colors=len(plot_df))
        
        # First subplot: Cell counts and percentages
        bars = ax1.barh(range(len(plot_df)), plot_df['percentage'], color=colors)
        
        # Add cell count annotations inside the bars
        for i, bar in enumerate(bars):
            width = bar.get_width()
            count = plot_df['cell_count'].iloc[i]
            count_text = f'{count:,}'
            text_position = max(width * 0.01, 0.5)
            
            ax1.text(text_position, i,
                    f'{count_text} ({plot_df["percentage"].iloc[i]:.1f}%)',
                    va='center',
                    color='black')
        
        # Second subplot: Dataset support
        bars2 = ax2.barh(range(len(plot_df)), plot_df['dataset_support'],
                        color=colors, alpha=0.7)
        
        # Add dataset support annotations
        for i, bar in enumerate(bars2):
            width = bar.get_width()
            dataset_count = plot_df['dataset_support'].iloc[i]
            ax2.text(width + 0.1, i,
                    str(dataset_count),
                    va='center')
        
        # Customize both subplots
        for ax in [ax1, ax2]:
            ax.set_yticks(range(len(plot_df)))
            ax.grid(True, axis='x', alpha=0.3)
        
        # Set y-labels
        ax1.set_yticklabels(plot_df['cell_type'])
        ax2.set_yticklabels([])
        
        # Set titles and labels
        ax1.set_title(f"Cluster {cluster_id}")
        ax1.set_xlabel('Percentage of Cells')
        ax2.set_xlabel('N Datasets')
        
        # Adjust x-axis limits for dataset support plot
        ax2.set_xlim(0, max(df['dataset_support']) * 1.15)  # Add 15% padding
    
    # Adjust layout with specific padding
    fig.tight_layout()  # Increase padding around the entire figure
    
    return fig

def convert_ensembl_to_symbol(adata: ad.AnnData) -> ad.AnnData:
    """
    Convert Ensembl IDs to gene symbols.
    """
    
    mg = mygene.MyGeneInfo()
    shape_before = adata.shape
    
    logging.info('Finding gene symbols for Ensembl IDs...')
    adata.var.index = adata.var.index.str.split('.').str[0]
    
    retry = 0
    while retry < 3:
        try:
            query_dict_list = mg.querymany(adata.var.index, scopes='ensembl.gene', fields='symbol')
            break
        except:
            retry += 1
            logging.warning('Failed to connect to MyGeneInfo. Retrying...')
    
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