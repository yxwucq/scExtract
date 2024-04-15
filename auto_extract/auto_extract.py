import anndata as ad
import logging
import os
import warnings
import time

from pyfiglet import Figlet
import colorama
from termcolor import colored

from .claude3 import Claude3, Openai
from .preprocess import filter, preprocess, clustering
from .annotation import get_marker_genes, annotate, query_datasets
from .parse_params import Params 
from .config import Config

from utils.uitls import convert_ensembl_to_symbol

def auto_extract(adata_path: str, 
                 pdf_path: str, 
                 output_dir: str,
                 output_name: str = 'processed.h5ad') -> None:
    """
    Extracts and processes single-cell data from literature.
    """
    
    logging.basicConfig(level=logging.INFO)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.exists(os.path.join(output_dir, 'auto_extract.log')):
        os.remove(os.path.join(output_dir, 'auto_extract.log'))
    file_handler = logging.FileHandler(os.path.join(output_dir, 'auto_extract.log'))
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)
    
    f = Figlet(font='slant')
    
    # Print "AutoExtract" in fancy font
    logging.info('\n'+f.renderText('scExtract')+ '\n')

    colorama.init()

    # Load AnnData object
    logging.info(f'Loading AnnData object from {adata_path}')
    adata = ad.read_h5ad(adata_path)
    
    if adata.shape[0] > 50000:
        logging.warning(colored(f"Large dataset detected. The dataset contains {adata.shape[0]} cells. \
                            The result extracted by Claude3 may not be accurate. Please consider manually \
                                subsetting the data to a smaller size.", color='light_red'))
    
    # Check if Ensembl IDs are present in AnnData object
    for genes in adata.var.index[:10]:
        if genes.startswith('ENSG'):
            logging.info('Ensembl IDs detected in AnnData object.')
            adata = convert_ensembl_to_symbol(adata)
            adata.var_names_make_unique()
            break
    
    logging.info(colored('1. Extracting information from literature', color='cyan', attrs=['bold']))
    
    if 'openai' in Config().TYPE:
        claude_agent = Openai(pdf_path)
    elif 'claude' in Config().TYPE:
        claude_agent = Claude3(pdf_path)
    else:
        raise ValueError(f"Model {Config().MODEL} not supported.")
    
    claude_agent.initiate_propmt()

    # Filter and preprocess data
    params = Params()
    
    logging.info(colored('2. Filtering and preprocessing data', color='cyan', attrs=['bold']))
    filter_response = claude_agent.chat(Config().get_prompt('FILTER_PROMPT'))
    # time.sleep(5)
    logging.info(filter_response)
    params.parse_response(filter_response)
    adata = filter(adata, params)

    # Save filtered AnnData object
    adata.raw = adata

    logging.info(colored('3. Preprocessing data', color='cyan', attrs=['bold']))
    preprocess_response = claude_agent.chat(Config().get_prompt('PREPROCESSING_PROMPT'))
    # time.sleep(5)
    logging.info(preprocess_response)
    params.parse_response(preprocess_response)
    adata = preprocess(adata, params)
    
    # Clustering
    logging.info(colored('4. Clustering data', color='cyan', attrs=['bold']))
    clustering_response = claude_agent.chat(Config().get_prompt('CLUSTERING_PROMPT'))
    # time.sleep(5)
    logging.info(clustering_response)
    params.parse_response(clustering_response)
    adata = clustering(adata, params)
    
    # Annotate
    adata, marker_genes = get_marker_genes(adata, params)
    logging.info(colored('Top 10 marker genes for each cluster:', color='yellow'))
    logging.info(colored(marker_genes, color='yellow'))
    
    starting_part = 'This is the output of the top 10 marker genes for each cluster:'
    annotate_prompt = Config().get_prompt('ANNOTATION_PROMPT').replace(f"{starting_part}", 
                                                                       f"{starting_part}\n{marker_genes}")
    
    # time.sleep(5)
    logging.info(colored('5. Annotating clusters', color='cyan', attrs=['bold']))
    annotate_response = claude_agent.chat(annotate_prompt, max_tokens=1500)
    logging.info(colored(annotate_response, color='yellow'))
    annotation_dict = params.parse_annotation_response(annotate_response)
    if not params['reannotation']:
        adata = annotate(adata, annotation_dict, params, final=True)
    else:
        adata = annotate(adata, annotation_dict, params, final=False)
        logging.info(colored('6. Reannotating clusters', color='cyan', attrs=['bold']))
        query_response = claude_agent.chat(Config().get_prompt('REVIEW_PROMPT'))
        logging.info(query_response)
        params.parse_response(query_response)
        if len(params['genes_to_query']) > 0:
            # query gene expression data
            query_genes_exp_dict = query_datasets(adata, params)
            logging.info(colored('Gene expression data queried:' + str(query_genes_exp_dict), color='yellow'))
            if len(query_genes_exp_dict) > 0:
                middle_start = 'and the order is the same as the cluster number:'
                reannotate_prompt = Config().get_prompt('REANNOTATION_PROMPT').replace(f"{middle_start}", 
                                                                                       f"{middle_start}\n{query_genes_exp_dict}")
                reannotate_response = claude_agent.chat(reannotate_prompt, max_tokens=2000)
                logging.info(reannotate_response)
                reannotation_dict = params.parse_annotation_response(reannotate_response)
                adata = annotate(adata, reannotation_dict, params, final=True)
            else:
                logging.info(colored('No genes found in the dataset to query.', color='yellow'))
        else:
            logging.info(colored('No genes to query.', color='yellow'))
                
    # Save processed AnnData object
    logging.info(colored('Saving processed AnnData object', color='cyan', attrs=['bold']))
    logging.info(f'Saving processed AnnData object to {output_dir}/{output_name}')
    adata.write(os.path.join(output_dir, output_name))

    
