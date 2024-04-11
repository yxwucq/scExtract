import anndata as ad
import logging
from pyfiglet import Figlet
import os
import warnings
import time

from .claude3 import Claude3
from .preprocess import filter, preprocess, clustering
from .annotation import get_marker_genes, annotate
from .parse_params import Params 
from .config import Config

def auto_extract(adata_path: str, 
                 pdf_path: str, 
                 output_dir: str,
                 output_name: str = 'processed.h5ad') -> None:
    """
    Extracts and processes single-cell data from literature.
    
    Args:
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

    # Load AnnData object
    logging.info(f'Loading AnnData object from {adata_path}')
    adata = ad.read_h5ad(adata_path)
    
    logging.info('Extracting information from literature')
    claude_agent = Claude3(pdf_path)
    claude_agent.initiate_propmt()
    # time.sleep(5)
    
    # Filter and preprocess data
    params = Params()
    
    filter_response = claude_agent.chat(Config().get_prompt('FILTER_PROMPT'))
    # time.sleep(5)
    logging.info(filter_response)
    params.parse_response(filter_response)
    adata = filter(adata, params)

    # Save filtered AnnData object
    adata.raw = adata

    preprocess_response = claude_agent.chat(Config().get_prompt('PREPROCESSING_PROMPT'))
    # time.sleep(5)
    logging.info(preprocess_response)
    params.parse_response(preprocess_response)
    adata = preprocess(adata, params)
    
    # Clustering
    clustering_response = claude_agent.chat(Config().get_prompt('CLUSTERING_PROMPT'))
    # time.sleep(5)
    logging.info(clustering_response)
    params.parse_response(clustering_response)
    adata = clustering(adata, params)
    
    # Annotate
    adata, marker_genes = get_marker_genes(adata, params)
    
    starting_part = 'This is the output of the top 10 marker genes for each cluster:'
    annotate_prompt = Config().get_prompt('ANNOTATION_PROMPT').replace(f"{starting_part}", 
                                                                       f"{starting_part}\n{marker_genes}")
    
    # time.sleep(5)
    annotate_response = claude_agent.chat(annotate_prompt)
    logging.info(annotate_response)
    annotation_dict = params.parse_annotation_response(annotate_response)
    adata = annotate(adata, annotation_dict, params)
    
    # Save processed AnnData object
    logging.info(f'Saving processed AnnData object to {output_dir}/{output_name}')
    adata.write(os.path.join(output_dir, output_name))

    
