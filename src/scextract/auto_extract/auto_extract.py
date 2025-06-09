import re
import anndata as ad
import logging
import os
import warnings
import time

import pickle
from typing import Optional
from pyfiglet import Figlet
import colorama
from termcolor import colored

from .agent import Claude3, Openai, get_cell_type_embedding_by_llm
from .preprocess import filter, preprocess, clustering
from .annotation import get_marker_genes, annotate, query_datasets, simple_annotate
from .parse_params import Params 

import configparser
from ..utils.utils import convert_ensembl_to_symbol, get_top_markers

from .. import __version__

def auto_extract(adata_path: str,
                 pdf_path: str,
                 marker_genes_excel_path: str,
                 output_dir: str,
                 config_path: str = 'config.ini',
                 output_name: str = 'processed.h5ad',
                 output_log: str = 'auto_extract.log',
                 output_config_pkl: str = 'config.pkl',
                 benchmark_no_context_key: Optional[str] = None,
                 ) -> None:
    """
    Extracts and processes single-cell data from literature.
    """
    
    logging.basicConfig(level=logging.INFO)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file {config_path} not found. Please run 'init' to create a config file.")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.exists(os.path.join(output_dir, output_log)):
        os.remove(os.path.join(output_dir, output_log))

    config = configparser.ConfigParser()
    config.read(config_path)

    file_handler = logging.FileHandler(os.path.join(output_dir, output_log))
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)
    
    f = Figlet(font='slant')
    logging.info('\n'+f.renderText('scExtract')+f" v{__version__}"+'\n')
    colorama.init()
    logging.info(colored(f"Using {config['API']['MODEL']} as extraction model", color='cyan'))
    logging.info(colored(f"Using {config['API']['TOOL_MODEL']} as tool model", color='cyan'))
    
    # Load AnnData object
    logging.info(f'Loading AnnData object from {adata_path}')
    adata = ad.read_h5ad(adata_path)
    
    if adata.shape[0] > 50000:
        logging.warning(colored(f"Large dataset detected. The dataset contains {adata.shape[0]} cells. \
The result extracted by LLM may not be accurate. Please consider manually \
subsetting the data to a smaller size.", color='light_red'))
    
    # Check if Ensembl IDs are present in AnnData object
    for genes in adata.var.index[:10]:
        if genes.startswith('ENSG'):
            logging.info('Ensembl IDs detected in AnnData object.')
            adata = convert_ensembl_to_symbol(adata)
            adata.var_names_make_unique()
            break
    
    logging.info(colored('1. Extracting information from literature', color='cyan', attrs=['bold']))
    
    if 'openai' in config['API']['TYPE']:
        claude_agent = Openai(pdf_path, config_path)
    elif 'claude' in config['API']['TYPE']:
        claude_agent = Claude3(pdf_path, config_path)
    else:
        raise ValueError(f"Model {config['API']['MODEL']} not supported.")
    
    claude_agent.initiate_propmt()

    # Filter and preprocess data
    params = Params(config_path)
    
    logging.info(colored('2. Filtering and preprocessing data', color='cyan', attrs=['bold']))
    filter_response = claude_agent.chat(params.get_prompt('FILTER_PROMPT'))
    logging.info(filter_response)
    params.parse_response(filter_response)
    adata = filter(adata, params)

    logging.info(colored('3. Preprocessing data', color='cyan', attrs=['bold']))
    preprocess_response = claude_agent.chat(params.get_prompt('PREPROCESSING_PROMPT'))
    logging.info(preprocess_response)
    params.parse_response(preprocess_response)
    adata = preprocess(adata, params)
    
    # Clustering
    logging.info(colored('4. Clustering data', color='cyan', attrs=['bold']))
    clustering_response = claude_agent.chat(params.get_prompt('CLUSTERING_PROMPT'))
    logging.info(clustering_response)
    params.parse_response(clustering_response)
    adata = clustering(adata, params)
    
    if config['OPTIONS'].getboolean('CLEAN_INTERMEDIATE_MESSAGES'):
        logging.info(colored('Cleaning up intermediate messages', color='cyan', attrs=['bold']))
        claude_agent.clear_intermediate_messages()
    
    # Annotate
    logging.info(colored('5. Getting marker genes', color='cyan', attrs=['bold']))
    adata, marker_genes = get_marker_genes(adata, params, fast_mode=config['OPTIONS'].getboolean('FAST_MODE'))
    logging.info(colored('Top 10 marker genes for each cluster:', color='yellow'))
    logging.info(colored(marker_genes, color='yellow'))
    
    if config['OPTIONS'].getboolean('BENCHMARK_GPTCELLTYPE'):
        logging.info(colored('Benchmarking GPTCellType', color='cyan', attrs=['bold']))
        tissue_name = claude_agent.chat(params.get_prompt('GET_TISSUE_NAME_PROMPT'))
        logging.info(colored(f'Tissue name: {tissue_name}', color='yellow'))
        benchmark_gptcelltype_prompt = params.get_tool_prompt('GPTCELLTYPE_ANNOTATION_PROMPT')
        benchmark_gptcelltype_prompt += "\n".join([f"{k}:{','.join(v)}" for k, v in marker_genes.items()])
        benchmark_gptcelltype_response = claude_agent._tool_retrieve(messages=[{"role": "user", "content": benchmark_gptcelltype_prompt.replace('{tissuename}', tissue_name)}])
        logging.info(colored(benchmark_gptcelltype_response, color='green'))
        benchmark_gptcelltype_response_list = [x for x in benchmark_gptcelltype_response.split('\n') if x and not x.startswith('Here are')]
        if any(re.match(r'^\d+:\s*', x) for x in benchmark_gptcelltype_response_list):
            benchmark_gptcelltype_response_list = [x for x in benchmark_gptcelltype_response_list if re.match(r'^\d+:\s*', x)]
            benchmark_gptcelltype_response_list = [re.sub(r'^\d+:\s*', '', string) for string in benchmark_gptcelltype_response_list]
        assert len(benchmark_gptcelltype_response_list) == len(marker_genes), 'Number of responses does not match number of clusters.'
        gptcelltype_annotation_dict = {k: v for k, v in enumerate(benchmark_gptcelltype_response_list)}
        adata = simple_annotate(adata, gptcelltype_annotation_dict, params, 'gptcelltype_annotation')

    if benchmark_no_context_key is not None:
        benchmark_no_context_prompt = params.get_tool_prompt('NO_CONTEXT_ANNOTATION_PROMPT').replace('Some can be a mixture of multiple cell types.',
                                                                                        'Some can be a mixture of multiple cell types.' + str(marker_genes))
        benchmark_no_context_summary = claude_agent._tool_retrieve(messages=[{"role": "user", "content": benchmark_no_context_prompt}])
        logging.info(colored(benchmark_no_context_summary, color='green'))
        no_context_annotation_dict = params.parse_annotation_response(benchmark_no_context_summary, simple_annotation=True)
        adata = simple_annotate(adata, no_context_annotation_dict, params, benchmark_no_context_key)
    
    if marker_genes_excel_path is not None:
        logging.info(colored('Loading marker genes from file', color='cyan', attrs=['bold']))
        top_k_markers = get_top_markers(marker_genes_excel_path, config_path=config_path)
        logging.info(colored('Top marker genes for each cluster:', color='yellow'))
        logging.info(colored('\n'.join(top_k_markers), color='yellow'))
    
    starting_part = 'This is the output of the top 10 marker genes for each cluster:'
    annotate_prompt = params.get_prompt('ANNOTATION_PROMPT').replace(f"{starting_part}", 
                                                                       f"{starting_part}\n{marker_genes}")
    
    verbose_annot = config['OPTIONS'].getboolean('ADD_VERBOSE_ANNOTATIONS')
    if verbose_annot:
        annotate_prompt = params.get_prompt('ANNOTATION_PROMPT_VERBOSE').replace(f"{starting_part}", 
                                                                       f"{starting_part}\n{marker_genes}")
        
    if marker_genes_excel_path is not None:    
        annotate_prompt = annotate_prompt.replace('{authors_defined_marker_genes}', "\nThe marker genes for each celltype provided by the author's supplementary data are:\n{}\n Remember to assign only standard *cell ontology types* to the clusters\n".format('\n'.join(top_k_markers)))
    else:
        annotate_prompt = annotate_prompt.replace('{authors_defined_marker_genes}', '')
    
    logging.info(colored('6. Annotating clusters', color='cyan', attrs=['bold']))
    annotate_response = claude_agent.chat(annotate_prompt, max_tokens=1500)
    logging.info(colored(annotate_response, color='yellow'))
    annotation_dict = params.parse_annotation_response(annotate_response)
    if not params['reannotation']:
        adata = annotate(adata, annotation_dict, params, final=True, verbose_annot=verbose_annot)
    else:
        adata = annotate(adata, annotation_dict, params, final=False, verbose_annot=verbose_annot)
        logging.info(colored('7. Reannotating clusters', color='cyan', attrs=['bold']))
        query_response = claude_agent.chat(params.get_prompt('REVIEW_PROMPT'))
        logging.info(query_response)
        params.parse_response(query_response)
        if len(params['genes_to_query']) > 0:
            # query gene expression data
            query_genes_exp_dict = query_datasets(adata, params)
            logging.info('Gene expression data queried:' + str(query_genes_exp_dict))
            if len(query_genes_exp_dict) > 0:
                smmary_message = params.get_tool_prompt('SUMMARY_QUERY_EXPRESSION')
                smmary_message += str(query_genes_exp_dict)
                query_genes_exp_dict_summary = claude_agent._tool_retrieve(messages=[{"role": "user", "content": smmary_message}])
                logging.info(colored(query_genes_exp_dict_summary, color='yellow'))
                middle_start = 're-annotate the clusters into cell types using a dictionary format:'
                reannotate_prompt = params.get_prompt('REANNOTATION_PROMPT').replace(f"{middle_start}", 
                                                                                       f"{middle_start}\n{query_genes_exp_dict_summary}")
                if verbose_annot:
                    reannotate_prompt = params.get_prompt('REANNOTATION_PROMPT_VERBOSE').replace(f"{middle_start}", 
                                                                                       f"{middle_start}\n{query_genes_exp_dict_summary}")
                
                reannotate_response = claude_agent.chat(reannotate_prompt, max_tokens=2000)
                logging.info(reannotate_response)
                reannotation_dict = params.parse_annotation_response(reannotate_response)
                adata = annotate(adata, reannotation_dict, params, final=True, verbose_annot=verbose_annot)
            else:
                logging.info(colored('No genes found in the dataset to query.', color='yellow'))
                adata = annotate(adata, annotation_dict, params, final=True, verbose_annot=verbose_annot)
        else:
            logging.info(colored('No genes to query.', color='yellow'))
            adata = annotate(adata, annotation_dict, params, final=True, verbose_annot=verbose_annot)

    if config['OPTIONS'].getboolean('ADD_CELLTYPE_DESCRIPTION'):
        logging.info(colored('Adding cell type description', color='cyan', attrs=['bold']))
        cluster_key = params.get_params['unsupervised_cluster_method']
        celltype_list = list(adata.obs[cluster_key].unique())
        celltype_description_prompt = params.get_prompt('ADD_CELLTYPE_DESCRIPTION_PROMPT').replace('{celltype_list}', str(celltype_list))
        # Description with no context
        celltype_description_response = claude_agent._tool_retrieve(messages=[{"role": "user", "content": celltype_description_prompt}], max_tokens=2000)
        logging.info(celltype_description_response)
        celltype_description_dict = params.parse_annotation_response(celltype_description_response, simple_annotation=True, annotation_type='celltype_descriptors')
        adata = simple_annotate(adata, celltype_description_dict, params, f'{cluster_key}_description', cluster_key)

    if config['OPTIONS'].getboolean('ADD_CELLTYPE_FUNCTION'):
        logging.info(colored('Adding cell type function', color='cyan', attrs=['bold']))
        cluster_key = params.get_params['unsupervised_cluster_method']
        celltype_list = list(adata.obs[cluster_key].unique())
        celltype_function_prompt = params.get_prompt('ADD_CELLTYPE_FUNCTION_PROMPT').replace('{celltype_list}', str(celltype_list))
        celltype_function_response = claude_agent.chat(celltype_function_prompt, max_tokens=2000)
        logging.info(celltype_function_response)
        celltype_function_dict = params.parse_annotation_response(celltype_function_response, simple_annotation=True, annotation_type='celltype_functions')
        adata = simple_annotate(adata, celltype_function_dict, params, f'{cluster_key}_function', cluster_key)

    if config['API'].getboolean('CONVERT_EMBEDDING'):
        params.embedding_dict = {}
        for key in adata.obs.columns:
            if key in ['leiden', 'louvain', 'tissue', 'disease', 'developmental_stage']:
                words_list = list(adata.obs[key].unique())
                embedding_list = get_cell_type_embedding_by_llm(words_list, config_path=config_path)
                for i in range(len(words_list)):
                    params.embedding_dict[words_list[i]] = embedding_list[i]
    
    # Save processed AnnData object
    logging.info(colored('Saving processed AnnData object and config file', color='cyan', attrs=['bold']))
    logging.info(f'Saving processed config file to {os.path.join(output_dir,output_config_pkl)}')
    with open(os.path.join(output_dir, output_config_pkl), 'wb') as f:
        pickle.dump(params, f)
    logging.info(f'Saving processed AnnData object to {os.path.join(output_dir,output_name)}')
    adata.write(os.path.join(output_dir, output_name))

    
