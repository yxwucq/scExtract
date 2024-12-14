import celltypist
import configparser
from celltypist import models
import anndata as ad
import re
import pandas as pd
from pathlib import Path
import tempfile
from ..auto_extract.agent import Claude3, Openai, get_cell_type_embedding_by_llm
from ..auto_extract.parse_params import Params 

def format_celltypist_model_list(df: pd.DataFrame) -> str:
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Invalid model list.")
    
    formatted_pairs = []
    for i in range(len(df)):
        formatted_pairs.append(f"{df.iloc[i]['model']}:{df.iloc[i]['description']}")
    
    return ', '.join(formatted_pairs)

def get_best_model(pdf_path: str, config_path: str = 'config.ini') -> str:
    config = configparser.ConfigParser()
    config.read(config_path)

    if 'openai' in config['API']['TYPE']:
        claude_agent = Openai(pdf_path, config_path)
    elif 'claude' in config['API']['TYPE']:
        claude_agent = Claude3(pdf_path, config_path)
    else:
        raise ValueError(f"Model {config['API']['MODEL']} not supported.")
    
    df_models = models.models_description()
    model_list_str = format_celltypist_model_list(df_models)
    
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

def add_celltypist_annotation(pdf_path: str,
                              adata_path: str,
                              config_path: str = 'config.ini',
                              key_added: str = 'celltypist_annotation',
                              output_path: str = None,
                                ):
    
    model_name = get_best_model(pdf_path, config_path)
    models.download_models(model = model_name)
    model = models.Model.load(model = model_name)
    
    adata = ad.read_h5ad(adata_path)
    # save raw counts as temoprary file
    temp_dir_obj = tempfile.TemporaryDirectory()
    tempdir = temp_dir_obj.name
    print(f"Created temporary directory at: {tempdir}")
    
    raw_count_df = pd.DataFrame(adata.layers['counts'].toarray(),
                                columns=adata.var.index,
                                index=adata.obs.index)
    
    input_file = f"{tempdir}/raw_counts.csv"
    raw_count_df.to_csv(input_file)
    
    predictions = celltypist.annotate(input_file, model = model)
    adata.obs[key_added] = predictions.predicted_labels['predicted_labels']
    
    temp_dir_obj.cleanup()
    
    # save adata
    if output_path is not None:
        adata.write_h5ad(output_path)
    else:
        adata.write_h5ad(adata_path)