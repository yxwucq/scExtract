from .config import Config
from copy import deepcopy
from typing import List, Dict
import re

class Params:
    def __init__(self):
        self.default_params = Config().DEFAULT_PARAMS
        self.list_type_params = Config().LIST_TYPE_PARAMS
        self.params = deepcopy(self.default_params)
    
    @property
    def get_params(self):
        return self.params
    
    def reset_params(self):
        self.params = deepcopy(self.default_params)
    
    def parse_annotation_response(self, annotation_response: str) -> Dict[int, List[str]]:
        annotation_response = annotation_response.replace('\n', '')
        annotation_dict = re.search(r'annotation_dict: {.*}', annotation_response)
        annotation_dict = annotation_dict.group(0).replace('annotation_dict: ', '')
        
        try:
            annotation_dict = eval(annotation_dict)
        except: # sometimes the dictionary is not in the correct format
            annotation_dict = annotation_dict.replace("{", "").replace("}", "").replace("'", "")
            annotation_dict = annotation_dict + ','
            annotation_list = annotation_dict.split('],')
            annotation_dict = {}
            for annotation in annotation_list:
                if ':' not in annotation:
                    continue
                cluster = int(annotation.strip().split(':')[0].strip())
                cluster_list = annotation.split('[')[1].split(',')
                cluster_list = [item.strip() for item in cluster_list]
                annotation_dict[cluster] = cluster_list
        
        return annotation_dict
    
    def parse_response(self, filter_response: str):
        filter_response = filter_response.split('\n')
        
        for idx, line in enumerate(filter_response):
            lines = line.strip().split(':')
            if len(lines) != 2:
                continue
            elif lines[0].strip() in self.params:
                if lines[0].strip() not in self.list_type_params:
                # sometimes , is not trimmed in the value
                    self._update_params(lines[0].strip(), lines[1].split(',')[0].strip())
                else:
                    self._update_params(lines[0].strip(), lines[1].strip())
    
    def _update_params(self, key: str, value: str) -> None:
        if value in ['default', 'Default', 'DEFAULT']:
            return
        elif value in ['null', 'Null', 'NULL']:
            self.params[key] = None
            return
        else:
            # int parameters
            if key in ['filter_cells_low', 'filter_genes_low', 'filter_n_gene_by_counts_high', 'filter_total_counts_high',
                          'normalize_total_target_sum', 'highly_variable_genes_num', 'pca_comps', 'find_neighbors_neighbors_num',
                            'find_neighbors_using_pcs', 'leiden_or_louvain_group_numbers']:
                if not value.isdigit():
                    raise ValueError(f'Invalid value for {key}: {value}')
                self.params[key] = int(value)
                
            # float parameters
            elif key in ['filter_mito_percentage_low', 'filter_ribo_percentage_low']:
                if not value.replace('.', '', 1).isdigit():
                    raise ValueError(f'Invalid value for {key}: {value}')
                self.params[key] = float(value)
                
            # bool parameters
            elif key in ['batch_correction', 'log1p_transform', 'scale']:
                if value in ['True', 'true', 'TRUE']:
                    self.params[key] = True
                elif value in ['False', 'false', 'FALSE']:
                    self.params[key] = False
                else:
                    raise ValueError(f'Invalid value for {key}: {value}')
                
            # categorical parameters
            elif key in ['unsupervised_cluster_method', 'visualize_method']:
                if key == 'unsupervised_cluster_method':
                    if value in ['louvain', 'Louvain', 'LOUVAIN']:
                        self.params[key] = 'louvain'
                    elif value in ['leiden', 'Leiden', 'LEIDEN']:
                        self.params[key] = 'leiden'
                    else:
                        raise ValueError(f'Invalid Clustering method: {value}')
                elif key == 'visualize_method':
                    if value in ['umap', 'UMAP', 'Umap']:
                        self.params[key] = 'umap'
                    elif value in ['t-SNE', 'tsne', 'TSNE', 'T-SNE', 'tSNE']:
                        self.params[key] = 'tsne'
                    else:
                        raise ValueError(f'Invalid Visualization method: {value}')
            
            # list parameters
            elif key in self.list_type_params:
                if value== '[]':
                    self.params[key] = []
                else:
                    value = value.replace('[', '').replace(']', '').replace("'", "").split(',')
                    self.params[key] = [item.strip() for item in value]
            
            else:
                raise ValueError(f'Invalid key queried: {key}')
            
    def __str__(self):
        return str(self.params)
    
    def __repr__(self):
        return str(self.params)
    
    def __getitem__(self, key):
        return self.params[key]
    
                
            