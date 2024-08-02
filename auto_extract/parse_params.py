from .config import Config
from copy import deepcopy
from typing import List, Dict
import re

class Params:
    def __init__(self):
        self.default_params = Config().DEFAULT_PARAMS
        self.list_type_params = Config().LIST_TYPE_PARAMS
        self.int_type_params = Config().INT_TYPE_PARAMS
        self.bool_type_params = Config().BOOL_TYPE_PARAMS
        self.categorical_params = Config().CATEGORICAL_PARAMS
        self.params = deepcopy(self.default_params)
    
    @property
    def get_params(self):
        return self.params
    
    def reset_params(self):
        self.params = deepcopy(self.default_params)
    
    def parse_annotation_response(self, 
                                  annotation_response: str,
                                  simple_annotation: bool = False,
                                  ) -> Dict[int, List[str]|str]:
        annotation_response = annotation_response.replace('\n', '')
        annotation_dict = re.search(r'annotation_dict: {.*?}', annotation_response)
        annotation_dict = annotation_dict.group(0).replace('annotation_dict: ', '')
        
        try:
            annotation_dict = eval(annotation_dict)
        except: # sometimes the dictionary is not in the correct format
            if not simple_annotation:
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
            else:
                annotation_dict = annotation_dict.replace("{", "").replace("}", "").replace("'", "")
                annotation_dict = annotation_dict.split(',')
                annotation_dict = {int(item.split(':')[0].strip()): item.split(':')[1].strip() for item in annotation_dict}
        
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
            if key in self.int_type_params:
                if 'percentage' in key:
                    if value.startswith('0.'):
                        value = str(int(float(value) * 100))
                if not value.isdigit():
                    raise ValueError(f'Invalid value for {key}: {value}')
                self.params[key] = int(value)
                
            # bool parameters
            elif key in self.bool_type_params:
                if value in ['True', 'true', 'TRUE']:
                    self.params[key] = True
                elif value in ['False', 'false', 'FALSE']:
                    self.params[key] = False
                else:
                    raise ValueError(f'Invalid value for {key}: {value}')
                
            # categorical parameters
            elif key in self.categorical_params:
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
    
                
            