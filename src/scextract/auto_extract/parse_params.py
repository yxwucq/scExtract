from copy import deepcopy
from typing import List, Dict
import functools
import re
import configparser
import logging
from termcolor import colored

from .agent import Claude3, Openai
from ..utils.prompts import Prompts

class Params:
    def __init__(self, config_path: str = 'config.ini'):
        config = configparser.ConfigParser()
        config.read(config_path)
        self.config = config
        self.config_path = config_path
        self.default_params = deepcopy(config['DEFAULT_PARAMS'])
        self.params = format_params(self.default_params)
        self.list_type_params = Prompts().LIST_TYPE_PARAMS
        self.int_type_params = Prompts().INT_TYPE_PARAMS
        self.bool_type_params = Prompts().BOOL_TYPE_PARAMS
        self.categorical_params = Prompts().CATEGORICAL_PARAMS
    
    @property
    def get_params(self):
        return self.params
    
    def reset_params(self):
        self.params = dict(self.default_params)
        
    def parse_retry_decorator(func):
        """
        Decorator to retry the function if it fails
        """
        @functools.wraps(func)
        def wrapper(self, response, *args, **kwargs):
            try:
                return func(self, response, *args, **kwargs)
            except Exception as e:
                logging.warning(colored(f"Response parsing failed: {e}, trying to reformat the response", 'magenta'))
                logging.warning(colored(f"Original response: {response}", 'magenta'))
                if 'openai' in self.config['API']['TYPE']:
                    claude_agent = Openai(pdf_path=None, config_path=self.config_path)
                elif 'claude' in self.config['API']['TYPE']:
                    claude_agent = Claude3(pdf_path=None, config_path=self.config_path)
                
                reformat_prompts = Prompts().get_tool_prompt('REFORMAT_RESPONSE_PROMPT') + "\n" + response + "\n" + "The Error is: " + str(e)
                new_response = claude_agent._tool_retrieve(messages=[{"role": "user", "content": reformat_prompts}])
                logging.warning(colored(f"Reformatted response:\n{new_response}", 'magenta'))
                if self.config['OPTIONS']['MANUAL_REFORMAT_RESPONSE']:
                    try:
                        return func(self, new_response, *args, **kwargs)
                    except Exception as e:
                        logging.warning(f"Automatic reformatting failed: {str(e)}")
                        logging.warning(colored(f"Original response: {response}", 'magenta'))
                        logging.warning("Please reformat the response manually")
                        logging.warning("Input 'END' on a new line to end the input")
                        
                        lines = []
                        while True:
                            line = input()
                            if line.strip() == 'END':
                                break
                            lines.append(line)
                        
                        logging.warning(colored(f"Submitted new response", 'magenta'))
                        new_response = '\n'.join(lines)
                        return func(self, new_response, *args, **kwargs)
                return func(self, new_response, *args, **kwargs)
        return wrapper
    
    @parse_retry_decorator
    def parse_annotation_response(self,
                                  annotation_response: str,
                                  simple_annotation: bool = False,
                                  annotation_type: str = 'annotation_dict',
                                  ) -> Dict[int, List[str]|str]:
        annotation_response = '\n'.join([x.split('#')[0] for x in annotation_response.split('\n')])
        annotation_response = annotation_response.replace('\n', '')
        try:
            annotation_dict = re.search(annotation_type+r'\s*[:=]\s*({.*?})', annotation_response, re.DOTALL)
        except:
            annotation_dict = re.search(r'({.*?})', annotation_response, re.DOTALL)
        annotation_dict = annotation_dict.group(1)
        
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
    
    @parse_retry_decorator
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
                if '.' in value:
                    self.params[key] = float(value)
                else:
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
                        Warning(f"Invalid Clustering method: {value}, using default: {self.params[key]}")
                elif key == 'visualize_method':
                    if value in ['umap', 'UMAP', 'Umap']:
                        self.params[key] = 'umap'
                    elif value in ['t-SNE', 'tsne', 'TSNE', 'T-SNE', 'tSNE']:
                        self.params[key] = 'tsne'
                    else:
                        Warning(f"Invalid Visualization method: {value}, using default: {self.params[key]}")
            
            # list parameters
            elif key in self.list_type_params:
                if value== '[]':
                    pass
                else:
                    value = value.replace('[', '').replace(']', '').replace("'", "").split(',')
                    self.params[key] = self.params[key] + [item.strip() for item in value]

            # string parameters
            else:
                self.params[key] = value
            
    def __str__(self):
        return str(self.params)
    
    def __repr__(self):
        return str(self.params)
    
    def __getitem__(self, key):
        return self.params[key]
    
    def get_prompt(self, prompt_name: str) -> str:
        return Prompts().get_prompt(prompt_name)
    
    def get_tool_prompt(self, prompt_name: str) -> str:
        return Prompts().get_tool_prompt(prompt_name)

def format_params(default_params) -> dict:
    formatted_params = {}
    for key, value in default_params.items():
        if value == 'None':
            formatted_params[key] = None
        elif value.lower() == 'true':
            formatted_params[key] = True
        elif value.lower() == 'false':
            formatted_params[key] = False
        elif value == '':
            formatted_params[key] = []
        elif ',' in value:
            formatted_params[key] = value.split(',')
        else:
            value = value.strip()
            try:
                formatted_params[key] = int(value)
            except ValueError:
                try:
                    formatted_params[key] = float(value)
                except ValueError:
                    formatted_params[key] = value
                
    formatted_params['title'] = 'N/A'
    formatted_params['author'] = 'N/A'
    formatted_params['journal_name'] = 'N/A'
    formatted_params['sample_description'] = 'N/A'
    formatted_params['total_cells'] = 'N/A'
    formatted_params['raw_data_source'] = 'N/A'
    formatted_params['summary'] = 'N/A'
    
    return formatted_params