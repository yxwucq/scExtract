import configparser
import os
from typing import Optional, List

import pandas as pd
from pyfiglet import Figlet
from tabulate import tabulate
from tqdm import tqdm

from .agent import Claude3, Openai
from .parse_params import Params

def wrap_string(long_string, max_line_length=20):
    words = long_string.split()
    lines = []
    current_line = ""

    for word in words:
        if len(current_line) + len(word) + 1 <= max_line_length:
            if current_line:
                current_line += " " + word
            else:
                current_line = word
        else:
            lines.append(current_line)
            current_line = word

    if current_line:
        lines.append(current_line)

    return "\n".join(lines)

def get_metadata_in_pdf(pdf_list: List[str],
                        output_dir: str = '.',
                        initiation_samples: bool = True,
                        config_path: str = 'config.ini',
                        output_name: str = 'metadata.csv',
                        ) -> None:
    """
    Extracts metadata from PDFs.
    """

    config = configparser.ConfigParser()
    config.read(config_path)

    sample_numbers = []
    article_titles = []
    authors = []
    journal_names = []
    sample_descriptions = []
    total_cells = []
    raw_data_sources = []
    abstracts = []

    for i, pdf_path in enumerate(pdf_list):
        print(f"Extracting {i}/{len(pdf_list)}...")
        print(f"Extracting metadata from {pdf_path}...")
        if 'openai' in config['API']['TYPE']:
            claude_agent = Openai(pdf_path, config_path)
        elif 'claude' in config['API']['TYPE']:
            claude_agent = Claude3(pdf_path, config_path)
        else:
            raise ValueError(f"Model {config['API']['MODEL']} not supported.")
        
        # Extract text from PDF
        params = Params(config_path)
        claude_agent.initiate_propmt()
        get_metadata_response = claude_agent.chat(params.get_prompt('GET_METADATA_PROMPT'))
        params.parse_response(get_metadata_response)
        summary_response = claude_agent.chat(params.get_prompt('SUMMARY_ARTICLE_PROMPT'))
        params.parse_response(summary_response)
        
        sample_numbers.append(i)
        article_titles.append(params['title'])
        authors.append(params['author'])
        journal_names.append(params['journal_name'])
        sample_descriptions.append(params['sample_description'])
        total_cells.append(params['total_cells'])
        raw_data_sources.append(params['raw_data_source'])
        abstracts.append(params['summary'])
        
        if initiation_samples:
            os.makedirs(os.path.join(output_dir, f"sample{i}/raw_data"), exist_ok=True)
            os.system(f'cp {pdf_path} {os.path.join(output_dir, f"sample{i}/raw_data", f"sample{i}.pdf")}')
        
    metadata = pd.DataFrame({'sample_number': sample_numbers,
                             'article_title': article_titles,
                             'author': authors,
                             'journal_name': journal_names,
                             'sample_description': sample_descriptions,
                             'total_cells': total_cells,
                             'raw_data_source': raw_data_sources,
                             'abstract': abstracts})
    
    print(f"Saving metadata to {os.path.join(output_dir, output_name)}...")
    metadata.to_csv(os.path.join(output_dir, output_name), index=False)
    metadata = metadata.map(lambda x: wrap_string(str(x), 20))
    metadata.drop('sample_number', axis=1, inplace=True)
    metadata.drop('abstract', axis=1, inplace=True)
    print(tabulate(metadata, headers='keys', tablefmt='fancy_grid'))