import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Automatically extract, process, and annotate single-cell data from literature.')
    
    parser.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the raw data in AnnData format.')
    parser.add_argument('--pdf_path', '-p', type=str, required=True, help='Path to the PDF file containing the article. should in pdf or txt format.')
    parser.add_argument('--output_dir', '-d', type=str, default='ExtractedData', help='Directory to save the processed data.')
    parser.add_argument('--output_name', '-o', type=str, default='processed.h5ad', help='Name of the output file.')

    return parser.parse_args()