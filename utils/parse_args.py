import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Automatically extract, process, and annotate single-cell data from literature.')
    subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
    
    auto_extract_parser = subparsers.add_parser('auto_extract', help='Automatically extract, process, and annotate single-cell data from literature.')
    
    auto_extract_parser.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the raw data in AnnData format.')
    auto_extract_parser.add_argument('--pdf_path', '-p', type=str, required=True, help='Path to the PDF file containing the article. should in pdf or txt format.')
    auto_extract_parser.add_argument('--output_dir', '-d', type=str, default='ExtractedData', help='Directory to save the processed data.')
    auto_extract_parser.add_argument('--output_name', '-o', type=str, default='processed.h5ad', help='Name of the output file.')
    auto_extract_parser.add_argument('--benchmark_no_context_key', '-b', type=str, default=None, help='If specified, Directly get annotation from marker genes without article context for benchmarking, \
                                    the result will be saved in adata.obs[benchmark_no_context_key].')
    
    benchmark_parser = subparsers.add_parser('benchmark', help='Benchmark annotation results using true labels.')
    
    benchmark_parser.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    benchmark_parser.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    benchmark_parser.add_argument('--true_group_key', '-t', type=str, required=True, help='Key of the true group in adata.obs.')
    benchmark_parser.add_argument('--predict_group_key', '-p', type=str, help='Key of the predicted group in adata.obs. Support multiple keys separated by comma.')
    benchmark_parser.add_argument('--ontology', '-l', type=str, default='cl', help='Ontology to use for annotation.')
    benchmark_parser.add_argument('--method', '-m', type=str, default='ols_api', help='Method to use for annotation.')
    benchmark_parser.add_argument('--similarity_key', '-s', type=str, default='similarity', help='Key to save the similarity results. Support multiple keys separated by comma. Order should be the same as predict_group_key.')

    add_singler_annotation = subparsers.add_parser('add_singler_annotation', help='Annotate single-cell data using py&c++ implementation of singler.')
    
    add_singler_annotation.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    add_singler_annotation.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    add_singler_annotation.add_argument('--ref_data', '-d', type=str, default='HumanPrimaryCellAtlasData', help='Reference data to use for annotation.')
    add_singler_annotation.add_argument('--ref_features', '-f', type=str, default='symbol', help='Reference features to use for annotation.')
    add_singler_annotation.add_argument('--ref_labels', '-l', type=str, default='main', help='Reference labels to use for annotation.')
    add_singler_annotation.add_argument('--cache_dir', '-c', type=str, default='_cache', help='Directory to save the cache files.')

    return parser.parse_args()