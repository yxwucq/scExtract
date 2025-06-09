import argparse

try:
    from scextract import __version__
except ImportError:
    __version__ = "unknown"

def parse_args():
    parser = argparse.ArgumentParser(
        description='Automatically extract, process, and integrate single-cell data from literature.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter    
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
    
    auto_extract_parser = subparsers.add_parser('auto_extract', help='Automatically extract, process, and annotate single-cell data from literature.',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    auto_extract_parser.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the raw data in AnnData format.')
    auto_extract_parser.add_argument('--pdf_path', '-p', type=str, required=True, help='Path to the PDF file containing the article. should in pdf or txt format.')
    auto_extract_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    auto_extract_parser.add_argument('--output_dir', '-d', type=str, default='Processed', help='Directory to save the processed data.')
    auto_extract_parser.add_argument('--output_name', '-o', type=str, default='processed.h5ad', help='Name of the output file.')
    auto_extract_parser.add_argument('--output_config_pkl', '-c', type=str, default='config.pkl', help='Name of the output config file, storing the config of the extraction and embeddings')
    auto_extract_parser.add_argument('--output_log', '-l', type=str, default='auto_extract.log', help='Name of the output log file.')
    auto_extract_parser.add_argument('--marker_genes_excel_path', '-m', type=str, help='Path to the Excel file containing the marker genes. should in xlsx or xls format.')
    auto_extract_parser.add_argument('--benchmark_no_context_key', '-b', type=str, default=None, help='If specified, Directly get annotation from marker genes without article context for benchmarking, \
                                    the result will be saved in adata.obs[benchmark_no_context_key].')
    
    init_parser = subparsers.add_parser('init', help='Initialize a project by creating a config file.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    init_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='Path to save the config file.')
    init_parser.add_argument('--overwrite', '-o', action='store_true', help='Overwrite the existing config file.')
    
    get_metadata_parser = subparsers.add_parser('get_metadata', help='Extract metadata from PDFs.',
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    get_metadata_parser.add_argument('--pdf_list', '-i', type=str, nargs='+', required=True, help='List of paths to the PDF files.')
    get_metadata_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    get_metadata_parser.add_argument('--output_dir', '-d', type=str, default='.', help='Directory to save the metadata.')
    get_metadata_parser.add_argument('--output_name', '-o', type=str, default='metadata.csv', help='Name of the output metadata file.')
    get_metadata_parser.add_argument('--initiation_samples', '-s', action='store_true', help='Whether to initiate projects')
    
    benchmark_parser = subparsers.add_parser('benchmark', help='Benchmark annotation results using true labels.',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    benchmark_parser.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    benchmark_parser.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    benchmark_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    benchmark_parser.add_argument('--output_config_pkl', '-c', type=str, default='config.pkl', help='Name of the output config file, storing the extracted processing parameters and cell type embeddings')
    benchmark_parser.add_argument('--result_metrics_path', '-r', type=str, help='Path to save the metrics of the benchmark results.')
    benchmark_parser.add_argument('--method', '-m', type=str, default='ols_api', help='Method to use for benchmarking annotation. Support ols_api and embedding.', choices=['ols_api', 'embedding'])
    benchmark_parser.add_argument('--predict_group_key', '-p', type=str, help='Key of the predicted group in adata.obs. Support multiple keys separated by comma.')
    benchmark_parser.add_argument('--true_group_key', '-t', type=str, required=True, help='Key of the true group in adata.obs.')
    benchmark_parser.add_argument('--similarity_key', '-s', type=str, default='similarity', help='Key to save the similarity results. Support multiple keys separated by comma. Order should be the same as predict_group_key.')
    benchmark_parser.add_argument('--ontology', '-l', type=str, default='cl', help='Ontology to use for annotation. Only valid for ols_api method.')

    add_singler_annotation = subparsers.add_parser('add_singler_annotation', help='Annotate single-cell data using py&c++ implementation of singler.',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    add_singler_annotation.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    add_singler_annotation.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    add_singler_annotation.add_argument('--ref_data', '-d', type=str, default='hpca', help='Reference data to use for annotation. See celldex for available references.')
    add_singler_annotation.add_argument('--database_version', '-v', type=str, default='2024-02-26', help='Version of the reference database to use for annotation.')
    add_singler_annotation.add_argument('--cache_dir', '-c', type=str, default='_cache', help='Directory to save the cache files.')
    add_singler_annotation.add_argument('--key_added', '-k', type=str, default='singler_annotation', help='Key to save the annotation results in adata.obs.')

    add_sctype_annotation = subparsers.add_parser('add_sctype_annotation', help='Annotate single-cell data using py implementation of scType.',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    add_sctype_annotation.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    add_sctype_annotation.add_argument('--pdf_path', '-p', type=str, required=True, help='Path to the PDF file containing the article. should in pdf or txt format.')
    add_sctype_annotation.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    add_sctype_annotation.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    add_sctype_annotation.add_argument('--key_added', '-k', type=str, default='sctype_annotation', help='Key to save the annotation results in adata.obs.')

    add_celltypist_annotation = subparsers.add_parser('add_celltypist_annotation', help='Annotate single-cell data using celltypist.',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    add_celltypist_annotation.add_argument('--pdf_path', '-p', type=str, required=True, help='Path to the PDF file containing the article. should in pdf or txt format.')
    add_celltypist_annotation.add_argument('--adata_path', '-i', type=str, required=True, help='Path to the processed data in AnnData format.')
    add_celltypist_annotation.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    add_celltypist_annotation.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    add_celltypist_annotation.add_argument('--key_added', '-k', type=str, default='celltypist_annotation', help='Key to save the annotation results in adata.obs.')

    extract_celltype_embedding_parser = subparsers.add_parser('extract_celltype_embedding', help='Extract cell type embeddings from the processed datasets.',
                                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    extract_celltype_embedding_parser.add_argument('--file_list', '-i', type=str, nargs='+', required=True, help='List of paths to the processed data in AnnData format.')
    extract_celltype_embedding_parser.add_argument('--output_embedding_pkl', '-o', type=str, required=True, help='Path to save the output embedding dictionary.')
    extract_celltype_embedding_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    extract_celltype_embedding_parser.add_argument('--cell_type_column', '-l', type=str, help='Column name of the cell type in the processed data.')
    extract_celltype_embedding_parser.add_argument('--output_individual_config_pkls', '-c', type=str, help='Name of the output config files of scExtract auot_extract, storing the config of the extraction and embeddings. \
                                            If specified, the cell type embeddings will be extracted from the specified config file.')
    integrate_parser = subparsers.add_parser('integrate', help='Integrate multiple processed datasets.',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    integrate_parser.add_argument('--file_list', '-i', type=str, nargs='+', required=True, help='List of paths to the processed data in AnnData format.')
    integrate_parser.add_argument('--config_path', '-f', type=str, default='config.ini', help='System config file path.')
    integrate_parser.add_argument('--method', '-m', type=str, default='cellhint_prior', choices=['scExtract', 'scanorama_prior', 'scanorama', 'cellhint_prior', 'cellhint'], help='Method to use for integration. Support scExtract, scanorama_prior, cellhint_prior and cellhint.')
    integrate_parser.add_argument('--output_path', '-o', type=str, help='Path to save the output file. If not specified, the input file will be overwritten.')
    integrate_parser.add_argument('--alignment_path', '-a', type=str, help='Path to the output alignment file.')
    integrate_parser.add_argument('--embedding_dict_path', '-e', type=str, help='Path to the cell type embedding dictionary. Only used for local. Can be generated by extract_celltype_embedding.')
    integrate_parser.add_argument('--downsample', '-d', action='store_true', help='Whether to downsample the cells.')
    integrate_parser.add_argument('--downsample_cells_per_label', '-c', type=int, default=1000, help='Number of cells to downsample per label.')
    integrate_parser.add_argument('--search_factor', '-s', type=int, default=5, help='Increased fold of the search space for the prior optimization. Only valid for scanorama_prior when approx is True.')
    integrate_parser.add_argument('--approx', '-x', action='store_true', help='Whether to use approximate optimization.')
    integrate_parser.add_argument('--use_gpu', '-g', action='store_true', help='Whether to use GPU for optimization. Cupy is required.')
    integrate_parser.add_argument('--batch_size', '-b', type=int, default=5000, help='Batch size for processing. Lower batch size will use less memory.')
    integrate_parser.add_argument('--use_pct', '-p', action='store_true', help='Whether to use pct for cellhint. Default is False to align with original cellhint.')
    integrate_parser.add_argument('--dimred', '-r', type=int, default=100, help='Dimensionality of the integrated embedding.')
    return parser.parse_args()