from .utils.parse_args import parse_args

def main():
    args = parse_args()
    if args.subcommand is None:
        print("Please specify a subcommand. Use -h for help.")
    if args.subcommand == 'auto_extract':
        from .auto_extract.auto_extract import auto_extract
        auto_extract(adata_path=args.adata_path,
                        pdf_path=args.pdf_path,
                        output_dir=args.output_dir,
                        marker_genes_excel_path=args.marker_genes_excel_path,
                        config_path=args.config_path,
                        output_name=args.output_name,
                        output_log=args.output_log,
                        output_config_pkl=args.output_config_pkl,
                        benchmark_no_context_key=args.benchmark_no_context_key)
    elif args.subcommand == 'init':
        from .utils.init_project import init_project
        init_project(config_path=args.config_path,
                     overwrite=args.overwrite)
    elif args.subcommand == 'get_metadata':
        from .auto_extract.get_metadata import get_metadata_in_pdf
        get_metadata_in_pdf(pdf_list=args.pdf_list,
                            output_dir=args.output_dir,
                            initiation_samples=args.initiation_samples,
                            config_path=args.config_path,
                            output_name=args.output_name)
    elif args.subcommand == 'benchmark':
        from .benchmark.benchmark import benchmark_annotation
        benchmark_annotation(adata_path=args.adata_path,
                             true_group_key=args.true_group_key,
                             config_path=args.config_path,
                             predict_group_key=args.predict_group_key,
                             ontology=args.ontology,
                             method=args.method,
                             output_config_pkl=args.output_config_pkl,
                             similarity_key=args.similarity_key,
                             result_metrics_path=args.result_metrics_path,
                             output_path=args.output_path)
    elif args.subcommand == 'add_singler_annotation':
        from .methods_comparison.singler_anno import add_singler_annotation
        add_singler_annotation(adata_path=args.adata_path,
                               output_path=args.output_path,
                               ref_data=args.ref_data,
                               database_version=args.database_version,
                               key_added=args.key_added,
                               cache_dir=args.cache_dir)
    elif args.subcommand == 'add_celltypist_annotation':
        from .methods_comparison.celltypist_anno import add_celltypist_annotation
        add_celltypist_annotation(pdf_path=args.pdf_path,
                                  adata_path=args.adata_path,
                                  config_path=args.config_path,
                                  key_added=args.key_added,
                                  output_path=args.output_path)
    elif args.subcommand == 'add_sctype_annotation':
        from .methods_comparison.sctype_anno import add_sctype_annotation
        add_sctype_annotation(pdf_path=args.pdf_path,
                             adata_path=args.adata_path,
                             output_path=args.output_path,
                             key_added=args.key_added,
                             config_path=args.config_path)
    elif args.subcommand == 'extract_celltype_embedding':
        from .integration.extract_celltype_embedding import extract_celltype_embedding
        extract_celltype_embedding(file_list=args.file_list,
                                   config_path=args.config_path,
                                   output_embedding_pkl=args.output_embedding_pkl,
                                   cell_type_column=args.cell_type_column,
                                   output_individual_config_pkls=args.output_individual_config_pkls,
                                   )
    elif args.subcommand == 'integrate':
        from .integration.integrate import integrate_processed_datasets
        integrate_processed_datasets(file_list=args.file_list,
                                     method=args.method,
                                     output_path=args.output_path,
                                     config_path=args.config_path,
                                     alignment_path=args.alignment_path,
                                     embedding_dict_path=args.embedding_dict_path,
                                     downsample=args.downsample,
                                     downsample_cells_per_label=args.downsample_cells_per_label,
                                     search_factor=args.search_factor,
                                     approx=args.approx,
                                     use_gpu=args.use_gpu,
                                     batch_size=args.batch_size,
                                     use_pct=args.use_pct,
                                     dimred=args.dimred,
                                     )
    
if __name__ == '__main__':
    main()