from auto_extract.auto_extract import auto_extract
from benchmark.benchmark import benchmark_annotation
from methods_comparison.singler_anno import add_singler_annotation
from integration.integrate import integrate_processed_datasets
from utils.parse_args import parse_args

def main():
    args = parse_args()
    if args.subcommand is None:
        print("Please specify a subcommand. Use -h for help.")
    if args.subcommand == 'auto_extract':
        auto_extract(adata_path=args.adata_path,
                        pdf_path=args.pdf_path,
                        output_dir=args.output_dir,
                        output_name=args.output_name,
                        output_log=args.output_log,
                        output_config_pkl=args.output_config_pkl,
                        benchmark_no_context_key=args.benchmark_no_context_key)
    elif args.subcommand == 'benchmark':
        benchmark_annotation(adata_path=args.adata_path,
                             true_group_key=args.true_group_key,
                             predict_group_key=args.predict_group_key,
                             ontology=args.ontology,
                             method=args.method,
                             similarity_key=args.similarity_key,
                             result_metrics_path=args.result_metrics_path,
                             output_path=args.output_path)
    elif args.subcommand == 'add_singler_annotation':
        add_singler_annotation(adata_path=args.adata_path,
                               output_path=args.output_path,
                               ref_data=args.ref_data,
                               ref_features=args.ref_features,
                               ref_labels=args.ref_labels,
                               key_added=args.key_added,
                               cache_dir=args.cache_dir)
    elif args.subcommand == 'integrate':
        integrate_processed_datasets(file_list=args.file_list,
                                     method=args.method,
                                     output_path=args.output_path,
                                     prior_weight=args.prior_weight,
                                     prior_method=args.prior_method,
                                     alignment_path=args.alignment_path,
                                     embedding_dict_path=args.embedding_dict_path)
    
if __name__ == '__main__':
    main()