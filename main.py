from auto_extract.auto_extract import auto_extract
from utils.parse_args import parse_args

def main():
    args = parse_args()
    auto_extract(args.adata_path, args.pdf_path, args.output_dir, args.output_name)
    
if __name__ == '__main__':
    main()