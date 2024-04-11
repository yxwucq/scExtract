from auto_extract.auto_extract import auto_extract

def main():
    adata_path = '/home/wu/datb1/AutoExtractSingleCell/raw_data/raw_data.h5ad'
    pdf_path = '/home/wu/datb1/AutoExtractSingleCell/extracted_text.txt'
    output_path = '/home/wu/datb1/AutoExtractSingleCell/processed_data'
    
    auto_extract(adata_path, pdf_path, output_path)
    
if __name__ == '__main__':
    main()