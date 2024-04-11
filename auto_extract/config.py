class Config:
    API_KEY = ''
    API_BASE_URL = "https://api.claude-plus.top/v1"
    
    DEFAULT_PARAMS = {
        'filter_cells_low': 300,
        'filter_genes_low': 3,
        'filter_n_gene_by_counts_high': 5000,
        'filter_total_counts_high': 100000,
        'filter_mito_percentage_low': 20,
        'filter_ribo_percentage_low': 50,
        'batch_correction': False,
        'normalize_total_target_sum': 10000,
        'log1p_transform': True,
        'highly_variable_genes_num': 3000,
        'scale': True,
        'pca_comps': 50,
        'find_neighbors_neighbors_num': 15,
        'find_neighbors_using_pcs': 50,
        'unsupervised_cluster_method': 'leiden',
        'leiden_or_louvain_group_numbers': 10,
        'visualize_method': 'UMAP',
    }
    
    PROMPTS = {
        'SYSTEM_PROMPT': """You are an assistant for automating the extraction, 
        processing, and annotation of single-cell data from literature. Next, 
        I will input an article and provide instructions. Please strictly 
        adhere to the format required by the instructions and avoid outputting any unnecessary content.""",
        
        'USER_ARTICLE_PROMPT': """Please extract the following information from the article, 
        OUTPUT_FORMAT: 'FINISHED' when you have finished extracting the information,
        Here is the article:\n""",
        
        'FILTER_PROMPT': """Following the processing workflow described in the article for single-cell datasets, 
        replace the placeholders {} in the following preprocessing parameters with the values used in the article. 
        If no specific value is described in the article, use 'default' at the corresponding position. 
        If the filtering is likely not used in the article, use 'null'. Be sure to provide reasoning for the filtering parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        filter_cells_low: {int|default|null, minimum number of genes expressed in a cell, usually around 300} 
        filter_genes_low: {int|default|null, minimum number of cells expressing a gene}
        filter_n_gene_by_counts_high: {int|default|null, maximum number of genes expressed in a cell}
        filter_total_counts_high: {int|default|null, maximum total counts per cell}
        filter_mito_percentage_low: {int|default|null, maximum mitochondrial gene percentage (0,100)}
        filter_ribo_percentage_low: {int|null, maximum ribosomal gene percentage (0,100)}
        batch_correction: {bool|default, whether to perform batch correction}
        
        reasoning: {str, reasoning for the filtering parameters}""",
        
        'PREPROCESSING_PROMPT': """Continuing with the processing workflow described in the article 
        for single-cell datasets, replace the placeholders {} in the following preprocessing 
        parameters with the values used in the article. If no specific value is described in the article, 
        use 'default' at the corresponding position.
        Be sure to provide reasoning for the preprocessing parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        normalize_total_target_sum: {int|default, usually default unless especially mentioned} 
        log1p_transform: {bool|default, mostly default unless especially mentioned} 
        highly_variable_genes_num: {int|default, mostly around 2k~3k} 
        scale(centralize): {bool|default, mostly default unless especially mentioned}
        pca_comps: {int|default, mostly default unless especially mentioned}
        find_neighbors_neighbors_num: {int|default, usually default}
        find_neighbors_using_pcs: {int|default, usually default}
        
        reasoning: {str, reasoning for the preprocessing parameters}""",
        
        'CLUSTERING_PROMPT': """Now, moving on to the crucial step of cell clustering and annotation. 
        Based on your understanding of the complexity of dataset clustering in the article, 
        please provide an approximate number of clusters to be identified, along with visualization parameters.
        Be sure to provide reasoning for the clustering parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        unsupervised_cluster_method: {str: leiden or louvain}
        leiden_or_louvain_group_numbers: {int, depends on the complexity of the dataset}
        visualize_method: {str: UMAP or t-SNE, only choose one for visualization}
        
        reasoning: {str, reasoning for the clustering parameters}""",
        
        'ANNOTATION_PROMPT': """This is the output of the top 10 marker genes for each cluster:
        
        Based on gene expression and information from the article, annotate these clusters into cell types using a dictionary format.
        Please provide the 'cell type', 'certainty', 'source', 'tissue', and reasoning for each cluster.
        You may annotate different groups with the same cell type. Be sure to provide reasoning for the annotation.
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                      each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                tissue], (value chosen from [Brain, Liver, Kidney, Heart, Lung...])
                ...}. "
                
        reasoning: {str, reasoning for the annotation}""",
    }
    
    def __init__(self):
        pass
    
    @classmethod
    def get_prompt(cls, prompt_name: str) -> str:
        return cls.PROMPTS[prompt_name]
        
        
