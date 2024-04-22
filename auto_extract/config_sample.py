class Config:
    API_KEY = 'YOUR_API_KEY'
    API_BASE_URL = "YOUR_API_BASE_URL" # for third party openai api
    TYPE = "openai" # claude or openai
    MODEL = "claude-3-sonnet-20240229" # model processing the article
    TOOL_MODEL = "claude-3-opus-20240229" # model for short messages
    
    LIST_TYPE_PARAMS = ['genes_to_query']
    INT_TYPE_PARAMS = ['filter_cells_min_genes', 'filter_cells_max_genes', 'filter_cells_min_counts', 'filter_cells_max_counts',
                          'filter_genes_low', 'filter_mito_percentage_low', 'filter_ribo_percentage_low',
                          'normalize_total_target_sum', 'highly_variable_genes_num', 'pca_comps',
                          'find_neighbors_neighbors_num', 'find_neighbors_using_pcs', 'leiden_or_louvain_group_numbers']
    BOOL_TYPE_PARAMS = ['batch_correction', 'log1p_transform', 'scale']
    CATEGORICAL_PARAMS = ['unsupervised_cluster_method', 'visualize_method']
    
    DEFAULT_PARAMS = {
        'filter_cells_min_genes': 300,
        'filter_cells_max_genes': 5000,
        'filter_cells_min_counts': None,
        'filter_cells_max_counts': 100000,
        'filter_genes_low': 3,
        'filter_mito_percentage_low': 20,
        'filter_ribo_percentage_low': 50,
        'batch_correction': True,
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
        'reannotation': True,
        'max_genes_to_query': 15,
        'genes_to_query': [],
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
        filter_cells_min_genes: {int|default|null, minimum number of genes expressed in a cell, usually around 300}
        filter_cells_max_genes: {int|default|null, maximum number of genes expressed in a cell, usually around 5000}
        filter_cells_min_counts: {int|default|null, minimum allowed total counts per cell usually null}
        filter_cells_max_counts: {int|default|null, maximum allowed total counts per cell}
        filter_genes_low: {int|default|null, minimum number of cells expressing a gene}
        filter_mito_percentage_low: {int|default|null, maximum mitochondrial gene percentage (0,100)}
        filter_ribo_percentage_low: {int|null, maximum ribosomal gene percentage (0,100)}
        
        reasoning: {str, reasoning for the filtering parameters}""",
        
        'PREPROCESSING_PROMPT': """Continuing with the processing workflow described in the article 
        for single-cell datasets, replace the placeholders {} in the following preprocessing 
        parameters with the values used in the article. If no specific value is described in the article, 
        use 'default' at the corresponding position.
        Be sure to provide reasoning for the preprocessing parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        normalize_total_target_sum: {int|default, usually default unless especially mentioned} 
        log1p_transform: {bool|default, mostly default unless especially mentioned}
        batch_correction: {bool|default, whether to perform batch correction, mostly default unless especially mentioned}
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
        
        Based on gene expression and the detailed discussion from the article, annotate these clusters into cell types using a dictionary format.
        Please provide the 'cell type', 'certainty', 'source', 'tissue', and reasoning for each cluster.
        You may annotate different groups with the same cell type. You should try to assign a **cell ontology** label to each cluster (e.g. B cell, T cell, etc.),
        with modification to make your annotations more concordant with the original paper  (e.g. 'CD4+ T cell' or 'T cell 2').
        If you cannot tell the cell type, name it as 'Unknown'. Be sure to provide reasoning for the annotation.
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                      each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                tissue], (value chosen from [Brain, Liver, Kidney, Heart, Lung...])
                ...}. "
                
        reasoning: {str, reasoning for the annotation}""",
        
        'REVIEW_PROMPT': """To refine your annotation for cluster, you can decide a list of genes to query their expression raw dataset. 
        These expression data should helps you decide your low confidence group annotation and increase certainty. Remind that:
        1. You don't need to change your high-confident annotation in most cases
        2. For those clusters you are unassure, you can query additional classical markers of the annotated cell type, 
        based on your biological knowledge. For those clusters you label as 'Unknown', you can query the most specific markers of remaining cell types.
        3. The max length of the gene list is """ + str(DEFAULT_PARAMS['max_genes_to_query']) + """ If there is no need to query, please output empty list
        
        OUTPUT_FORMAT:
        genes_to_query: list of genes to query, e.g. ['gene1', 'gene2', 'gene3', ...]
        
        reasoning: {str, reasoning for the genes to query}""",
        
        'REANNOTATION_PROMPT': """Based on the gene expression data queried and previous annotation, please re-annotate the clusters into cell types using a dictionary format:
    
        If you are still unsure about the cell type, you can mark it as 'Unknown'. Be sure to provide reasoning for the re-annotation. You can also change
        the previous annotation if you think it is necessary, even if you are confident about the previous annotation.
    
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                tissue], (value chosen from [Brain, Liver, Kidney, Heart, Lung...])
                ...}. "
                
        reasoning: {str, reasoning for the re-annotation}""",
    }
    
    TOOL_PROMPTS = {
        'SUMMARY_QUERY_EXPRESSION': """This is the expression data of certain genes in a single-cell dataset, with the format {gene: [exp_in_cluster_i, ...]}.
        Please indicate the expression level of each gene in cluster_i in order, and summarize their expression states, focusing on specific outlier expressions. 
        For example: geneA: highest expression in cluster 0, not expressed in other clusters. geneB: not expressed in cluster 4, relatively high expression in cluster 2. 
        Please directly output your summary by gene. """,
        
        'NO_CONTEXT_ANNOTATION_PROMPT': """Identify cell types of cells using the following markers separately for each row. Only provide the cell type name. 
        Do not show numbers before the name. Some can be a mixture of multiple cell types.
         
        OUTPUT_FORMAT (description of the parameters is in the curly braces, do not include the description in the output):
        annotation_dict: {0: 'cell_type1', 1: 'cell_type2', ...}
        
        reasoning: {str, reasoning for the annotation}""",
    }
    
    def __init__(self):
        pass
    
    @classmethod
    def get_prompt(cls, prompt_name: str) -> str:
        return cls.PROMPTS[prompt_name]
    
    @classmethod
    def get_tool_prompt(cls, prompt_name: str) -> str:
        return cls.TOOL_PROMPTS[prompt_name]