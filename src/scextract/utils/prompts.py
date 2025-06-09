class Prompts:
    # Parameters
    LIST_TYPE_PARAMS = ['genes_to_query']
    INT_TYPE_PARAMS = ['filter_cells_min_genes', 'filter_cells_max_genes', 'filter_cells_min_counts', 'filter_cells_max_counts',
                          'filter_genes_min_cells', 'filter_mito_percentage_min', 'filter_mito_percentage_max', 'filter_ribo_percentage_min', 'filter_ribo_percentage_max',
                          'normalize_total_target_sum', 'highly_variable_genes_num', 'pca_comps',
                          'find_neighbors_neighbors_num', 'find_neighbors_using_pcs', 'leiden_or_louvain_group_numbers']
    BOOL_TYPE_PARAMS = ['batch_correction', 'log1p_transform', 'scale']
    CATEGORICAL_PARAMS = ['unsupervised_cluster_method', 'visualize_method']
    
    PROMPTS = {
        'SYSTEM_PROMPT': """You are an assistant for automating the extraction, 
        processing, and annotation of single-cell data from literature. Next, 
        I will input an article and provide instructions. Please strictly 
        adhere to the format required by the instructions and avoid outputting any unnecessary content.""",
        
        'USER_ARTICLE_PROMPT': """Please extract the following information from the article, 
        OUTPUT_FORMAT: 'FINISHED' when you have finished extracting the information,
        Here is the article:\n""",
        
        'GET_METADATA_PROMPT': """Please extract the basic information and metadata of study samples from the article, 
        including the title, author, magazine name, number of samples and their characteristics, total number of single cells analyzed, raw data source,
        replace the placeholders {} in the following parameters with the values used in the article. If no specific value is described in the article,
        use 'N/A' at the corresponding position. 
        
        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        title: {str, title of the article}
        author: {str, family name of the first author plus 'et al.'}
        journal_name: {str, name of the magazine}
        sample_description: {str, short description of the samples used in scRNA-seq and their characteristics}
        total_cells: {int, total number of single cells analyzed}
        raw_data_source: {str, accession number or link to the raw data, e.g., GEO, SRA, etc.}
        
        This is an example response:
        <example>        
        <response>
        title: Single-cell RNA-seq analysis of human liver immune cells reveals a novel subset of liver-resident natural killer cells
        author: Wu et al.
        journal_name: Nature Communications
        sample_description: liver immune cells from 5 healthy donors and 3 patients with chronic hepatitis B
        total_cells: 10000
        raw_data_source: GSE123456
        </response>
        </example>""",
        
        'SUMMARY_ARTICLE_PROMPT': """Please summarize the article in a few sentences, including the main findings,
        the methodology used, and any key insights or conclusions. Summarize the article within 5 sentences.
        
        OUTPUT_FORMAT: 
        summary: {str, summary of the article}
        
        This is an example response:
        <example>
        <response>
        summary: The article presents a single-cell RNA-seq analysis of human liver immune cells...
        </response>
        </example>""",
        
        'FILTER_PROMPT': """Following the processing workflow described in the article for single-cell datasets, 
        replace the placeholders {} in the following preprocessing parameters with the values used in the article. 
        If no specific value is described in the article, use 'default' at the corresponding position. 
        If the filtering is likely not used in the article, use 'null'. Be sure to provide reasoning for the filtering parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        filter_cells_min_genes: {int|default|null, minimum number of genes expressed in a cell, usually around 300}
        filter_cells_max_genes: {int|default|null, maximum number of genes expressed in a cell, usually around 5000}
        filter_cells_min_counts: {int|default|null, minimum allowed total counts per cell usually null}
        filter_cells_max_counts: {int|default|null, maximum allowed total counts per cell}
        filter_genes_min_cells: {int|default|null, minimum number of cells expressing a gene}
        filter_mito_percentage_min: {int|default|null, minimum mitochondrial gene percentage (0,100), usually null}
        filter_mito_percentage_max: {int|default|null, maximum mitochondrial gene percentage (0,100)}
        filter_ribo_percentage_min: {int|null, minimum ribosomal gene percentage (0,100), usually null}
        filter_ribo_percentage_max: {int|null, maximum ribosomal gene percentage (0,100), usually null}
        
        reasoning: {str, reasoning for the filtering parameters}
        
        This is an example of how to extract filter_cells_min_genes, other arguments are extracted in the same way:
        <example>
        <text>
        We filter out cells expressing fewer than 300 genes
        </text>
        The output should be:
        <response>
        filter_cells_min_genes: 300
        
        reasoning: The article mentions that 'We filter out cells expressing fewer than 300 genes'
        </response>
        </example>""",
        
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
        
        reasoning: {str, reasoning for the preprocessing parameters}
        
        This is an example of how to extract highly_variable_genes_num, other arguments are extracted in the same way:
        <example>
        <text>
        We use function FindVariableFeatures to identify 2000 highly variable genes
        </text>
        The output should be:
        <response>
        highly_variable_genes_num: 2000
        </response>
        </example>""",
        
        'CLUSTERING_PROMPT': """Now, moving on to the crucial step of cell clustering and annotation. 
        Based on your understanding of the complexity of dataset clustering in the article, 
        please provide an approximate number of clusters to be identified, along with visualization parameters.
        Be sure to provide reasoning for the clustering parameters.

        OUTPUT_FORMAT(description of the parameters is after the colon, do not include the description in the output):
        unsupervised_cluster_method: {str, choose from leiden | louvain}
        leiden_or_louvain_group_numbers: {int, depends on the complexity of the datasetm, you can infer the number from the article if not mentioned}
        visualize_method: {str: UMAP or t-SNE, only choose one for visualization}
        
        reasoning: {str, reasoning for the clustering parameters}""",
        
        'ANNOTATION_PROMPT': """This is the output of the top 10 marker genes for each cluster:
        {authors_defined_marker_genes}
        Based on gene expression and the detailed discussion from the article, annotate these clusters into cell types using a dictionary format.
        Please provide the 'cell type', 'certainty', 'source' and reasoning for each cluster.
        You may annotate different groups with the same cell type. You should try to assign a **cell ontology** label to each cluster (e.g. B cell, T cell, etc.),
        with modification to make your annotations more concordant with the original paper  (e.g. 'CD4+ T cell' or 'T cell 2'). Some can be a mixture of multiple cell types.
        If you cannot tell the cell type, name it as 'Unknown'. Be sure to provide reasoning for the annotation.
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                ...}. "
                
        reasoning: {str, reasoning for the re-annotation}
        
        <example>
        <response>
        annotation_dict: {0: ['T cell', 'High', 'Article-defined'],
                            1: ['B cell', 'Medium', 'Knowledge-based']}
                            
        reasoning: The expression of CD3E and CD3D is high in cluster 0, which is a typical marker of T cells. The expression of CD19 is medium high in cluster 1, which is a typical marker of B cells, but not as high as in cluster 2
        </response>
        </example>""",
        
        'REVIEW_PROMPT': """To refine your annotation for cluster, you can decide a list of genes to query their expression raw dataset. 
        These expression data should helps you decide your low confidence group annotation and increase certainty. Remind that:
        1. You don't need to change your high-confident annotation in most cases
        2. For those clusters you are unassure, you can query additional classical markers of the annotated cell type, 
        based on your biological knowledge. For those clusters you label as 'Unknown', you can query the most specific markers of remaining cell types.
        3. The max length of the gene list is 20. If there is no need to query, please output empty list
        
        OUTPUT_FORMAT:
        genes_to_query: list of genes to query, e.g. ['gene1', 'gene2', 'gene3', ...]
        
        reasoning: {str, reasoning for the genes to query}
        
        <example>
        <response>
        genes_to_query: ['CD8A', 'CD14', 'CD34']

        reasoning: These genes are classical markers for CD8+ T cells, monocytes/macrophages, and hematopoietic stem cells, respectively. Querying their expression levels can help refine the annotation of clusters showing ambiguous marker profiles, such as Cluster 2 labeled as 'Unknown'.
        </response>
        </example>""",
        
        'REANNOTATION_PROMPT': """Based on the gene expression data queried and previous annotation, please re-annotate the clusters into cell types using a dictionary format:
    
        If you are still unsure about the cell type, you can mark it as 'Unknown'. Be sure to provide reasoning for the re-annotation. You can also change
        the previous annotation if you think it is necessary.
    
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                ...}. "
                
        reasoning: {str, reasoning for the re-annotation}
        
        <example>
        <response>
        annotation_dict: {0: ['T cell', 'High', 'Article-defined'],
                            1: ['B cell', 'Medium', 'Knowledge-based']}
                            
        reasoning: The expression of CD3E and CD3D is high in cluster 0, which is a typical marker of T cells. The expression of CD19 is medium high in cluster 1, which is a typical marker of B cells, but not as high as in cluster 2
        </response>
        </example>""",

        'ANNOTATION_PROMPT_VERBOSE': """This is the output of the top 10 marker genes for each cluster:
        {authors_defined_marker_genes}
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
                tissue, (value chosen from [Brain, Liver, Kidney, Heart, Lung...])
                disease, (value chosen from [Healthy, Cancer, Alzheimer, Parkinson...])
                developmental_stage], (value chosen from [Embryonic, Fetal, Neonatal, Adult...])
                ...}. "
                
        reasoning: {str, reasoning for the re-annotation}
        
        <example>
        <response>
        annotation_dict: {0: ['T cell', 'High', 'Article-defined', 'Blood', 'Healthy', 'Adult'],
                            1: ['B cell', 'Medium', 'Knowledge-based', 'Blood', 'Healthy', 'Adult']}
                            
        reasoning: The expression of CD3E and CD3D is high in cluster 0, which is a typical marker of T cells. The expression of CD19 is medium high in cluster 1, which is a typical marker of B cells, but not as high as in cluster 2
        </response>
        </example>""",

        'REANNOTATION_PROMPT_VERBOSE': """Based on the gene expression data queried and previous annotation, please re-annotate the clusters into cell types using a dictionary format:
    
        If you are still unsure about the cell type, you can mark it as 'Unknown'. Be sure to provide reasoning for the re-annotation. You can also change
        the previous annotation if you think it is necessary.
    
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output,
                each value in the [] should be quoted so that it is clear that it is a string value):
        annotation_dict: {0: [cell_type, 
                certainty, (value chosen from [Low, Medium, High])
                source, (value chosen from [Article-defined, Knowledge-based])
                tissue, (value chosen from [Brain, Liver, Kidney, Heart, Lung...])
                disease, (value chosen from [Healthy, Cancer, Alzheimer, Parkinson...])
                developmental_stage], (value chosen from [Embryonic, Fetal, Neonatal, Adult...])
                ...}. "
                
        reasoning: {str, reasoning for the re-annotation}
        
        <example>
        <response>
        annotation_dict: {0: ['T cell', 'High', 'Article-defined', 'Blood', 'Healthy', 'Adult'],
                            1: ['B cell', 'Medium', 'Knowledge-based', 'Blood', 'Healthy', 'Adult']}
                            
        reasoning: The expression of CD3E and CD3D is high in cluster 0, which is a typical marker of T cells. The expression of CD19 is medium high in cluster 1, which is a typical marker of B cells, but not as high as in cluster 2
        </response>
        </example>""",
        
        'CHOOSE_MODEL_PROMPT': """You are celltypist, an automated cell type annotation tool to annotate your single-cell data.
        Based on the following models, please choose the model that best fits your data:
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output):
        model: {str, model name}
        reasoning: {str, reasoning for the model choice}
        
        <example>
        <response>
        model: 'Immune_All_Low.pkl'
        reasoning: This model is trained on immune cells and is suitable for the immune cell dataset in the article.
        </response>
        </example>""",
        
        'GET_TISSUE_NAME_PROMPT': """Please provide the tissue name the dataset is derived from. Only provide the tissue name without any additional information.""",
        
        'ADD_CELLTYPE_DESCRIPTION_PROMPT': """Please provide short descriptions for cell types you have annotated: 
        
        {celltype_list}
        
        Descriptions should be concise and informative, and can include the function, location, or any other relevant information about the cell type in the context of the dataset.
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output):
        celltype_descriptors: {cell_type1: 'description1', cell_type2: 'description2', ...}
        
        <example>
        <response>
        celltype_descriptors: {'T cell': 'T cells block the growth of breast cancer cells by releasing cytotoxic granules',
                                'B cell': 'B cells produce antibodies to neutralize pathogens and prevent infection', ...}
        </response>
        </example>""",
        
        'ADD_CELLTYPE_FUNCTION_PROMPT': """Please provide the short function discussed in the article for each of the cell types you have annotated:
        
        {celltype_list}
        
        If the function is not mentioned in the article, use 'N/A' at the corresponding position.
        
        OUTPUT_FORMAT(description of the parameters is in the curly braces, do not include the description in the output):
        celltype_functions: {cell_type1: 'function1', cell_type2: 'function2', ...}
        
        <example>
        <response>
        celltype_functions: {'T cell': 'T cells block the growth of breast cancer cells by releasing cytotoxic granules',
                            'B cell': 'B cells produce antibodies to neutralize pathogens and prevent infection', ...}
        </response>
        </example>""",
    }
    
    TOOL_PROMPTS = {
        'SUMMARY_QUERY_EXPRESSION': """This is the expression data of certain genes in a single-cell dataset, with the format {gene: [exp_in_cluster_i, ...]}.
        Please indicate the expression level of each gene in cluster_i in order, and summarize their expression states, focusing on specific outlier expressions. 
        For example: geneA: highest expression in cluster 0, not expressed in other clusters. geneB: not expressed in cluster 4, relatively high expression in cluster 2. 
        Please directly output your summary by gene. """,
        
        'SUMMARY_MARKER_GENES_FROM_SHEET_PROMPT': """This is the marker genes data from Excel sheet provided by the authors
        
        {marker_genes_sheet}
        
        Please summarize the marker genes in the following format:
        celltype1: marker1, marker2, ...
        celltype2: marker1, marker2, ...
        ...
        """,
        
        'GPTCELLTYPE_ANNOTATION_PROMPT': """Identify cell types of {tissuename} cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. \n""",
        
        'NO_CONTEXT_ANNOTATION_PROMPT': """Identify cell types of cells using the following markers separately for each row. Only provide the cell type name.
        You should try to assign a **cell ontology** label to each cluster (e.g. B cell, T cell, etc.),  
        Do not show numbers before the name. Some can be a mixture of multiple cell types. 
         
        OUTPUT_FORMAT (description of the parameters is in the curly braces, do not include the description in the output):
        annotation_dict: {0: 'cell_type1', 1: 'cell_type2', ...}
        
        reasoning: {str, reasoning for the annotation}""",
        
        'MAJOR_VOTE_ANNOTATION_PROMPT': """I will give you the top 5 most abundant cell types annotation of each cluster. Based on the majority vote,
        assign the cell type to each cluster. If you are unsure, mark it as 'Unknown'. You should try to assign a **cell ontology** label to each cluster (e.g. T helper cell, etc.),

        OUTPUT_FORMAT (description of the parameters is in the curly braces, do not include the description in the output):
        annotation_dict: {0: 'cell_type1', 1: 'cell_type2', ...}
        reasoning: {str, reasoning for the annotation}

        <example>
        Cluster 0:
        cell_type                         cell_count    percentage    dataset_support
        Fibroblast (FB)                      46573        13.2              3
        Secretory-papillary fibroblast        4977         7.0               1
        Fibroblast                            3834         5.4               2
        Secretory-reticular fibroblast        2847         4.0               1
        Mesenchymal fibroblast                 936         1.3               1

        Cluster 1:
        cell_type                         cell_count    percentage    dataset_support
        T cell                               46573        45.5              4
        T helper cell                         4977         8.4              1
        T cytotoxic cell                      3834         6.5               1
        T regulatory cell                     2847         4.8               1
        T cell 2                              936         1.6               1
        ...
        <response>
        annotation_dict: {0: 'Fibroblast', 1: 'T cell', ...}

        reasoning: {Cluster 0 is annotated as 'Fibroblast' because the majority of cells (65.2%) are labeled as Fibroblast (FB), and other top annotations are also fibroblast subtypes. This annotation is well-supported across 12 datasets. Cluster 1 is clearly a T cell cluster with 78.5% of cells annotated as T cells and remaining populations being T cell subtypes, with strong support across 15 datasets...}
        </response>
        </example>

        This is the output of the top 5 marker genes for each cluster:
        """,
        'REFORMAT_RESPONSE_PROMPT': """The response format is incorrect. Please reformat the response to the correct format. The response should be strictly like one of the <response> blocks in the example below:
        <response>
        annotation_dict: {0: ['T cell', 'High', 'Article-defined'],
                            1: ['B cell', 'Medium', 'Knowledge-based']}
                            
        reasoning: The expression of CD3E and CD3D is high in cluster 0, which is a typical marker of T cells. The expression of CD19 is medium high in cluster 1, which is a typical marker of B cells, but not as high as in cluster 2
        </response>
        <response>
        filter_cells_min_genes: 300
        ...
        reasoning: The article mentions that 'We filter out cells expressing fewer than 300 genes'
        </response>
        <response>
        genes_to_query: ['CD8A', 'CD14', 'CD34']

        reasoning: These genes are classical markers for CD8+ T cells, monocytes/macrophages, and hematopoietic stem cells, respectively. Querying their expression levels can help refine the annotation of clusters showing ambiguous marker profiles, such as Cluster 2 labeled as 'Unknown'.
        </response>
        """,
    }
    
    def __init__(self):
        pass
    
    @classmethod
    def get_prompt(cls, prompt_name: str) -> str:
        return cls.PROMPTS[prompt_name]
    
    @classmethod
    def get_tool_prompt(cls, prompt_name: str) -> str:
        return cls.TOOL_PROMPTS[prompt_name]