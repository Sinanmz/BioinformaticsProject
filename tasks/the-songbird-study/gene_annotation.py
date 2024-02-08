import mygene
import pandas as pd


def query_gene(df, gene):
    """
    Query the gene using MyGeneInfo if not found in the DataFrame.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing gene information.
    - gene (str): Gene name or identifier to be queried.
    
    Returns:
    - str or None: Ensembl gene ID if available, None if the gene is not found.
    
    This function checks if the specified gene is present in the DataFrame's columns.
    If the gene is not found, it queries MyGeneInfo API to retrieve the Ensembl gene ID.
    
    Example:
    ```python
    import pandas as pd
    
    # Assuming 'mygene' library is imported and installed
    
    # Create a DataFrame with gene information
    gene_df = pd.DataFrame({
        'Gene_Name': ['TP53', 'BRCA1', 'EGFR'],
        'Ensembl_ID': ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000146648']
    })
    
    # Query a gene using the function
    result = query_gene(gene_df, 'KRAS')
    print(result)  # Output: ENSG00000133703 (example Ensembl ID for KRAS)
    ```
    
    Note: Ensure the 'mygene' library is installed before using this function.
    """
    if gene not in df.columns:
        mg = mygene.MyGeneInfo()
        query_result = mg.query(gene, species=59729, fields='ensembl.gene')
        try:
            gene = query_result['hits'][0]['ensembl']['gene']
        except:
            print('Gene not found')
            return None
    return gene



def annotate_gene_ids(ensembl_ids: list, species_taxid: int):
    """
    Annotate Ensembl gene IDs with corresponding gene symbols.

    Parameters:
    - ensembl_ids (list): List of Ensembl gene IDs to annotate.
    - species_taxid (int): Taxonomic ID of the species for annotation.

    Returns:
    - gene_dict (dict): A dictionary mapping Ensembl gene IDs to their corresponding gene symbols.
    """

    species_taxid = str(species_taxid)

    mg = mygene.MyGeneInfo()

    try:
        query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', species=species_taxid, fields='symbol', returnall=True)
    except Exception as e:
        print(f"An error occurred: {e}")
        return {}
    
    gene_dict = {}
    for gene in query_result['out']:
        if 'symbol' in gene:
            gene_dict[gene['query']] = gene['symbol']
        else:
            gene_dict[gene['query']] = 'NA'

    return gene_dict



def add_gene_names_to_df(df: pd.DataFrame, species_taxid: int):
    """
    Add gene names to DataFrame columns based on Ensembl gene IDs.

    Parameters:
    - df (pd.DataFrame): DataFrame containing columns with Ensembl gene IDs.
    - species_taxid (int): Taxonomic ID of the species for annotation.

    Returns:
    - df_c (pd.DataFrame): A copy of the input DataFrame with columns renamed to include gene names.
    """

    # Extract Ensembl gene IDs from DataFrame columns
    ensembl_ids = [col for col in df.columns if col.startswith('ENST')]

    # Annotate Ensembl gene IDs with gene names using the annotate_gene_ids function
    gene_names = annotate_gene_ids(ensembl_ids, species_taxid)

    # Create a copy of the input DataFrame
    df_c = df.copy()

    # Rename columns to include gene names
    df_c.columns = [f'{gene_names[col]}' if col in gene_names else col for col in df_c.columns]

    return df_c