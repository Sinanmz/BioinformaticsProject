from scipy.stats import ttest_ind
from scipy.stats import f_oneway

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import mygene


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



def t_test_in_social_setting(df, social_setting, gene):
    """
    Perform an independent two-sample t-test for a specific gene and social setting.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - social_setting (str): Social setting to compare ('KFC', 'KF', 'NF', 'ISO').
    - gene (str): Gene identifier or name.

    This function conducts an independent two-sample t-test to compare the gene expression levels
    between two study groups ('S' and 'L') within a specific social setting.

    Example:
    ```python
    import pandas as pd

    # Assuming 't_test_social_setting' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'social_setting': ['KFC', 'KF', 'NF', 'ISO', 'KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Perform t-test for a specific gene and social setting
    t_test_social_setting(gene_expression_df, 'KFC', 'ENSG00000133703')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', and 'ttest_ind' functions are defined
    and available before using this function. The specified gene and social setting must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    
    if gene is None:
        print("Gene Not Found")
        return

    subset = df[(df['social_settting'] == social_setting) & (df[gene].notna())]

    group_s = subset[subset['study_group'] == 'S'][gene]
    group_l = subset[subset['study_group'] == 'L'][gene]

    t_statistic, p_value = ttest_ind(group_s, group_l)

    print(f"Social Setting: {social_setting}")
    print(f"Gene: {annotate_gene_ids([gene], 59729).get(gene, 'NA')}")
    print(f"Ensemble ID: {gene}")
    print(f"T-Statistic: {t_statistic}")
    print(f"P-Value: {p_value}")
    
    if p_value < 0.05:
        print("The difference in means is statistically significant (p < 0.05).")
    else:
        print("The difference in means is not statistically significant (p >= 0.05).")

    effect_size = (group_s.mean() - group_l.mean()) / group_s.std()
    print(f"Effect Size: {effect_size}")
    
    if t_statistic > 0:
        print("The mean of group 'S' is greater than the mean of group 'L'.")
    elif t_statistic < 0:
        print("The mean of group 'L' is greater than the mean of group 'S'.")
    else:
        print("The means of group 'S' and group 'L' are equal.")
    print("-" * 30)




def t_test_between_social_settings(df, social_setting1, social_setting2, gene, study_group='both'):
    """
    Perform an independent two-sample t-test between two social settings for a specific gene and study group.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - social_setting1 (str): First social setting to compare.
    - social_setting2 (str): Second social setting to compare.
    - gene (str): Gene identifier or name.
    - study_group (str, optional): Study group to consider ('both', 'L', or 'S'). Default is 'both'.

    This function conducts an independent two-sample t-test to compare the gene expression levels
    between two specified social settings for a specific gene and study group.

    Example:
    ```python
    import pandas as pd

    # Assuming 't_test_between_social_settings' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'social_setting': ['KFC', 'KF', 'NF', 'ISO', 'KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Perform t-test between two social settings for a specific gene and study group
    t_test_between_social_settings(gene_expression_df, 'KFC', 'KF', 'ENSG00000133703', study_group='both')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', and 'ttest_ind' functions are defined
    and available before using this function. The specified gene and social settings must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    
    if gene is None:
        print("Gene Not Found")
        return
    
    subset1 = df[(df['social_settting'] == social_setting1) & (df[gene].notna())]
    subset2 = df[(df['social_settting'] == social_setting2) & (df[gene].notna())]

    if study_group == 'both':
        data1 = subset1[gene]
        data2 = subset2[gene]
    elif study_group == 'L':
        data1 = subset1[subset1['study_group'] == 'L'][gene]
        data2 = subset2[subset2['study_group'] == 'L'][gene]
    elif study_group == 'S':
        data1 = subset1[subset1['study_group'] == 'S'][gene]
        data2 = subset2[subset2['study_group'] == 'S'][gene]
    else:
        raise ValueError("Invalid study_group. Choose 'both', 'L', or 'S'.")

    t_statistic, p_value = ttest_ind(data1, data2)

    print(f"Social Setting 1: {social_setting1}")
    print(f"Social Setting 2: {social_setting2}")
    print(f"Study Group: {study_group}")
    print(f"Gene: {annotate_gene_ids([gene], 59729).get(gene, 'NA')}")
    print(f"Ensemble ID: {gene}")
    print(f"T-Statistic: {t_statistic}")
    print(f"P-Value: {p_value}")
    
    if p_value < 0.05:
        print("The difference in means is statistically significant (p < 0.05).")
    else:
        print("The difference in means is not statistically significant (p >= 0.05).")

    effect_size = (data1.mean() - data2.mean()) / data1.std()
    print(f"Effect Size: {effect_size}")
    
    if t_statistic > 0:
        print(f"The mean of {study_group} study group in {social_setting1} is greater than the mean in {social_setting2}.")
    elif t_statistic < 0:
        print(f"The mean of {study_group} study group in {social_setting2} is greater than the mean in {social_setting1}.")
    else:
        print(f"The means of {study_group} study group in {social_setting1} and {social_setting2} are equal.")
    print("-" * 30)



def anova_between_social_settings(df, gene, social_settings=['ISO', 'KFC', 'KF', 'NF'], study_group='both'):
    """
    Perform an Analysis of Variance (ANOVA) test for a specific gene across different social settings and study groups.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.
    - social_settings (list, optional): List of social settings to include. Default is ['ISO', 'KFC', 'KF', 'NF'].
    - study_group (str, optional): Study group to consider ('both', 'L', or 'S'). Default is 'both'.

    This function conducts an ANOVA test to assess the statistical significance of the mean differences
    in gene expression levels across different social settings and study groups.

    Example:
    ```python
    import pandas as pd

    # Assuming 'anova_test_for_gene' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'social_setting': ['KFC', 'KF', 'NF', 'ISO', 'KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Perform ANOVA test for a specific gene across different social settings and study groups
    anova_test_for_gene(gene_expression_df, 'ENSG00000133703', social_settings=['KFC', 'KF', 'NF', 'ISO'], study_group='both')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', and 'f_oneway' functions are defined
    and available before using this function. The specified gene, social settings, and study group must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    
    if gene is None:
        print("Gene Not Found")
        return

    subsets = [df[(df['social_settting'] == setting) & (df[gene].notna())] for setting in social_settings]

    if study_group == 'both':
        data = [subset[gene] for subset in subsets]
    elif study_group == 'L':
        data = [subset[subset['study_group'] == 'L'][gene] for subset in subsets]
    elif study_group == 'S':
        data = [subset[subset['study_group'] == 'S'][gene] for subset in subsets]
    else:
        raise ValueError("Invalid study_group. Choose 'both', 'L', or 'S'.")

    # Perform ANOVA test
    f_statistic, p_value = f_oneway(*data)

    # Output results
    print(f"Social Settings: {', '.join(social_settings)}")
    print(f"Study Group: {study_group}")
    print(f"Gene: {annotate_gene_ids([gene], 59729).get(gene, 'NA')}")
    print(f"Ensemble ID: {gene}")
    print(f"F-Statistic: {f_statistic}")
    print(f"P-Value: {p_value}")
    
    if p_value < 0.05:
        print("The difference in means is statistically significant (p < 0.05).")
    else:
        print("The difference in means is not statistically significant (p >= 0.05).")

    effect_size = None
    if len(data) == 2:
        mean_diff = data[0].mean() - data[1].mean()
        pooled_std = ((data[0].var() + data[1].var()) / 2) ** 0.5
        effect_size = mean_diff / pooled_std
    print(f"Effect Size: {effect_size}")
    
    print("-" * 30)


def t_test_between_study_groups(df, gene):
    """
    Perform an independent two-sample t-test between study groups for a specific gene.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.

    This function conducts an independent two-sample t-test to compare the gene expression levels
    between two study groups ('S' and 'L') for a specific gene.

    Example:
    ```python
    import pandas as pd

    # Assuming 't_test_between_study_groups' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Perform t-test between study groups for a specific gene
    t_test_between_study_groups(gene_expression_df, 'ENSG00000133703')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', and 'ttest_ind' functions are defined
    and available before using this function. The specified gene must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    if gene is None:
        print("Gene Not Found")
        return
    group_s = df[(df['study_group'] == 'S') & (df[gene].notna())][gene]
    group_l = df[(df['study_group'] == 'L') & (df[gene].notna())][gene]

    t_statistic, p_value = ttest_ind(group_s, group_l)

    print(f"Gene: {annotate_gene_ids([gene], 59729).get(gene, 'NA')}")
    print(f"Ensemble ID: {gene}")
    print(f"T-Statistic: {t_statistic}")
    print(f"P-Value: {p_value}")
    
    if p_value < 0.05:
        print("The difference in means is statistically significant (p < 0.05).")
    else:
        print("The difference in means is not statistically significant (p >= 0.05).")

    effect_size = (group_s.mean() - group_l.mean()) / group_s.std()
    print(f"Effect Size: {effect_size}")
    
    if t_statistic > 0:
        print("The mean of study group 'S' is greater than the mean of study group 'L'.")
    elif t_statistic < 0:
        print("The mean of study group 'L' is greater than the mean of study group 'S'.")
    else:
        print("The means of study group 'S' and 'L' are equal.")
    print("-" * 30)


