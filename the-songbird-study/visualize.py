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


def calculate_heatmap_values(array):
    """
    Calculate heatmap values and counts based on social settings and study groups.
    
    Parameters:
    - array (numpy.ndarray): 2D array containing data rows, each representing a social setting and study group.
    
    Returns:
    - numpy.ndarray: Heatmap values normalized by counts.
    
    This function takes a 2D array as input, where each row represents data for a specific social setting and study group.
    The function calculates the heatmap values by summing up the values in the fourth column of the input array based on
    the specified social setting and study group. The counts of occurrences are also tracked.
    
    The resulting heatmap is normalized by the corresponding counts to obtain average values for each combination of
    social setting and study group.
    
    Example:
    ```python
    import numpy as np
    
    # Example input array
    data_array = np.array([
        ['Person1', 'KFC', 'L', 10],
        ['Person2', 'KF', 'S', 5],
        ['Person3', 'NF', 'L', 8],
        # ... more rows
    ])
    
    # Calculate heatmap values using the function
    result_heatmap = calculate_heatmap_values(data_array)
    print(result_heatmap)
    ```
    
    Note: The input array is expected to have specific columns representing social settings, study groups, and values.
    Ensure that the input array adheres to the expected structure.
    """
    heatmap = np.zeros((4, 2))
    counts = np.zeros((4, 2))

    for row in array:
        social_index = {'KFC': 0, 'KF': 1, 'NF': 2, 'ISO': 3}[row[1]]
        group_index = {'L': 0, 'S': 1}[row[2]]

        heatmap[social_index][group_index] += row[4]
        counts[social_index][group_index] += 1

    return heatmap / counts

def plot_heatmap(heatmap, gene, gene_name):
    """
    Plot the heatmap using Seaborn.
    
    Parameters:
    - heatmap (pd.DataFrame): DataFrame containing heatmap values.
    - gene (str): Gene identifier or name.
    - gene_name (str): Optional gene name for enhanced plot title.

    This function takes a DataFrame containing heatmap values and plots the heatmap using Seaborn.
    The social settings and study groups are labeled and mapped for better visualization.
    
    Parameters:
    - heatmap (pd.DataFrame): DataFrame containing heatmap values.
    - gene (str): Gene identifier or name.
    - gene_name (str): Optional gene name for enhanced plot title.

    The input heatmap DataFrame is transposed and labeled with social setting and study group names.
    Seaborn's heatmap function is used to visualize the gene expression across study groups and social settings.

    Example:
    ```python
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    # Example heatmap values and gene information
    example_heatmap = np.array([[10, 5, 8, 15], [7, 3, 12, 9]])
    example_gene = 'ENSG00000133703'
    example_gene_name = 'KRAS'

    # Plot the heatmap using the function
    plot_heatmap(example_heatmap, example_gene, example_gene_name)
    ```
    
    Note: Ensure Seaborn and Matplotlib libraries are installed before using this function.
    """
    heatmap = pd.DataFrame(heatmap).T
    heatmap.columns = ['KFC', 'KF', 'NF', 'ISO']

    social_setting_map = {
        'ISO': 'Isolated',
        'KFC': 'Known Female Continuous',
        'KF': 'Known Female Reunited',
        'NF': 'Novel Female'
    }

    heatmap.columns = [social_setting_map[x] for x in heatmap.columns]
    heatmap.index = ['L', 'S']

    sns.heatmap(heatmap, cmap='coolwarm', annot=False)
    plt.title(f'Gene Expression of {gene} {"or" if gene_name else ""} {gene_name} Across Study Groups and Social Settings')
    plt.ylabel('Study Group')
    plt.xlabel('Social Setting')
    plt.show()

def plot_violin(df, gene, gene_name):
    """
    Plot the violin plot using Seaborn.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.
    - gene_name (str): Optional gene name for enhanced plot title.

    This function generates a violin plot using Seaborn to visualize the distribution of gene expression
    across different study groups and social settings.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.
    - gene_name (str): Optional gene name for enhanced plot title.

    The social settings are mapped for better visualization, and a violin plot is created using Seaborn.
    The plot shows the distribution of gene expression for each combination of social setting and study group.

    Example:
    ```python
    import pandas as pd

    # Assuming 'plot_violin' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'social_setting': ['KFC', 'KF', 'NF', 'ISO', 'KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Plot the violin plot using the function
    plot_violin(gene_expression_df, 'ENSG00000133703', 'KRAS')
    ```

    Note: Ensure Seaborn and Matplotlib libraries are installed before using this function.
    """
    social_setting_map = {
        'ISO': 'Isolated',
        'KFC': 'Known Female Continuous',
        'KF': 'Known Female Reunited',
        'NF': 'Novel Female'
    }

    array = df[['social_settting', 'study_group', 'tissue_id', gene]]
    df['social_settting'] = df['social_settting'].map(social_setting_map)

    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.violinplot(x="social_settting", y=gene, hue="study_group", data=array, palette="muted", split=True)

    plt.title(f'Violin Plot of {gene} {"or" if gene_name != "" else ""} {gene_name} Across Study Groups and Social Settings')
    plt.ylabel(f'{gene} {"or" if gene_name != "" else ""} {gene_name} Expression')
    plt.xlabel('Social Setting')
    plt.legend(title='Study Group', loc='upper right', bbox_to_anchor=(1.2, 1))
    plt.show()

def prepare_scatter_data(df, gene1, gene2):
    """
    Prepare data for a scatter plot from the input DataFrame.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene1 (str): Identifier or name of the first gene.
    - gene2 (str): Identifier or name of the second gene.

    Returns:
    - pd.DataFrame: DataFrame with relevant columns for the scatter plot.

    This function extracts relevant columns from the input DataFrame to prepare data for a scatter plot.
    It creates a new DataFrame with columns representing study group, social setting, and the expression levels
    of the specified genes (gene1 and gene2) for each sample.

    Example:
    ```python
    import pandas as pd

    # Assuming 'prepare_scatter_data' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'sample_id': [1, 2, 3, 4],
        'social_setting': ['KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4'],
        'Gene1': [10, 5, 8, 15],
        'Gene2': [7, 3, 12, 9]
    })

    # Prepare data for a scatter plot
    scatter_data = prepare_scatter_data(gene_expression_df, 'Gene1', 'Gene2')
    print(scatter_data)
    ```

    Note: Ensure the specified genes (gene1 and gene2) are present in the DataFrame.
    """
    array = df[['sample_id', 'social_settting', 'study_group', 'tissue_id', gene1, gene2]]
    array = np.array(array)

    data = {'Study Group': [], 'Social Setting': [], f'Gene Expression - {gene1}': [], f'Gene Expression - {gene2}': []}

    for row in array:
        study_group = row[2]
        social_setting = row[1]
        expression_gene1 = row[4]
        expression_gene2 = row[5]

        data['Study Group'].append(study_group)
        data['Social Setting'].append(social_setting)
        data[f'Gene Expression - {gene1}'].append(expression_gene1)
        data[f'Gene Expression - {gene2}'].append(expression_gene2)

    return pd.DataFrame(data)


def plot_scatter(df_scatter, gene1, gene2, gene_names):
    """
    Plot the scatter plot using Seaborn.

    Parameters:
    - df_scatter (pd.DataFrame): DataFrame containing scatter plot data.
    - gene1 (str): Identifier or name of the first gene.
    - gene2 (str): Identifier or name of the second gene.
    - gene_names (dict): Dictionary mapping gene identifiers to gene names.

    This function generates a scatter plot using Seaborn to visualize the relationship between
    the expressions of two genes across different study groups and social settings.

    Parameters:
    - df_scatter (pd.DataFrame): DataFrame containing scatter plot data.
    - gene1 (str): Identifier or name of the first gene.
    - gene2 (str): Identifier or name of the second gene.
    - gene_names (dict): Dictionary mapping gene identifiers to gene names.

    The social settings are mapped for better visualization, and a scatter plot is created using Seaborn.
    Each point represents a sample, and the color and style indicate the social setting and study group, respectively.

    Example:
    ```python
    import pandas as pd

    # Assuming 'plot_scatter' function is available

    # Create a DataFrame with scatter plot data
    scatter_data = pd.DataFrame({
        'Study Group': ['L', 'S', 'L', 'S'],
        'Social Setting': ['KFC', 'KF', 'NF', 'ISO'],
        'Gene Expression - Gene1': [10, 5, 8, 15],
        'Gene Expression - Gene2': [7, 3, 12, 9]
    })

    # Gene names mapping
    gene_names_mapping = {'Gene1': 'ENSG00000133703', 'Gene2': 'ENSG00000146648'}

    # Plot the scatter plot using the function
    plot_scatter(scatter_data, 'Gene1', 'Gene2', gene_names_mapping)
    ```

    Note: Ensure Seaborn and Matplotlib libraries are installed before using this function.
    """

    social_setting_map = {
        'ISO': 'Isolated',
        'KFC': 'Known Female Continuous',
        'KF': 'Known Female Reunited',
        'NF': 'Novel Female'
    }

    gene1_name = gene_names.get(gene1, 'NA')
    gene2_name = gene_names.get(gene2, 'NA')
    if gene1_name == 'NA':
        gene1_name = ''
    if gene2_name == 'NA':
        gene2_name = ''

    df_scatter['Social Setting'] = df_scatter['Social Setting'].map(social_setting_map)
    df_scatter['Study Group'] = pd.Categorical(df_scatter['Study Group'], categories=['L', 'S'], ordered=True)
    df_scatter['Social Setting'] = pd.Categorical(df_scatter['Social Setting'], categories=social_setting_map.values(), ordered=True)

    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_scatter, x=f'Gene Expression - {gene1}', y=f'Gene Expression - {gene2}', hue='Social Setting', style='Study Group', palette='viridis', s=100)

    plt.title(f'Scatter Plot of Gene Expressions ({gene1} vs {gene2}) Across Study Groups and Social Settings')
    plt.ylabel(f'Gene Expression - {gene2} {"or" if gene2_name != "" else ""} {gene2_name}')
    plt.xlabel(f'Gene Expression - {gene1} {"or" if gene1_name != "" else ""} {gene1_name}')
    
    plt.legend(title='Social Setting')
    plt.show()



def make_heatmap(df, gene):
    """
    Creates a heatmap to visualize gene expression across different study groups and social settings.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.

    This function integrates various functionalities to create a heatmap visualizing gene expression.
    It first queries the gene using the query_gene function to ensure it is present in the DataFrame.
    If the gene is not found, the function returns without plotting the heatmap.

    The function then extracts relevant data from the DataFrame, calculates heatmap values using
    the calculate_heatmap_values function, and finally plots the heatmap using the plot_heatmap function.

    Example:
    ```python
    import pandas as pd

    # Assuming 'query_gene', 'calculate_heatmap_values', 'annotate_gene_ids', and 'plot_heatmap' functions are available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'sample_id': [1, 2, 3, 4],
        'social_setting': ['KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4'],
        'ENSG00000133703': [10, 5, 8, 15]
    })

    # Create a heatmap for gene expression
    make_heatmap(gene_expression_df, 'ENSG00000133703')
    ```

    Note: Ensure the 'query_gene', 'calculate_heatmap_values', 'annotate_gene_ids', and 'plot_heatmap' functions are defined
    and available before using this function. The gene identifier or name must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    if gene is None:
        print("Gene Not Found")
        return

    array = df[['sample_id', 'social_settting', 'study_group', 'tissue_id', gene]].to_numpy()
    heatmap = calculate_heatmap_values(array)

    gene_name = annotate_gene_ids([gene], 59729).get(gene, 'NA')
    if gene_name == 'NA':
        gene_name = ''

    plot_heatmap(heatmap, gene, gene_name)


def make_violin_plot(df, gene):
    """
    Creates a violin plot to visualize the distribution of gene expression across different study groups and social settings.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene (str): Gene identifier or name.

    This function integrates various functionalities to create a violin plot visualizing the distribution of gene expression.
    It first queries the gene using the query_gene function to ensure it is present in the DataFrame.
    If the gene is not found, the function returns without plotting the violin plot.

    The function then annotates the gene name using the annotate_gene_ids function and finally
    plots the violin plot using the plot_violin function.

    Example:
    ```python
    import pandas as pd

    # Assuming 'make_violin_plot' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'social_setting': ['KFC', 'KF', 'NF', 'ISO', 'KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S', 'L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8'],
        'ENSG00000133703': [10, 5, 8, 15, 12, 7, 9, 14]
    })

    # Create a violin plot for gene expression
    make_violin_plot(gene_expression_df, 'ENSG00000133703')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', and 'plot_violin' functions are defined
    and available before using this function. The gene identifier or name must be present in the DataFrame.
    """
    gene = query_gene(df, gene)
    if gene is None:
        print("Gene Not Found")
        return

    gene_name = annotate_gene_ids([gene], 59729).get(gene, 'NA')
    if gene_name == 'NA':
        gene_name = ''

    plot_violin(df, gene, gene_name)


def make_scatter_plot(df, gene1, gene2):
    """
    Creates a scatter plot to visualize the distribution of gene expressions across different study groups and social settings.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gene expression data.
    - gene1 (str): Identifier or name of the first gene.
    - gene2 (str): Identifier or name of the second gene.

    This function integrates various functionalities to create a scatter plot visualizing the relationship between
    the expressions of two genes across different study groups and social settings.
    It queries the genes using the query_gene function to ensure they are present in the DataFrame.
    If either gene is not found, the function returns without plotting the scatter plot.

    The function then annotates the gene names using the annotate_gene_ids function,
    prepares the data for the scatter plot using the prepare_scatter_data function,
    and finally plots the scatter plot using the plot_scatter function.

    Example:
    ```python
    import pandas as pd

    # Assuming 'make_scatter_plot' function is available

    # Create a DataFrame with gene expression data
    gene_expression_df = pd.DataFrame({
        'sample_id': [1, 2, 3, 4],
        'social_setting': ['KFC', 'KF', 'NF', 'ISO'],
        'study_group': ['L', 'S', 'L', 'S'],
        'tissue_id': ['T1', 'T2', 'T3', 'T4'],
        'ENSG00000133703': [10, 5, 8, 15],
        'ENSG00000146648': [7, 3, 12, 9]
    })

    # Create a scatter plot for gene expressions
    make_scatter_plot(gene_expression_df, 'ENSG00000133703', 'ENSG00000146648')
    ```

    Note: Ensure the 'query_gene', 'annotate_gene_ids', 'prepare_scatter_data', and 'plot_scatter' functions are defined
    and available before using this function. The specified genes (gene1 and gene2) must be present in the DataFrame.
    """
    gene1 = query_gene(df, gene1)
    gene2 = query_gene(df, gene2)

    if gene1 is None or gene2 is None:
        print("Gene Not Found")
        return

    gene_names = annotate_gene_ids([gene1, gene2], 59729)
    df_scatter = prepare_scatter_data(df, gene1, gene2)
    plot_scatter(df_scatter, gene1, gene2, gene_names)