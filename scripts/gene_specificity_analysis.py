import matplotlib.pyplot as plt
import csv
import pandas as pd


def divide_to_tissue_dataframes(adata):
    tissue_dict = {}
    for tissue in adata.obs['tissue'].unique():
        # Select rows for the current tissue
        mask = adata.obs['tissue'] == tissue
        tissue_dict[tissue] = adata[mask, :].copy()

    return tissue_dict

def transform_matrix_to_binary(tissue_dfs):
    tissue_groups_binary = {}

    for tissue, adata_group in tissue_dfs.items():
        # Extract the expression matrix (X)
        expression_matrix = adata_group.X

        # Ensure we're working with dense matrices (if necessary)
        if hasattr(expression_matrix, 'toarray'):  # If it's a sparse matrix
            expression_matrix = expression_matrix.toarray()

        # Apply the transformation: convert values > 0 to 1 and <= 0 to 0
        binary_matrix = (expression_matrix > 0).astype(int)

        # Create a new AnnData object with the binary matrix
        binary_adata = adata_group.copy()  # Copy to retain metadata
        binary_adata.X = binary_matrix  # Update the expression matrix with the binary version

        tissue_groups_binary[tissue] = binary_adata

    return tissue_groups_binary

def gene_count_per_tissue(tissue_groups_binary, percentage):
    gene_tissue_count = {}
    gene_binary_tissue = {}

    for tissue, adata in tissue_groups_binary.items():
        # Sum across cells (rows) for each gene (columns) to get the number of cells expressing each gene
        ones_count = adata.X.sum(axis=0)  # Sum over cells (axis=0)

        # Calculate the threshold for each gene based on the percentage of cells
        num_cells = adata.shape[0]
        threshold = percentage * num_cells

        # Get the genes that are expressed in at least the required percentage of cells
        expressed_genes = (ones_count >= threshold)
        flattened_gene_names = adata.var.index.to_list()


        df_expressed_genes = pd.DataFrame([expressed_genes], columns=flattened_gene_names)
        gene_binary_tissue[tissue] = df_expressed_genes

        # For each gene expressed in this tissue, increment its tissue count
        for idx, is_expressed in enumerate(expressed_genes):
            if is_expressed:
                gene = adata.var.index[idx]  # Get the gene name
                if gene in gene_tissue_count:
                    gene_tissue_count[gene] += 1
                else:
                    gene_tissue_count[gene] = 1

    # create_tissues_list_for_specific_genes(gene_binary_tissue)

    return gene_tissue_count

def create_tissues_list_for_specific_genes(gene_binary_tissue):
    tissues_with_one = []
    genes = []
    with open("specific_analysis/specific_genes_0.05.csv", mode='r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            genes.append(row[0])

    for gene in genes:
        for tissue, adata_group in gene_binary_tissue.items():
            if gene_binary_tissue[tissue][gene].iloc[0]:
                tissues_with_one.append(tissue)
                break
    unique_tissues = list(set(tissues_with_one))
    print(unique_tissues)
    write_list_to_csv(unique_tissues, "specific_analysis/unique_specific_tissues_0.05.csv")

    write_list_to_csv(tissues_with_one, "specific_analysis/specific_tissues_0.05.csv")


def gene_counter_to_number_of_tissues(tissue_genes):
    counter_tissue = {}

    for gene, num_tissues in tissue_genes.items():
        if num_tissues in counter_tissue:
            counter_tissue[num_tissues] += 1
        else:
            counter_tissue[num_tissues] = 1
    sorted_dict = dict(sorted(counter_tissue.items()))

    gene_in_one_tissue = [key for key, value in tissue_genes.items() if value==1]
    # gene_in_all_tissues = [key for key, value in tissue_genes.items() if value==14]
    # write_list_to_csv(gene_in_all_tissues, "specific_analysis/general_genes_0.05.csv")
    write_list_to_csv(gene_in_one_tissue, "specific_analysis/specific_genes_0.05.csv")
    return sorted_dict

def write_list_to_csv(list, output_path):
    with open(output_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        for word in list:
            writer.writerow([word])

def plot_count_gene_per_num_of_tissues(gene_num_for_tissue_num, percentage):

    plt.bar(gene_num_for_tissue_num.keys(), gene_num_for_tissue_num.values(), color='#9ce7c9', edgecolor='black')
    plt.title(f"Histogram of genes expressed in a number of tissues "
              f"\nwith {percentage}% expressed genes")
    plt.xlabel('Number of tissues a gene is expressed in')
    plt.ylabel('Number of genes')
    plt.xticks(range(1, len(gene_num_for_tissue_num.keys()) + 1))
    plt.savefig("plots/specific_genes_per_num_of_tissues.png")


