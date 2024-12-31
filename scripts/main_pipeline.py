import numpy as np
import scanpy as sc
from pandas.core.computation.expr import intersection
from scipy.sparse import issparse

import scipy.sparse as sp


import pre_processing as pp
import gene_specificity_analysis as gsa
import batch_effect as be
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns

import anndata
import anndata as ad



def process_adata():
    adata = sc.read("data/raw_data/fibro_data.h5ad")
    tissues_to_keep = [
        "cardiovascular system", "liver", "tongue", "trachea", "thymus",
        "bladder organ", "adipose", "prostate gland",
        "skin", "eye", "large intestine", "duodenum", "ascending colon", "exocrine pancreas, combine_and_filter_tissues"
    ]
    adata = pp.combine_and_filter_tissues(adata, tissues_to_keep)
    adata = pp.filter_out_tissues_with_small_sampleset(adata)
    adata = pp.sample_1000_cells_per_tissue(adata)

    # data_df = pd.DataFrame(adata.X.toarray(),  # Convert sparse matrix to dense if necessary
    #                        index=adata.obs.index,  # Cell names
    #                        columns=adata.var.index)  # Gene names
    return adata

def add_new_tissues(adata, tissues_to_keep):
    adata = pp.combine_and_filter_tissues(adata, tissues_to_keep)
    adata = pp.filter_out_tissues_with_small_sampleset(adata)
    adata = pp.sample_1000_cells_per_tissue(adata)
    return adata


def normalization_and_log_transform(adata):
    adata = adata.copy()
    is_sparse = issparse(adata.X)

    if is_sparse:
        row_sums = adata.X.sum(axis=1).A1  # Convert sparse result to dense array
        row_sums[row_sums == 0] = 1  # Avoid division by zero for rows with all zeros
        adata.X = adata.X.multiply(1 / row_sums[:, np.newaxis])
    else:
        row_sums = adata.X.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero for rows with all zeros
        adata.X = adata.X / row_sums

    # Apply log(1 + n)
    if is_sparse:
        adata.X = adata.X.log1p()
    else:
        adata.X = np.log1p(adata.X)

    return adata


def transpose_for_sanity():
    adata = process_adata()
    matrix = adata.X.T.toarray()
    df = pd.DataFrame(
        matrix,
        index=adata.var.index,  # Gene names
        columns=adata.obs.index  # Cell names
    )
    return df

def print_shape_csv(file_path):
    # Open the CSV file to determine the shape
    # with open(file_path, newline='') as csvfile:
    #     reader = csv.reader(csvfile)
    #     rows = list(reader)  # Convert reader to a list of rows
    #     num_rows = len(rows)  # Total number of rows
    #     num_cols = len(rows[0]) if rows else 0  # Number of columns in the first row (if any)
    #
    # print(f"Shape: ({num_rows}, {num_cols})")

    # Reopen the file to iterate through rows
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)

        # Print the first row
        first_row = next(reader, None)  # Safely handle empty files
        if first_row:
            print('First Row:', first_row)
        else:
            print("The file is empty.")
            return

        # Print the first column for all subsequent rows
        print("First Column Values:")
        for row in reader:
            if row:  # Check for empty rows
                print(row[0])


def gene_specificity(adata):
    # creates dictionary for each tissue with a matching adata
    tissue_dict = gsa.divide_to_tissue_dataframes(adata)
    # transform each tissue adata to binary - 1 for gene is expressed in that cell, 0 otherwise
    tissue_groups_binary = gsa.transform_matrix_to_binary(tissue_dict)
    # creates dictionary for each gene with the number os tissues that gene is expressed in
    gene_tissue_count = gsa.gene_count_per_tissue(tissue_groups_binary, 0.05)
    #
    counter_tissue = gsa.gene_counter_to_number_of_tissues(gene_tissue_count)
    # gsa.plot_count_gene_per_num_of_tissues(counter_tissue, 0.05)



def filter_mean_std_genes(adata):
    # df = df.iloc[:, 1:]
    # df = df.T
    # be.plot_histogram(df)
    df = pp.filter_adata(adata)
    # df.to_csv("normalized_data/log_transcription_quotients_filtered.csv", index=True)
    return df

def plot_PCA(df_with_meta):

    be.plot_virtual_3D_PCA(df_with_meta)

def convert_gene_names(input_file_path, output_file_path):
    df_genes = pd.read_csv(input_file_path)
    df_genes.set_index(df_genes.columns[0], inplace=True)
    df_genes = df_genes[['name']]
    df_genes.to_csv(output_file_path, index=True, sep='\t')

def add_cols_of_meta_to_df(df, meta):
    df = df.set_index(df.columns[0])
    meta = meta.set_index(meta.columns[0])

    # meta.to_csv("data/raw_data/metadata.csv", index=True, sep="\t")
    df[['tissue', 'donor_id']] = meta[['tissue', 'donor_id']]
    return df

def filter_specific_genes_expression_matrix():
    gene_names = pd.read_csv("specific_analysis/specific_genes_0.05.csv", sep='\t')

    df = pd.read_csv("data/normalized_data/expression_matrix_normalized.csv")

    cols=df.columns
    intersection = df.columns.intersection(gene_names)
    non_intersecting = df.columns.difference(gene_names)

    filtered_df = df.loc[:, df.columns.intersection(gene_names)]
    # filtered_df.to_csv("data/normalized_data/specific_genes_expression_matrix.csv", index=True)
    return filtered_df

def filter_common_genes_expression_matrix():
    gene_name = pd.read_csv("specific_analysis/general_genes_names.csv", sep='\t')
    gene_names = gene_name['initial_alias']
    df = pd.read_csv("data/normalized_data/expression_matrix_normalized.csv")
    cols=df.columns
    intersection = df.columns.intersection(gene_names)
    non_intersecting = df.columns.difference(gene_names)

    filtered_df = df.loc[:, df.columns.intersection(gene_names)]
    filtered_df.to_csv("data/normalized_data/common_genes_expression_matrix.csv", index=True)
    return filtered_df

def print_all_tissues(adata):
    for tissue in adata.obs['tissue'].unique():
        print(tissue)

if __name__ == '__main__':
    adata = sc.read("data/raw_data/fibro_data.h5ad")
    adata = process_adata()
    print("shape of raw data: ", adata.shape)

    adata = normalization_and_log_transform(adata)

    print("finished normalization")

    # meta = adata.obs[['tissue', 'donor_id']]
    # df = adata.to_df()
    # df = pp.concat_meta_to_df(df, meta)
    plot_PCA(adata)

    print("filtering")
    adata = filter_mean_std_genes(adata)
    print("shape of filtered data: ", adata.shape)
    print("finished filtering")
    adata.write_h5ad("data/normalized_data/manual_norm/fibro_data_norm_filtered.h5ad")




    # df.to_csv("data/normalized_data/expression_matrix_normalized_smallIntestine.csv", index=True)

    # meta = pd.read_csv("data/normalized_data/metadata_smallIntestine.csv")
    # df = pd.read_csv("data/normalized_data/log_transcription_quotients.txt", sep='\t')
    # df = df.set_index(df.columns[0])
    # df = df.T
    # df.to_csv("data/normalized_data/expression_matrix_normalized_smallIntestine_for_Sanity.csv")
    # print(df.shape)
    # df = pp.concat_meta_to_df(df, meta)
    # df = filter_mean_std_genes(df)
    # print(df.shape)
    #
    # df.to_csv("data/normalized_data/expression_matrix_normalized_smallIntestine_filtered.csv")
    # df = pd.read_csv("data/normalized_data/log_transcription_quotients_filtered.csv")






