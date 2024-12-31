import numpy as np
import  pandas as pd
# import anndata
from anndata import AnnData
from scipy.sparse import issparse




def convert_format(adata, output_file, output_meta, index=True):
    adata.to_df().to_csv(output_file, index=index)
    metadata = adata.obs[["tissue", "donor_id"]]  # Select metadata columns
    metadata.to_csv(output_meta, index=True)

def filter_mean(adata):
    # Compute mean and standard deviation along cells (axis=0)
    gene_means = adata.X.mean(axis=0)


    # Define thresholds
    mean_threshold = 0.09

    # Create filter based on mean and std
    genes_to_keep = (gene_means >= mean_threshold)

    # Apply the filter to the adata object (filter out genes based on conditions)
    adata = adata[:, genes_to_keep]

    return adata

def sample_1000_cells_per_tissue(combined_dataset):
    n_samples_per_tissue = 1000
    np.random.seed(42)
    unique_tissues = combined_dataset.obs['tissue'].unique()
    sampled_indices = []
    for tissue in unique_tissues:
        tissue_indices = combined_dataset.obs[combined_dataset.obs['tissue'] == tissue].index
        if len(tissue_indices) > n_samples_per_tissue:
            sampled_indices.extend(np.random.choice(tissue_indices, n_samples_per_tissue, replace=False))
        else:
            sampled_indices.extend(tissue_indices)
    sampled_dataset = combined_dataset[sampled_indices]
    return sampled_dataset




def filter_adata(adata: AnnData) -> AnnData:
    """
    Filters the AnnData object by removing features (columns) based on
    mean and standard deviation thresholds.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.

    Returns
    -------
    AnnData
        A new AnnData object with filtered features.
    """
    mean_threshold = 0.00001
    std_threshold = 0.0001

    # Ensure `.X` is dense for computations
    X = adata.X.toarray() if issparse(adata.X) else adata.X

    # Compute means and standard deviations
    col_means = X.mean(axis=0)
    col_stds = X.std(axis=0)

    # Filter based on thresholds
    valid_cols = (col_means >= mean_threshold) & (col_stds >= std_threshold)

    # Create a new AnnData object with the filtered data
    filtered_adata = AnnData(
        X=X[:, valid_cols],
        obs=adata.obs.copy(),
        var=adata.var[valid_cols].copy(),
        uns=adata.uns.copy()
    )

    return filtered_adata


def combine_and_filter_tissues(adata, tissues_to_keep):
    tissue_mapping = {
        'duodenum' : 'small intestine',


        'adipose tissue': 'fat',
        'subcutaneous adipose tissue': 'fat',
        'brown adipose tissue': 'fat',

        'skin of body': 'skin',
        'skin of chest': 'skin',
        'skin of abdomen': 'skin',
        'eyelid': 'skin',


        'aorta': 'cardiovascular system',
        'vasculature': 'cardiovascular system',
        'coronary artery': 'cardiovascular system',
        'inferior vena cava': 'cardiovascular system',
        'left coronary artery': 'cardiovascular system',
        'abdominal aorta': 'cardiovascular system',


        'anterior part of tongue': 'tongue',
        'posterior part of tongue': 'tongue',

        'chorioretinal region': 'eye',
        'ocular surface region': 'eye',
        'sclera': 'eye',
        'anterior segment of eyeball': 'eye',
        'conjunctiva': 'eye',
        'posterior segment of eyeball': 'eye',
        'cornea': 'eye',

    }

    adata.obs['tissue'] = adata.obs['tissue'].replace(tissue_mapping)
    adata = adata[adata.obs['tissue'].isin(tissues_to_keep)].copy()

    return adata

def filter_out_tissues_with_small_sampleset(adata):
    min_cells = 130
    tissue_count = adata.obs['tissue'].value_counts()
    non_valid_tissue = tissue_count[tissue_count < min_cells]

    valid_tissue = tissue_count[tissue_count >= min_cells].index
    filtered_adata = adata[adata.obs['tissue'].isin(valid_tissue)].copy()
    return filtered_adata

def concat_meta_to_df(df, meta):
    df['tissue'] = meta['tissue']
    df['donor_id'] = meta['donor_id']
    return df




    # df = adata.to_df()
    # df_meta = adata.obs[['tissue', 'donor_id']]
    # df_meta.to_csv("raw_data/fibro_metadata.csv", index=False)
    # df.to_csv("raw_data/fibro_data.csv", index=False)
    file_path = "sampled_data/fibro_data_t.csv"
    # df = pd.read_csv(file_path)
    # df.insert(0, 'new_column', 0)
    # df.to_csv("raw_data/fibro_data_t.csv", index=True)
    # df = df.T

    #
    #     # Print the second row
    #     second_row = next(reader)
    #     print('Second Row:', second_row)




            # df=pd.read_csv("sampled_data/gene_expression_matrix_fibro.csv")
    # df_meta = pd.read_csv("sampled_data/metadata.csv")
    # df = concat_meta_to_df(df, df_meta)
    # tissue_df = divide_to_tissue_dataframes(df)
    # binary_tissue_df = transform_matrix_to_binary(tissue_df)
    # expression_percent = [0.01, 0.05, 0.1, 0.15]
    # for percentage in expression_percent:
    #     gene_count = gene_count_per_tissue(binary_tissue_df, percentage)
    #     plot_count_gene_per_num_of_tissues(gene_count, percentage)



