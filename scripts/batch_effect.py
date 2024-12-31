import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import plotly.express as px
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.datasets import load_iris
import kaleido



def plot_histogram(df):
    print("plotting histogram")
    df = df.drop(columns=['tissue', 'donor_id'])

    gene_means = df.mean(axis=0)
    gene_stds = df.std(axis=0)
    print("shape of std vector", gene_stds.shape)
    sns.histplot(gene_means, bins=30, kde=True, color='blue')
    plt.title('Histogram of Gene Expression Means', fontsize=16)
    plt.xlabel('Mean Expression', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    # plt.tight_layout()
    plt.savefig('plots/histogram_mean.png')
    plt.clf()

    # Plot the histogram of gene std
    sns.histplot(gene_stds, bins=50, kde=True, color='blue')  # Use kde=True for a smooth curve overlay
    plt.title('Histogram of Gene Expression std', fontsize=16)
    plt.xlabel('Std', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.savefig('plots/histogram_std.png')

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D

def plot_3D_PCA(df):
    # Select features for PCA
    features = df.drop(columns=['tissue', 'donor_id', 'assay'])

    # Standardize the data
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    # Perform PCA
    pca = PCA(n_components=50)
    pca_result = pca.fit_transform(scaled_features)

    explained_variance = pca.explained_variance_ratio_
    pca_result=pca_result[:, :3]

    # Create a DataFrame for PCA results
    df_pca = pd.DataFrame(pca_result, columns=['PCA1', 'PCA2', 'PCA3'])
    df_pca['tissue'] = df['tissue'].values
    df_pca['assay'] = df['assay'].values
    df_pca['donor_id'] = df['donor_id'].values


    # Generate color mappings for assays and donors
    tissue_colors = {
        tissue: color for tissue, color in
        zip(df_pca['tissue'].unique(), plt.cm.tab20.colors[:len(df_pca['tissue'].unique())])
    }
    assay_colors = {
        assay: color for assay, color in zip(df_pca['assay'].unique(), plt.cm.tab20.colors[:len(df_pca['assay'].unique())])
    }
    donor_colors = {
        donor: color for donor, color in zip(df_pca['donor_id'].unique(), plt.cm.tab20.colors[:len(df_pca['donor_id'].unique())])
    }

    # 3D PCA Plot Colored by tissue
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for tissue, color in tissue_colors.items():
        subset = df_pca[df_pca['tissue'] == tissue]
        ax.scatter(subset['PCA1'], subset['PCA2'], subset['PCA3'], label=tissue, color=color, s=50, alpha=0.7)

    ax.set_title('3D PCA of Gene Expression (Colored by Tissue)', fontsize=14)
    ax.set_xlabel('PCA1', fontsize=12)
    ax.set_ylabel('PCA2', fontsize=12)
    ax.set_zlabel('PCA3', fontsize=12)
    ax.legend(title="Tissue", fontsize=8, title_fontsize=12)
    plt.tight_layout()
    plt.savefig("plots/PCA_tissues_3D.png")

    # 3D PCA Plot Colored by Assay
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for assay, color in assay_colors.items():
        subset = df_pca[df_pca['assay'] == assay]
        ax.scatter(subset['PCA1'], subset['PCA2'], subset['PCA3'], label=assay, color=color, s=50, alpha=0.7)

    ax.set_title('3D PCA of Gene Expression (Colored by Assay)', fontsize=14)
    ax.set_xlabel('PCA1', fontsize=12)
    ax.set_ylabel('PCA2', fontsize=11)
    ax.set_zlabel('PCA3', fontsize=12)
    ax.legend(title="Assay", fontsize=8, title_fontsize=12)
    plt.tight_layout()
    plt.savefig("plots/PCA_assays_3D.png")

    # 3D PCA Plot Colored by Donor ID
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for donor, color in donor_colors.items():
        subset = df_pca[df_pca['donor_id'] == donor]
        ax.scatter(subset['PCA1'], subset['PCA2'], subset['PCA3'], label=donor, color=color, s=50, alpha=0.7)

    ax.set_title('3D PCA of Gene Expression (Colored by Donor ID)', fontsize=14)
    ax.set_xlabel('PCA1', fontsize=12)
    ax.set_ylabel('PCA2', fontsize=12)
    ax.set_zlabel('PCA3', fontsize=12)
    ax.legend(title="Donor ID", fontsize=8, title_fontsize=12)
    plt.tight_layout()
    plt.savefig("plots/PCA_donors_3D.png")

    plt.figure(figsize=(8, 6))
    plt.plot(range(1, 51), explained_variance, marker='o', linestyle='-', color='blue')
    plt.title('Explained Variance by PCA Components', fontsize=14)
    plt.xlabel('Number of Principal Components', fontsize=12)
    plt.ylabel('Explained Variance Ratio', fontsize=12)
    plt.xticks(range(1, 51), fontsize=6)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("plots/PCA_explained_variance.png")


def plot_PCA_with_explained_variance_curve(df):
    # Select features for PCA
    features = df.drop(columns=['tissue', 'donor_id', 'assay'])

    # Standardize the data
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_features)


    # Create a DataFrame for PCA results
    pca_columns = ["PC1", "PC2"]
    df_pca = pd.DataFrame(pca_result, columns=pca_columns)
    df_pca['tissue'] = df['tissue'].values
    df_pca['donor_id'] = df['donor_id'].values
    df_pca['assay'] = df['assay'].values

    # Map each tissue to a distinct color
    tissue_colors = {
        tissue: color for tissue, color in zip(df_pca['tissue'].unique(), plt.cm.tab10.colors[:len(df_pca['tissue'].unique())])
    }
    donors_colors = {
        tissue: color for tissue, color in
        zip(df_pca['tissue'].unique(), plt.cm.tab10.colors[:len(df_pca['tissue'].unique())])
    }
    assays_colors = {
        tissue: color for tissue, color in
        zip(df_pca['tissue'].unique(), plt.cm.tab10.colors[:len(df_pca['tissue'].unique())])
    }
    df_pca['TissueColor'] = df_pca['tissue'].map(tissue_colors)
    df_pca['DonorColor'] = df_pca['donor_id'].map(tissue_colors)
    df_pca['AssayColor'] = df_pca['assay'].map(tissue_colors)


    plt.figure(figsize=(8, 6))
    for tissue, color in tissue_colors.items():
        subset = df_pca[df_pca['tissue'] == tissue]
        plt.scatter(subset['PCA1'], subset['PCA2'], label=tissue, color=color, s=50, alpha=0.7)

    plt.legend(title="Tissue", loc='lower left', fontsize=8, title_fontsize=12)
    plt.title('2D PCA of Gene Expression (Colored by Tissue)', fontsize=14)
    plt.xlabel('PCA1', fontsize=12)
    plt.ylabel('PCA2', fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("PCA_tissues_2D.png")

    plt.figure(figsize=(8, 6))
    for donor, color in donors_colors.items():
        subset = df_pca[df_pca['donor_id'] == donor]
        plt.scatter(subset['PCA1'], subset['PCA2'], label=donor, color=color, s=50, alpha=0.7)

    plt.legend(title="Donors", loc='lower left', fontsize=8, title_fontsize=12)
    plt.title('2D PCA of Gene Expression (Colored by donor)', fontsize=14)
    plt.xlabel('PCA1', fontsize=12)
    plt.ylabel('PCA2', fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("PCA_donors_2D.png")

    plt.figure(figsize=(8, 6))
    for assay, color in assays_colors.items():
        subset = df_pca[df_pca['assay'] == assay]
        plt.scatter(subset['PCA1'], subset['PCA2'], label=assay, color=color, s=50, alpha=0.7)

    plt.legend(title="Assay", loc='lower left', fontsize=8, title_fontsize=12)
    plt.title('2D PCA of Gene Expression (Colored by Assay)', fontsize=14)
    plt.xlabel('PCA1', fontsize=12)
    plt.ylabel('PCA2', fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("PCA_assays_2D.png")


def plot_virtual_3D_PCA(adata, tissue_filter=False):
    print("Plotting 3D PCA")

    if tissue_filter:
        adata = adata[adata.obs['tissue'].isin(['trachea', 'liver', 'bladder organ', 'large intestine', 'cardiovascular system', 'adipose']), :]

    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(adata.X)  # adata.X is the data matrix
    adata.obsm['X_pca'] = X_pca  # Store the PCA results in 'obsm' for easy access

    # tissue coloring
    fig = px.scatter_3d(
        adata.obs,
        x=adata.obsm['X_pca'][:, 0],
        y=adata.obsm['X_pca'][:, 1],
        z=adata.obsm['X_pca'][:, 2],
        color='tissue',
        title='3D PCA with Tissue Type'
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(
        legend=dict(
            itemclick='toggleothers',
            itemdoubleclick='toggle',  # These ensure interactivity with the legend
            tracegroupgap=5,  # Optional, adjusts spacing between legend items
            itemsizing='constant'  # Keeps the legend symbol size constant regardless of point size
        )
    )
    if tissue_filter:
        fig.write_html("plots/pca_3d_tissues_filter.html")
    else:
        fig.write_html("plots/pca_3d_tissues_no_norm.html")

    # donor coloring
    fig = px.scatter_3d(
        adata.obs,
        x=adata.obsm['X_pca'][:, 0],
        y=adata.obsm['X_pca'][:, 1],
        z=adata.obsm['X_pca'][:, 2],
        color='donor_id',
        title='3D PCA with Donor ID'
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(
        legend=dict(
            itemclick='toggleothers',
            itemdoubleclick='toggle',  # These ensure interactivity with the legend
            tracegroupgap=5,  # Optional, adjusts spacing between legend items
            itemsizing='constant'  # Keeps the legend symbol size constant regardless of point size
        )
    )
    if tissue_filter:
        fig.write_html("plots/pca_3d_donors_filter.html")
    else:
        fig.write_html("plots/pca_3d_donors.html")



