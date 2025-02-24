import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import os
import datetime
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
from PCA_plot import pca_with_optional_clustering


def apply_tsne_to_genes(df, output_dir, output_file, perplexity=None):
    """
    Applies t-SNE to a DataFrame of genes (rows) and batches (columns).

    Parameters:
    - df: pandas.DataFrame, genes as rows and batches as columns.

    Returns:
    - tsne_results: numpy.ndarray, the two-dimensional t-SNE representation of the data.
    """
    # Transpose df to get batches as samples and genes as features
    transposed_df = df.transpose()

    # Set perplexity based on the number of samples if not specified
    if perplexity is None:
        perplexity = min(30, transposed_df.shape[0] / 5)

    # Apply t-SNE with adjusted perplexity
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    tsne_results = tsne.fit_transform(transposed_df)

    # Optionally, plot the results
    plt.figure(figsize=(8, 8))
    plt.scatter(tsne_results[:, 0], tsne_results[:, 1])
    plt.xlabel('t-SNE feature 1')
    plt.ylabel('t-SNE feature 2')
    plt.title('t-SNE visualization of gene expression data')
    plt.savefig(output_dir + output_file)
    plt.close()

    return tsne_results


# def apply_tsne_and_cluster(df, output_dir, output_file, n_clusters, perplexity=None):
#     """
#         Applies t-SNE to a DataFrame of genes (rows) and batches (columns), clusters the t-SNE results,
#         calculates the Silhouette Score, plots the data points colored by their cluster labels, and annotates
#         each point with its batch label. Includes the average Silhouette Score in the plot title.
#
#         Parameters:
#         - df: pandas.DataFrame, genes as rows and batches as columns.
#         - n_clusters: int, the number of clusters to find.
#         - perplexity: int or None, the perplexity parameter for t-SNE. If None, sets to min(30, n_samples/5).
#
#         Returns:
#         - tsne_results: numpy.ndarray, the two-dimensional t-SNE representation of the data.
#         - cluster_labels: numpy.ndarray, the cluster labels for each point.
#         - silhouette_avg: float, average Silhouette Score for the clustering.
#         """
#     # Transpose df to get batches as samples and genes as features
#     transposed_df = df.transpose()
#
#     # Set perplexity based on the number of samples if not specified
#     if perplexity is None:
#         perplexity = min(30, transposed_df.shape[0] / 5)
#
#     # Apply t-SNE
#     tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
#     tsne_results = tsne.fit_transform(transposed_df)
#
#     # Cluster the t-SNE results
#     kmeans = KMeans(n_clusters=n_clusters, random_state=42)
#     cluster_labels = kmeans.fit_predict(tsne_results)
#
#     # Calculate Silhouette Score
#     silhouette_avg = silhouette_score(tsne_results, cluster_labels)
#
#     # Create a color map for clusters
#     cmap = plt.cm.get_cmap('viridis', n_clusters)
#
#     # Plot the results with clusters colored and annotated with batch labels
#     plt.figure(figsize=(10, 10))
#     for i in range(tsne_results.shape[0]):
#         plt.scatter(tsne_results[i, 0], tsne_results[i, 1], color=cmap(cluster_labels[i]), alpha=0.7)
#         plt.text(tsne_results[i, 0], tsne_results[i, 1], transposed_df.index[i], fontsize=9, ha='right', va='bottom')
#
#     plt.xlabel('t-SNE feature 1')
#     plt.ylabel('t-SNE feature 2')
#     plt.title(
#         f't-SNE visualization with K-means clustering and batch labels\nAverage Silhouette Score: {silhouette_avg:.2f}')
#
#     # Create a legend for the clusters
#     handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10, label=f'Cluster {i}') for i in range(n_clusters)]
#     plt.legend(handles=handles, title="Clusters", loc='upper left')
#     plt.savefig(output_dir + output_file)
#     plt.close()
#
#     return tsne_results, cluster_labels


def apply_tsne_and_cluster(df, output_dir, output_file, n_clusters, perplexity=None):
    """
    Applies t-SNE to a DataFrame of genes (rows) and batches (columns), clusters the t-SNE results,
    calculates the Silhouette Score and the final Kullback-Leibler divergence from the t-SNE optimization,
    plots the data points colored by their cluster labels, and annotates each point with its batch label.
    Includes the average Silhouette Score and KL divergence in the plot title.

    Parameters:
    - df: pandas.DataFrame, genes as rows and batches as columns.
    - output_dir: str, the directory to save the output plot.
    - output_file: str, the filename for the output plot.
    - n_clusters: int, the number of clusters to find.
    - perplexity: int or None, the perplexity parameter for t-SNE. If None, sets to min(30, n_samples/5).

    Returns:
    - tsne_results: numpy.ndarray, the two-dimensional t-SNE representation of the data.
    - cluster_labels: numpy.ndarray, the cluster labels for each point.
    - silhouette_avg: float, average Silhouette Score for the clustering.
    - kl_divergence: float, the final Kullback-Leibler divergence from the t-SNE optimization.
    """
    # Transpose df to get batches as samples and genes as features
    transposed_df = df.transpose()

    # Set perplexity based on the number of samples if not specified
    if perplexity is None:
        perplexity = min(30, transposed_df.shape[0] / 5)
        print(perplexity)

    # Apply t-SNE
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    tsne_results = tsne.fit_transform(transposed_df)
    kl_divergence = tsne.kl_divergence_

    # Cluster the t-SNE results
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(tsne_results)

    # Calculate Silhouette Score
    silhouette_avg = silhouette_score(tsne_results, cluster_labels)

    # Plotting
    plt.figure(figsize=(10, 10))
    cmap = plt.cm.get_cmap('viridis', n_clusters)
    for i in range(tsne_results.shape[0]):
        plt.scatter(tsne_results[i, 0], tsne_results[i, 1], color=cmap(cluster_labels[i]), alpha=0.7)
        plt.text(tsne_results[i, 0], tsne_results[i, 1], transposed_df.index[i], fontsize=9, ha='right', va='bottom')

    plt.xlabel('t-SNE feature 1')
    plt.ylabel('t-SNE feature 2')
    plt.title(
        f't-SNE visualization with K-means clustering and batch labels\n'
        f'Average Silhouette Score: {silhouette_avg:.2f}, KL Divergence: {kl_divergence:.2f}')

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10, label=f'Cluster {i}')
               for i in range(n_clusters)]
    plt.legend(handles=handles, title="Clusters", loc='upper left')
    plt.savefig(output_dir + output_file)
    plt.close()

    return tsne_results, cluster_labels, silhouette_avg, kl_divergence


def main():
    sns.set_context("paper")
    # sns.set(font_scale=1.2)
    sns.set_style("ticks")

    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/20240228R_outputs_with_groups/"
    today = datetime.date.today().strftime("%Y%m%d")
    clustering_method = "t-sne"
    output_file = "{0}_28_samples_clusters.png".format(clustering_method)
    output_dir = input_dir + "{0}_{1}/".format(today, clustering_method)
    n_clusters = 4

    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))

    """With Filtering"""
    # data = pd.read_csv(groups_dir + "final_data_svf_combat_seq.csv")
    # data = data.rename(columns={"Name": "Gene"})
    # data = data.set_index("Gene")
    # data = data.drop(["AD371", "AD371.2", "AD376.2", "AD384.2", "AD387.2"], axis=1)
    # data["var"] = data.var(axis=1)
    # data["sum"] = data.sum(axis=1)
    # data = data[data["sum"] > len(data.columns)*10]
    # data = data.sort_values(by="var", ascending=False)
    # data = data.drop(["var", "sum"], axis=1)
    """W/O Filtering"""
    data = pd.read_csv(input_dir + "final_data_svf_combat_seq.csv")
    data = data.rename(columns={"Name": "Gene"})
    data = data.set_index("Gene")
    data = data.drop(["AD371", "AD371.2", "AD376.2", "AD384.2", "AD387.2", "AD374"], axis=1)
    inf = np.log10(0)
    data = np.log10(data)
    data = data.loc[(data != inf).all(axis=1), :]
    data = data.astype(float)

    apply_tsne_and_cluster(data, output_dir, output_file, n_clusters=n_clusters)
    pca_with_optional_clustering(data, output_dir, output_file="pca_28_samples_clusters.png", perform_clustering=True,
                                 clustering_method="hierarchical", n_clusters=n_clusters, label_points=True)


if __name__ == '__main__':
    main()