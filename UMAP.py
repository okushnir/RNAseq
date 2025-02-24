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
import umap
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

def apply_umap_and_cluster(df, output_dir, output_file, n_clusters, n_neighbors=15, min_dist=0.1):
    """
    Applies UMAP to a DataFrame of genes (rows) and batches (columns), clusters the UMAP results,
    calculates the Silhouette Score, plots the data points colored by their cluster labels,
    and annotates each point with its batch label. Includes the average Silhouette Score in the plot title.

    Parameters:
    - df: pandas.DataFrame, genes as rows and batches as columns.
    - output_dir: str, the directory to save the output plot.
    - output_file: str, the filename for the output plot.
    - n_clusters: int, the number of clusters to find.
    - n_neighbors: int, the number of neighbors considered for each point in UMAP (controls local vs global structure preservation).
    - min_dist: float, the minimum distance between points in the low-dimensional space.

    Returns:
    - umap_results: numpy.ndarray, the two-dimensional UMAP representation of the data.
    - cluster_labels: numpy.ndarray, the cluster labels for each point.
    - silhouette_avg: float, average Silhouette Score for the clustering.
    """
    # Transpose df to get batches as samples and genes as features
    transposed_df = df.transpose()

    # Apply UMAP
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
    umap_results = reducer.fit_transform(transposed_df)

    # Cluster the UMAP results
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = kmeans.fit_predict(umap_results)

    # Calculate Silhouette Score
    silhouette_avg = silhouette_score(umap_results, cluster_labels)

    # Plotting
    plt.figure(figsize=(10, 10))
    cmap = plt.cm.get_cmap('viridis', n_clusters)
    for i in range(umap_results.shape[0]):
        plt.scatter(umap_results[i, 0], umap_results[i, 1], color=cmap(cluster_labels[i]), alpha=0.7)
        plt.text(umap_results[i, 0], umap_results[i, 1], transposed_df.index[i], fontsize=9, ha='right', va='bottom')

    plt.xlabel('UMAP feature 1')
    plt.ylabel('UMAP feature 2')
    plt.title(
        f'UMAP visualization with K-means clustering and batch labels\n'
        f'Average Silhouette Score: {silhouette_avg:.2f}')

    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10, label=f'Cluster {i}')
               for i in range(n_clusters)]
    plt.legend(handles=handles, title="Clusters", loc='upper left')
    plt.savefig(output_dir + output_file)
    plt.close()

    return umap_results, cluster_labels, silhouette_avg


def main():
    sns.set_context("paper")
    # sns.set(font_scale=1.2)
    sns.set_style("ticks")

    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/20240228R_outputs_with_groups/"
    today = datetime.date.today().strftime("%Y%m%d")
    clustering_method = "umap"
    output_file = "{0}_28_samples_clusters.png".format(clustering_method)
    output_dir = input_dir + "{0}_{1}/".format(today, clustering_method)
    n_clusters = 4

    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))

    data = pd.read_csv(input_dir + "final_data_svf_combat_seq.csv")
    data = data.drop(["AD371", "AD371.2", "AD376.2", "AD384.2", "AD387.2", "AD374"], axis=1)
    data = data.set_index("Name")
    inf = np.log10(0)
    data = np.log10(data)
    data = data.loc[(data != inf).all(axis=1), :]
    data = data.astype(float)
    data = data.reset_index()

    genes = data["Name"].tolist()
    data["ensembl_transcript"] = [id.split(".")[0] for id in genes]
    ensembl_df = pd.read_csv(input_dir + "ensembl_hgnc_symbol.csv")
    data = data.merge(ensembl_df, on="ensembl_transcript", how="outer")
    data["hgnc_symbol"] = data["hgnc_symbol"].fillna("None")
    data["hgnc_symbol"] = np.where(data["hgnc_symbol"] == "None", data["Name"], data["hgnc_symbol"])
    data = data.set_index("hgnc_symbol")
    data = data.drop(["Name", "ensembl_transcript"], axis=1)

    umap_results, cluster_labels, silhouette_avg = apply_umap_and_cluster(data, output_dir, output_file, n_clusters=n_clusters)
    data = data.T
    # for feature in data.columns:
    #     for cluster in set(cluster_labels):
    #         print(f"Statistics for {feature} in cluster {cluster}:")
    #         print(data.loc[cluster_labels == cluster, feature].describe())

    # clf = RandomForestClassifier(random_state=42)
    # clf.fit(data, cluster_labels)
    #
    # # Get and plot feature importances
    # importances = clf.feature_importances_
    # indices = np.argsort(importances)[::-1]
    #
    # print("Feature ranking:")
    # for f in range(data.shape[1]):
    #     if importances[indices[f]] > 0:
    #         print(f"{f + 1}. feature {data.columns[indices[f]]} ({importances[indices[f]]})")
    clf = RandomForestClassifier(random_state=42)
    clf.fit(data, cluster_labels)

    # Get feature importances
    feature_importances = clf.feature_importances_

    # Create a DataFrame with features and their importance scores
    features_df = pd.DataFrame({
        'Feature': data.columns,
        'Importance': feature_importances
    })

    # Filter the DataFrame to only include features with importance > 0
    important_features_df = features_df[features_df['Importance'] > 0].copy()

    # Sort the DataFrame by importance in descending order
    important_features_df.sort_values(by='Importance', ascending=False, inplace=True)

    # Reset index for better readability
    important_features_df.reset_index(drop=True, inplace=True)

    important_features_df.to_csv(output_dir + "important_features_df.csv", index=False)

if __name__ == '__main__':
    main()