#!/usr/bin/env python
import os
import random
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA  # <-- Import PCA
from sklearn.cluster import AgglomerativeClustering  # <-- Import AgglomerativeClustering for hierarchical clustering
from sklearn.metrics import silhouette_score
from umap import UMAP  # <-- Import UMAP for dimensionality reduction

# ----------------------------------------------
# 0. Set Random Seeds for Reproducibility
# ----------------------------------------------
def set_seeds(seed=42):
    """
    Sets the random seeds for Python, NumPy, and TensorFlow 
    to ensure reproducible behavior.
    """
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)

# ----------------------------------------------
# One-hot encoding for DNA k-mers
# ----------------------------------------------
def one_hot_encode(seq):
    """
    Convert a DNA sequence to a 1D NumPy array using one-hot encoding for A, C, G, T.
    Unknown or unexpected characters are mapped to [0,0,0,0].
    """
    mapping = {
        'A': [1, 0, 0, 0],
        'C': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'T': [0, 0, 0, 1]
    }
    return np.array([mapping.get(base, [0, 0, 0, 0]) for base in seq]).flatten()

def main(args):
    # 1. Set seeds for reproducibility
    set_seeds(42)

    # 2. Load the data
    data = pd.read_csv(args.input)
    if 'kmer' not in data.columns:
        raise ValueError("Input CSV must contain a column named 'kmer'.")

    # 3. Extract k-mer sequences and one-hot encode
    kmers = data['kmer'].values
    numeric_kmers = np.array([one_hot_encode(kmer) for kmer in kmers])

    # 4. Dimensionality reduction with PCA (now uses args.pca_dim)
    pca = PCA(n_components=args.pca_dim, random_state=42)
    reduced_kmers = pca.fit_transform(numeric_kmers)
    print(f"Explained variance by PCA ({args.pca_dim} components): {pca.explained_variance_ratio_}")

    # 5. Apply Hierarchical Clustering on PCA-reduced data
    hierarchical = AgglomerativeClustering(n_clusters=args.clusters, linkage='ward')
    clusters = hierarchical.fit_predict(reduced_kmers)

    # 6. Add cluster labels to the DataFrame
    data['Cluster'] = clusters

    # 7. Compute the silhouette score
    sil_score = silhouette_score(reduced_kmers, clusters)
    print(f"Silhouette Score: {sil_score}")

    # 8. Add the same silhouette score in a new column for each row
    data['silhouette_score'] = sil_score

    # 9. (Optional) Reorder columns so that silhouette_score appears right after 'Cluster'
    cols = list(data.columns)
    cluster_index = cols.index('Cluster') + 1
    cols.remove('silhouette_score')
    cols.insert(cluster_index, 'silhouette_score')
    data = data[cols]

    # 10. Save the DataFrame to a new CSV file
    data.to_csv(args.output, index=False)
    print(f"Clustering complete. Results (including silhouette scores) saved to {args.output}")

    # 11. Reduce dimensionality using UMAP for visualization only (still 2D)
    umap_model = UMAP(n_components=2, random_state=42)
    z_umap = umap_model.fit_transform(reduced_kmers)

    # 12. Plot clusters with a dynamic color map
    cluster_number = args.clusters
    plt.figure(figsize=(10, 8))

    if cluster_number <= 20:
        cmap = plt.get_cmap('tab20', cluster_number)  # Get first 'cluster_number' colors
    else:
        cmap = plt.get_cmap('Spectral')

    scatter = plt.scatter(z_umap[:, 0], z_umap[:, 1], c=clusters, cmap=cmap, s=5)

    # Create colorbar with correct ticks and labels
    cbar = plt.colorbar(scatter, ticks=range(cluster_number))
    cbar.set_label('Cluster')
    cbar.set_ticks(range(cluster_number))
    cbar.set_ticklabels(range(cluster_number))

    # 13. Optional cluster highlighting
    if args.highlight_clusters is not None and len(args.highlight_clusters) > 0:
        highlight_cmap = plt.get_cmap('tab10', len(args.highlight_clusters))
        for i, cid in enumerate(args.highlight_clusters):
            if cid < 0 or cid >= cluster_number:
                print(f"Warning: cluster ID {cid} is out of range (0..{cluster_number-1}). Skipping.")
                continue

            idxs = np.where(np.array(clusters) == cid)[0]
            color = highlight_cmap(i)
            plt.scatter(
                z_umap[idxs, 0], z_umap[idxs, 1],
                color=color, s=50,
                marker='x', linewidths=1.5,
                label=f"Highlighted cluster {cid}"
            )
        plt.legend(loc='best')

    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.title(f'UMAP Projection of k-mer Space ({cluster_number} Clusters)')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster k-mers from a CSV file using PCA, Hierarchical Clustering, compute silhouette score, and visualize via UMAP."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input CSV file containing a 'kmer' column."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output CSV file where cluster assignments will be saved."
    )
    parser.add_argument(
        "--clusters",
        type=int,
        default=10,
        help="Number of hierarchical clusters to form (default: 10)."
    )
    parser.add_argument(
        "--pca_dim",
        type=int,
        default=2,
        help="Number of principal components to use for PCA (default: 2)."
    )
    parser.add_argument(
        "--highlight-clusters",
        nargs='*',
        type=int,
        default=None,
        help="List of cluster IDs to highlight (e.g., 0 2)."
    )
    args = parser.parse_args()
    
    main(args)
