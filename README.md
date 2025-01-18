# RNA-seq Multi-dimensional RNA Structure Sequencing Analysis Toolkit
This repository provides a VAE–GMM–based pipeline for exploring how RNA structure impacts reverse transcription (RT) efficiency in RNA-seq experiments, complemented by GC content-based and MFE-based Gaussian distribution models for additional structure-related analyses. By combining Variational Autoencoders (VAEs) with Gaussian Mixture Models (GMMs), the pipeline effectively captures high-dimensional k-mer patterns reflecting both RNA secondary and tertiary structures, revealing how they influence sequencing biases during RT. 

Usage
The tools within this toolkit can be used individually or in combination depending on the specific needs of your RNA-seq data analysis pipeline. Ensure to have the required dependencies installed for each script, including R, Python, and MATLAB environments, as necessary.

## Python Script: `kmer_clustering_VAEs_KL_GMM.py`

**Description:**  
This Python script performs Variational Autoencoder (VAE)–based dimensionality reduction followed by Gaussian Mixture Model (GMM) clustering on k-mer sequences. By first learning a low-dimensional latent representation of one-hot–encoded k-mers, the script can either use a user-specified number of clusters or automatically determine an optimal number of clusters by scanning for the best Bayesian Information Criterion (BIC). It also offers UMAP visualization of clustered k-mers and supports highlighting specific clusters of interest.

**Usage:**  
running VAE-GMM script
```python
python kmer_clustering_VAEs_KL_GMM.py \
    --input <input_file> \
    --output <output_file> \
    [--clusters <cluster_number>] \
    [--highlight-clusters <list_of_clusters>] \
    [--bic-output <bic_info_file>]
```
- --input: Path to the input CSV or TSV file containing a kmer column.
- --output: Name (and path) of the output CSV where cluster assignments and scores will be saved.
- --clusters (optional): If provided, the GMM will use exactly that many components for clustering; if omitted, the script will scan k=2..30 and pick the best k by BIC.
- --highlight-clusters (optional): Space-separated cluster IDs to highlight in the UMAP plot (e.g., 0 3 5).
- --bic-output (optional): If provided (and no fixed --clusters specified), writes per-cluster BIC scores and the chosen best cluster size to a text file.

**Key Features:**  
- VAE Dimensionality Reduction: Learns compressed representations of k-mers, capturing key structural or sequence patterns.
- Flexible or Fixed Clustering: Either specify --clusters directly or let the script automatically find the best number of clusters via BIC.
- UMAP Visualization: Generates a 2D projection of the latent space, colored by cluster labels, with optional highlighting of specific clusters.
- Silhouette Score: Outputs a silhouette score to gauge the quality of clustering.
- BIC Logging: If scanning for clusters, logs the BIC values for each tested cluster number and notes the optimal choice.


## Python Script: `kmer_clustering_GMM_PCA.py`

**Description:**  
This Python script clusters DNA k-mers using Principal Component Analysis (PCA) for dimensionality reduction and a Gaussian Mixture Model (GMM) for clustering. Optionally, it visualizes the clusters in UMAP space, computes the silhouette score, and highlights specified cluster IDs.

**Usage:**  
running GMM-PCA script
```python
python kmer_clustering_GMM_PCA.py \
    --input <input_csv> \
    --output <output_csv> \
    [--clusters <n_clusters>] \
    [--pca_dim <n_components>] \
    [--highlight-clusters <list_of_clusters>]
```
- --input: Path to the input CSV or TSV file containing a kmer column.
- --output: Name (and path) of the output CSV where cluster assignments and scores will be saved.
- --clusters (optional): If provided, the GMM will use exactly that many components for clustering; if omitted, the script will scan k=2..30 and pick the best k by BIC.
- --pca_dim (optional): Number of PCA components (default = 2).
- --highlight-clusters (optional): Space-separated list of cluster IDs to highlight on the UMAP plot (e.g., 0 1 5).

**Key Features:**  
- PCA-Based Dimensionality Reduction
- GMM Clustering
- Silhouette Score
- UMAP Visualization
 

## Python Script: `kmer_clustering_Kmeans_PCA.py`

**Description:**  
This Python script clusters DNA k-mers using Principal Component Analysis (PCA) for dimensionality reduction, followed by K-Means clustering. It then computes the silhouette score to measure clustering quality and generates an optional UMAP visualization of k-mer distributions, with the ability to highlight selected clusters.

**Usage:**  
running Kmeans-PCA script
```python
python kmer_clustering_Kmeans_PCA.py \
    --input <input_csv> \
    --output <output_csv> \
    [--clusters <n_clusters>] \
    [--pca_dim <n_components>] \
    [--highlight-clusters <list_of_clusters>]
```
- --input: Path to the input CSV or TSV file containing a kmer column.
- --output: Name (and path) of the output CSV where cluster assignments and scores will be saved.
- --clusters (optional): If provided, the GMM will use exactly that many components for clustering; if omitted, the script will scan k=2..30 and pick the best k by BIC.
- --pca_dim (optional): Number of PCA components (default = 2).
- --highlight-clusters (optional): Space-separated list of cluster IDs to highlight on the UMAP plot (e.g., 0 1 5).

**Key Features:**  
- PCA-Based Dimensionality Reduction
- Kmeans Clustering
- Silhouette Score
- UMAP Visualization


## Python Script: `kmer_clustering_Hierarchical_PCA.py`

**Description:**  
This Python script clusters DNA k-mers using Principal Component Analysis (PCA) for dimensionality reduction, followed by Hierarchical clustering. It then computes the silhouette score to measure clustering quality and generates an optional UMAP visualization of k-mer distributions, with the ability to highlight selected clusters.

**Usage:**  
running Hierarchical-PCA script
```python
python kmer_clustering_Hierarchical_PCA.py \
    --input <input_csv> \
    --output <output_csv> \
    [--clusters <n_clusters>] \
    [--pca_dim <n_components>] \
    [--highlight-clusters <list_of_clusters>]
```
- --input: Path to the input CSV or TSV file containing a kmer column.
- --output: Name (and path) of the output CSV where cluster assignments and scores will be saved.
- --clusters (optional): If provided, the GMM will use exactly that many components for clustering; if omitted, the script will scan k=2..30 and pick the best k by BIC.
- --pca_dim (optional): Number of PCA components (default = 2).
- --highlight-clusters (optional): Space-separated list of cluster IDs to highlight on the UMAP plot (e.g., 0 1 5).

**Key Features:**  
- PCA-Based Dimensionality Reduction
- Hierarchical Clustering
- Silhouette Score
- UMAP Visualization 


## Python Script: `kmer_clustering_GMM.py`

**Description:**  
This Python script applies Gaussian Mixture Model (GMM) clustering to one-hot–encoded DNA k-mers, computes a silhouette score to evaluate clustering quality, and creates a UMAP visualization for intuitive inspection. It also supports highlighting selected clusters of interest in the visualization.

**Usage:**  
running GMM script
```python
python kmer_clustering_GMM.py \
    --input <input_csv> \
    --output <output_csv> \
    [--clusters <n_clusters>] \
    [--highlight-clusters <list_of_clusters>]
```
- --input: Path to the input CSV or TSV file containing a kmer column.
- --output: Name (and path) of the output CSV where cluster assignments and scores will be saved.
- --clusters (optional): If provided, the GMM will use exactly that many components for clustering; if omitted, the script will scan k=2..30 and pick the best k by BIC.
- --highlight-clusters (optional): Space-separated list of cluster IDs to highlight on the UMAP plot (e.g., 0 1 5).

**Key Features:**  
- Hierarchical Clustering
- Silhouette Score
- UMAP Visualization

  
## R Script: `UMAP_Splitted_kmers.R`

**Description:**  
`UMAP_Splitted_kmers.R` is an R script geared towards dimensionality reduction of k-mer data from specific transcripts. It employs the Uniform Manifold Approximation and Projection (UMAP) technique to reduce the high-dimensional space into two main dimensions, UMAP1 and UMAP2. These dimensions are further used to visualize the distribution and relationships among k-mers, facilitating easier interpretation and analysis.

**Usage:**  
Input your dataset containing k-mers and the script will handle the dimensionality reduction, culminating in a plot displaying the k-mers in the two-dimensional UMAP space.

**Key Features:**  
- Utilizes UMAP for efficient dimensionality reduction
- Outputs a 2D visualization of k-mers
- Easy interpretation of complex k-mer relationships

## R Script: `sorting_cluster.R`

**Description:**  
The `sorting_cluster.R` script is crafted to organize and categorize k-mers based on their respective cluster assignment, specifically targeting 50-mers. By sorting the k-mers according to the cluster names, this script enhances data management and visualization readiness, making it an essential tool for researchers analyzing clustered sequence data.

**Usage:**  
Feed the script with clustered k-mer data, and it will sort these sequences by cluster, preparing them for subsequent graphical representation or further analysis.

**Key Features:**  
- Efficient sorting of k-mers by cluster names
- Specifically handles 50-mer sequences
- Prepares data for visualization and further analysis

## `Shannon_entropy_cal.py`
Description:
Shannon_entropy_cal.py is a Python script tailored to calculate Shannon's entropy for each k-mer extracted from a given transcript. Shannon's entropy, a measure of the randomness or uncertainty inherent in a data set, is crucial for analyzing the complexity and diversity of k-mers in genetic sequences. This script provides a valuable tool for bioinformaticians and geneticists looking to assess the informational content of sequences within genomic data.

Usage:
To utilize this script, provide it with the k-mer data from your transcript. The script will compute and return the Shannon entropy for each k-mer, allowing you to evaluate the variability within the sequence data.

Key Features:
Calculates Shannon entropy for detailed sequence analysis
Supports analysis of any k-mer length
Efficient and accurate entropy computation

Components
## `model_sequence-N8-0708.R`
This script generates RNA templates with varying N sequences for RNA-seq sequencing data analysis. For each template, counts are computed based on all incorporated spike-in sequences. This tool is instrumental in preparing and validating the spike-in controls.

## `pipeline.sh`
A shell script optimized for cleaning the RNA-seq data by removing contaminants from both spike-in sequences and natural samples. It serves as a preliminary step in ensuring data quality before in-depth analysis.

## `kmer_counting_loop.py`
This Python script is designed for efficient calculation of k-mer counts from spike-in RNA sequences or natural transcripts. By leveraging multiprocessing, it can handle large datasets in a batch mode, significantly reducing computation time.

## `RNA-extract-fragment.R`
This R script is focused on extracting k-mers from natural transcripts, facilitating the analysis of sequence motifs and their characteristics in a given RNA-seq dataset.

## `MFE_cal_k-mer.m`
A MATLAB script for calculating the Minimum Free Energy (MFE) values of k-mers. This is crucial for understanding the stability and structural propensity of RNA sequences under study.



License
This project is licensed under the MIT License - see the LICENSE file for details.
