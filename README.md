# RNA-seq Multi-dimensional RNA Structure Sequencing Analysis Toolkit
This toolkit is designed for comprehensive analysis and processing of RNA-seq sequencing data with a focus on handling spike-in sequences, calculating various metrics, and correcting GC bias in the data. It consists of several components, each tailored for specific tasks within the RNA-seq data analysis pipeline.

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

## `GC_based_smoothing.py`
Implements a cubic function smoothing model for processing RNA-seq data. This Python script generates a predictive model based on defined parameters for each GC content category, helping in the normalization of sequencing data.

## `GC_based_count_predicting.py`
Predicts unbiased counts for each GC category utilizing a self-benchmarking Gaussian process. This Python script is key to adjusting for GC-content bias in RNA-seq datasets, ensuring more accurate quantification.

## `GC_based_bias_calibrating.py`
Adjusts for transcript position-specific bias to achieve unbiased counts. This script is critical for accurate transcript quantification, correcting for biases introduced during the RNA-seq workflow.

Usage
The tools within this toolkit can be used individually or in combination depending on the specific needs of your RNA-seq data analysis pipeline. Ensure to have the required dependencies installed for each script, including R, Python, and MATLAB environments, as necessary.

## Python Script: `kmeans_cluster_kmer.py`

**Description:**  
This Python script is designed to perform clustering on all k-mers extracted from a specified transcript. It utilizes the unsupervised machine learning algorithm, k-means, with a user-defined number of clusters. This script is particularly useful in genomics for grouping similar sequences, aiding in pattern recognition and data reduction.

**Usage:**  
To use this script, simply provide the sequence data and the desired number of clusters. The output will include the k-mers grouped into the specified number of clusters.

**Key Features:**  
- Implements the k-means clustering algorithm
- Customizable number of clusters
- Optimized for genomic k-mer data

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

License
This project is licensed under the MIT License - see the LICENSE file for details.
