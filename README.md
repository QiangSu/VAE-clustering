# RNA-seq Analysis Toolkit
This toolkit is designed for comprehensive analysis and processing of RNA-seq sequencing data with a focus on handling spike-in sequences, calculating various metrics, and correcting GC bias in the data. It consists of several components, each tailored for specific tasks within the RNA-seq data analysis pipeline.

Components
model_sequence-N8-0708.R
This script generates RNA templates with varying N sequences for RNA-seq sequencing data analysis. For each template, counts are computed based on all incorporated spike-in sequences. This tool is instrumental in preparing and validating the spike-in controls.

## pipeline.sh
A shell script optimized for cleaning the RNA-seq data by removing contaminants from both spike-in sequences and natural samples. It serves as a preliminary step in ensuring data quality before in-depth analysis.

## quaternary_count_multiprocess_batch_chunk.py
This Python script is designed for efficient calculation of k-mer counts from spike-in RNA sequences or natural transcripts. By leveraging multiprocessing, it can handle large datasets in a batch mode, significantly reducing computation time.

## RNA-extract-fragment.R
This R script is focused on extracting k-mers from natural transcripts, facilitating the analysis of sequence motifs and their characteristics in a given RNA-seq dataset.

## MFE_cal_k-mer.m
A MATLAB script for calculating the Minimum Free Energy (MFE) values of k-mers. This is crucial for understanding the stability and structural propensity of RNA sequences under study.

## GC_based_smoothing.py
Implements a cubic function smoothing model for processing RNA-seq data. This Python script generates a predictive model based on defined parameters for each GC content category, helping in the normalization of sequencing data.

## GC_based_count_predicting.py
Predicts unbiased counts for each GC category utilizing a self-benchmarking Gaussian process. This Python script is key to adjusting for GC-content bias in RNA-seq datasets, ensuring more accurate quantification.

## GC_based_bias_calibrating.py
Adjusts for transcript position-specific bias to achieve unbiased counts. This script is critical for accurate transcript quantification, correcting for biases introduced during the RNA-seq workflow.

Usage
The tools within this toolkit can be used individually or in combination depending on the specific needs of your RNA-seq data analysis pipeline. Ensure to have the required dependencies installed for each script, including R, Python, and MATLAB environments, as necessary.



Contributing
Contributions to this toolkit are welcome. Please open an issue to discuss proposed changes or submit a pull request.

License
This project is licensed under the MIT License - see the LICENSE file for details.
