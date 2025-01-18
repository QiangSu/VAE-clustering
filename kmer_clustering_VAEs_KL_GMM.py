#!/usr/bin/env python
import os
import random
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from umap import UMAP  # Import UMAP for dimensionality reduction

from tensorflow import keras
from tensorflow.keras.layers import Input, Dense, Lambda, Reshape, Flatten, Conv1D
from tensorflow.keras.models import Model
import tensorflow.keras.backend as K
import csv

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
# Utility: One-Hot Encode DNA Sequence
# ----------------------------------------------
def one_hot_encode(seq):
    """
    One-hot encode a DNA sequence (A, C, G, T).
    Unknown bases or other characters -> [0, 0, 0, 0].
    """
    mapping = {'A': [1, 0, 0, 0],
               'C': [0, 1, 0, 0],
               'G': [0, 0, 1, 0],
               'T': [0, 0, 0, 1]}
    return np.array([mapping.get(base.upper(), [0, 0, 0, 0]) for base in seq]).flatten()

# ----------------------------------------------
# Utility: Sampling for VAE
# ----------------------------------------------
def sampling(args):
    """
    Reparametrization trick: z = z_mean + exp(0.5 * z_log_var) * eps
    """
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=K.shape(z_mean), seed=42)
    return z_mean + K.exp(0.5 * z_log_var) * epsilon

# ------------------------------------------------------------
# Custom VAE Model
# ------------------------------------------------------------
class VAE(keras.Model):
    def __init__(self, encoder, decoder, input_dim, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.input_dim = input_dim
        # Tracking metrics
        self.total_loss_tracker = keras.metrics.Mean(name='total_loss')
        self.reconstruction_loss_tracker = keras.metrics.Mean(name='reconstruction_loss')
        self.kl_loss_tracker = keras.metrics.Mean(name='kl_loss')

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker
        ]

    def call(self, inputs):
        # Encode
        z_mean, z_log_var, z = self.encoder(inputs)
        # Decode
        reconstruction = self.decoder(z)
        
        # Reconstruction loss
        reconstruction_loss = tf.reduce_mean(
            keras.losses.binary_crossentropy(inputs, reconstruction)
        ) * self.input_dim
        
        # KL loss
        kl_loss = -0.5 * tf.reduce_mean(
            tf.reduce_sum(1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var), axis=1)
        )

        total_loss = reconstruction_loss + kl_loss
        self.add_loss(total_loss)
        
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        
        return reconstruction

# ------------------------------------------------------------
# Main Pipeline
# ------------------------------------------------------------
def run_vae_clustering(
    input_filepath,
    output_filepath,
    cluster_number=None,
    latent_dim=2,      
    epochs=100,        
    batch_size=32,     
    highlight_clusters=None,
    bic_output_path=None
):
    """
    1. Loads data from `input_filepath` (must have 'kmer' column).
    2. Trains VAE on one-hot-encoded sequences to get a latent representation.
    3. If cluster_number is None => scan from k=2..30 (flexible BIC-based).
       Otherwise => use that fixed cluster_number.
    4. Perform final GMM clustering with the chosen k, add cluster labels,
       compute silhouette, etc.
    5. Optionally write BIC file, highlight clusters, produce UMAP plot.
    """
    # 1. Load data with correct delimiter
    # Attempt to read with tab delimiter first
    print(f"[INFO] Attempting to read the input file with tab (',') delimiter.")
    try:
        data = pd.read_csv(input_filepath, delimiter=',', encoding='utf-8')
        print("[INFO] Successfully read the file with tab delimiter.")
    except Exception as e_tab:
        print(f"[WARNING] Failed to read with tab delimiter: {e_tab}")
        # Attempt to read with comma delimiter
        print(f"[INFO] Attempting to read the input file with comma ('\t') delimiter.")
        try:
            data = pd.read_csv(input_filepath, delimiter='\t', encoding='utf-8')
            print("[INFO] Successfully read the file with comma delimiter.")
        except Exception as e_comma:
            print(f"[ERROR] Failed to read with comma delimiter: {e_comma}")
            raise ValueError("Failed to read the input file with both tab and comma delimiters.")

    print(f"[INFO] Columns detected: {data.columns.tolist()}")
    
    # Normalize column names to lowercase and strip whitespace
    data.columns = data.columns.str.strip().str.lower()
    
    if 'kmer' not in data.columns:
        raise ValueError("Input CSV must contain a 'kmer' column.")
    
    kmers = data['kmer'].values
    numeric_kmers = np.array([one_hot_encode(k) for k in kmers], dtype=np.float32)
    input_dim = numeric_kmers.shape[1]
    
    # 2. Split into train/validation
    x_train, x_val = train_test_split(numeric_kmers, test_size=0.2, random_state=42)
    
    # 3. Define and train VAE
    # Encoder
    inputs = Input(shape=(input_dim,), name='encoder_input')
    reshaped = Reshape((input_dim // 4, 4))(inputs)  # e.g., for a 50-mer => (50, 4)

    x = Conv1D(filters=64, kernel_size=3, padding='same', activation='relu')(reshaped)
    x = Conv1D(filters=32, kernel_size=3, padding='same', activation='relu')(x)
    x = Conv1D(filters=16, kernel_size=3, padding='same', activation='relu')(x)
    x = Flatten()(x)

    z_mean = Dense(latent_dim, name='z_mean')(x)
    z_log_var = Dense(latent_dim, name='z_log_var')(x)
    z = Lambda(sampling, name='z')([z_mean, z_log_var])

    encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder')

    # Decoder
    latent_inputs = Input(shape=(latent_dim,), name='z_sampling')
    x = Dense((input_dim // 4) * 16, activation='relu')(latent_inputs)
    x = Reshape((input_dim // 4, 16))(x)
    x = keras.layers.Conv1DTranspose(filters=16, kernel_size=3, padding='same', activation='relu')(x)
    x = keras.layers.Conv1DTranspose(filters=32, kernel_size=3, padding='same', activation='relu')(x)
    x = keras.layers.Conv1DTranspose(filters=64, kernel_size=3, padding='same', activation='relu')(x)
    decoded = keras.layers.Conv1DTranspose(filters=4, kernel_size=3, padding='same', activation='sigmoid')(x)
    outputs = Flatten()(decoded)

    decoder = Model(latent_inputs, outputs, name='decoder')

    vae = VAE(encoder, decoder, input_dim=input_dim)
    vae.compile(optimizer='adam')

    callbacks = [
        keras.callbacks.EarlyStopping(
            monitor='val_loss', patience=5, restore_best_weights=True
        )
    ]
    print("[INFO] Starting VAE training...")
    vae.fit(
        x_train, epochs=epochs, batch_size=batch_size,
        validation_data=(x_val, None),
        callbacks=callbacks,
        verbose=1
    )
    print("[INFO] VAE training completed.")

    # 4. Encode all data => latent space
    print("[INFO] Encoding data into latent space using the trained encoder.")
    z_mean_vals, z_log_var_vals, _ = encoder.predict(numeric_kmers)

    # Reparam (sample) => final latent points for GMM
    rng = np.random.RandomState(42)
    epsilon = rng.randn(*z_mean_vals.shape)
    z_samples = z_mean_vals + np.exp(0.5 * z_log_var_vals) * epsilon

    # 5. Decide cluster approach
    if cluster_number is not None:
        # Fixed cluster mode
        print(f"[INFO] Using fixed cluster number: {cluster_number}")
        best_k = cluster_number
        best_gmm = GaussianMixture(n_components=best_k, covariance_type='full', random_state=42)
        best_gmm.fit(z_samples)
        bic_dict = {}
    else:
        # Flexible cluster mode: scan k=2..30 for best BIC
        print("[INFO] cluster_number not provided => scanning k=2..30 for best BIC")
        bic_dict = {}
        best_bic = float('inf')
        best_k = None
        best_gmm = None

        for k in range(2, 31):
            print(f"[INFO] Fitting GMM with k={k}")
            gmm_temp = GaussianMixture(n_components=k, covariance_type='full', random_state=42)
            gmm_temp.fit(z_samples)
            bic_val = gmm_temp.bic(z_samples)
            bic_dict[k] = bic_val
            print(f"[INFO] BIC for k={k}: {bic_val}")

            if bic_val < best_bic:
                best_bic = bic_val
                best_k = k
                best_gmm = gmm_temp

        print(f"[INFO] Found best_k={best_k} with BIC={best_bic}")

    # 6. Final GMM clustering with best_k
    print(f"[INFO] Performing final GMM clustering with k={best_k}")
    clusters = best_gmm.predict(z_samples)

    # 7. Add cluster labels & compute silhouette
    print("[INFO] Computing silhouette score.")
    data['cluster_label'] = clusters
    sil_score = silhouette_score(z_samples, clusters)
    data['silhouette_score'] = sil_score
    print(f"[INFO] Silhouette Score for k={best_k}: {sil_score}")

    # 8. Write the final CSV
    data.to_csv(output_filepath, index=False)
    print(f"[INFO] Clustering results (k={best_k}) saved to {output_filepath}")

    # 9. Optionally write BIC info if we did scanning (cluster_number=None)
    if bic_output_path and (cluster_number is None):
        with open(bic_output_path, 'w') as f:
            f.write("ClusterNumber\tBIC\n")
            for k in sorted(bic_dict.keys()):
                f.write(f"{k}\t{bic_dict[k]}\n")
            f.write("\n")
            f.write(f"Best cluster number: {best_k}\n")
            f.write(f"Lowest BIC: {bic_dict[best_k]}\n")
        print(f"[INFO] BIC info written to {bic_output_path}")

    # 10. UMAP Visualization
    print("[INFO] Reducing dimensionality with UMAP for visualization.")
    umap_model = UMAP(n_components=2, random_state=42)
    z_umap = umap_model.fit_transform(z_samples)

    plt.figure(figsize=(10, 8))
    # Adjust colormap based on #clusters
    if best_k <= 20:
        cmap = plt.get_cmap('tab20', best_k)
    else:
        cmap = plt.get_cmap('Spectral')

    scatter = plt.scatter(z_umap[:, 0], z_umap[:, 1], c=clusters, cmap=cmap, s=5)

    # 11. Optional highlighting
    if highlight_clusters is not None and len(highlight_clusters) > 0:
        print("[INFO] Highlighting specified clusters.")
        highlight_cmap = plt.get_cmap('tab10', len(highlight_clusters))
        for i, cluster_id in enumerate(highlight_clusters):
            if cluster_id < 0 or cluster_id >= best_k:
                print(f"Warning: cluster ID {cluster_id} is out of range (0..{best_k-1}). Skipping.")
                continue

            idxs = np.where(clusters == cluster_id)[0]
            color = highlight_cmap(i)
            plt.scatter(
                z_umap[idxs, 0], z_umap[idxs, 1],
                color=color, s=50,
                marker='x', linewidths=1.5,
                label=f"Highlighted cluster {cluster_id}"
            )
        plt.legend(loc='best')

    # Colorbar
    cbar = plt.colorbar(scatter, ticks=range(best_k))
    cbar.set_label('Cluster')
    cbar.set_ticks(range(best_k))
    cbar.set_ticklabels(range(best_k))

    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.title(f'UMAP Projection of k-mer Space (k={best_k} Clusters)')
    plt.show()

# ------------------------------------------------------------
# CLI: Argparse
# ------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run VAE + GMM clustering on kmers from a CSV file."
    )
    parser.add_argument(
        "--input", 
        required=True,
        help="Path to the input CSV or TSV file with a 'kmer' column."
    )
    parser.add_argument(
        "--output", 
        required=True,
        help="Path to the output CSV where cluster assignments will be saved."
    )
    parser.add_argument(
        "--clusters", 
        type=int, 
        help=(
            "If provided, GMM uses this fixed cluster number. "
            "If omitted, do a flexible BIC-based approach scanning k=2..30."
        )
    )
    parser.add_argument(
        "--highlight-clusters", 
        type=int, 
        nargs="+", 
        default=None,
        help="List of cluster IDs to highlight, e.g., '--highlight-clusters 0 3'."
    )
    parser.add_argument(
        "--bic-output", 
        default=None,
        help="If provided (and clusters is omitted), write BIC info for k=2..30 to this .txt file."
    )

    args = parser.parse_args()

    # Set seeds for reproducibility
    set_seeds(42)

    run_vae_clustering(
        input_filepath=args.input,
        output_filepath=args.output,
        cluster_number=args.clusters,
        highlight_clusters=args.highlight_clusters,
        bic_output_path=args.bic_output
    )
