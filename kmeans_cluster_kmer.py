import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Load the dataset
file_path = 'USF2_201_kmer_clustering_output.csv'
data = pd.read_csv(file_path)

# Assuming the first column is 'index' and the rest are the mer base columns
X = data.iloc[:, 1:51]  # Assuming there are exactly 50 base columns as per the structure

# Optional: Standardize the dataset
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Perform k-means clustering
n_clusters = 1650
kmeans = KMeans(n_clusters=n_clusters, random_state=0)
kmeans.fit(X_scaled)

# Add cluster labels to the original DataFrame
data['cluster_label'] = kmeans.labels_

# Optionally, save the DataFrame with cluster labels to a new CSV file
output_file_path = 'USF2_201_kmer_clustering_1650group_output.csv'
data.to_csv(output_file_path, index=False)

print(f"K-means clustering complete. Results saved to {output_file_path}.")
