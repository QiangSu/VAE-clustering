
# Read the CSV file into a dataframe
data <- read.csv("C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/umap_GAPDH_201_kmer_splitting_output.csv", header = TRUE, fileEncoding = "UTF-8")
nrow(data)
head(data)
# Sort the dataframe based on the specified column (e.g., "Cluster_traditox")
sorted_data <- data[order(data$cluster_label), ]

# Determine the unique groups in the sorting column
groups <- unique(sorted_data$cluster_label)

# Create a list to store split dataframes
split_data <- vector("list", length(groups))

# Split the data based on the different groups
for (i in seq_along(groups)) {
  split_data[[i]] <- sorted_data[sorted_data$cluster_label == groups[i], ]
}

# Define the maximum number of rows among all groups
max_rows <- max(sapply(split_data, nrow))

# Ensure consistent column names across all dataframes
col_names <- colnames(split_data[[1]])

# Pad the dataframes with NA values and ensure consistent column names
padded_data <- lapply(split_data, function(df) {
  if (nrow(df) < max_rows) {
    missing_rows <- max_rows - nrow(df)
    pad <- data.frame(matrix(NA, nrow = missing_rows, ncol = ncol(df)))
    colnames(pad) <- col_names
    rbind(df, pad)
  } else {
    df
  }
})

# Merge the split dataframes into one dataframe
merged_data <- do.call(cbind, padded_data)

head(merged_data)
# Write the final dataframe to a new CSV file
write.csv(merged_data, "C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/umap_GAPDH_201_kmer_splitting_output_sorted.csv")

# Load the CSV file
df <- merged_data

start_index <- 2
end_index <- 200

column_indices <- c()
for (n in seq(1, end_index)) {
  col1_index <- 4*n - 2
  col2_index <- 4*n - 1
  if(col1_index <= end_index) {
    column_indices <- c(column_indices, col1_index, col2_index)
  }
}
column_indices <- column_indices[column_indices >= start_index]

selected_columns <- df[, column_indices]

# Print the selected columns
print(selected_columns)


write.csv(selected_columns, "C:/new_pipeline/drug_Injury/sorted_tsne-GSE92742_Broad_LINCS_pert_info_10210-kmeans-cluster_21labels-selectedcol.csv")
