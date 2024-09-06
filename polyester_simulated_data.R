# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install necessary packages safely, checking for library write permissions
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  if (file.access(.libPaths()[1], 2) == 0) {
    install.packages("BiocManager")
  } else {
    warning("Default library path is not writable. Try setting a personal library path or run R as an administrator.")
  }
}

# Install Bioconductor packages with forced re-installation if needed
BiocManager::install(c("polyester", "Biostrings"), version = "3.18", force = TRUE)

# Install other necessary packages
install.packages("profvis")
install.packages("htmlwidgets")

# Load necessary libraries
library(BiocManager)
library(polyester)
library(Biostrings)
library(profvis)
library(htmlwidgets)

# Define file paths
reference_transcriptome <- "/home/data/qs/data/reference_isoform/Homo_sapiens.GRCh38.cdna.all.fa"
abundance_file <- "/home/data/qs/data/reference_isoform/99272N_results_file_GC.csv"

# Check if the reference transcriptome exists
if (!file.exists(reference_transcriptome)) {
  stop("Reference transcriptome file does not exist.")
}

# Load the reference transcriptome
fasta_info <- readDNAStringSet(reference_transcriptome)
all_transcript_ids <- names(fasta_info)

# Modify transcript IDs to remove version numbers
all_transcript_ids_base <- sapply(all_transcript_ids, function(id) sub("\\..*", "", id))

# Load abundance data
if (!file.exists(abundance_file)) {
  stop("Abundance file does not exist.")
}
abundance_data <- read.csv(abundance_file)

# Match abundance data with transcript IDs
matched_abundance <- abundance_data$Abundance[match(all_transcript_ids_base, abundance_data$Transcript_ID, nomatch = NA_integer_)]
matched_abundance[is.na(matched_abundance)] <- 0  # Handle unmatched transcripts by setting their abundance to 0

# Create the abundance matrix
abundance_matrix <- matrix(matched_abundance, ncol = 1)
rownames(abundance_matrix) <- all_transcript_ids_base
colnames(abundance_matrix) <- "sample"

# Set up the output directory
output_dir <- "simulated_data_99272N_828isoform"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define fold changes (if needed, otherwise using 1 for all)
fold_changes <- rep(1, length(all_transcript_ids_base))

# Check for Pandoc installation
if (Sys.which("pandoc") == "") {
  stop("Pandoc is not installed. Please install Pandoc from https://pandoc.org/installing.html and ensure it's in your PATH.")
}

# Simulation function wrapped with performance profiling
simulate_with_profiling <- function() {
  simulate_experiment(
    fasta = reference_transcriptome,
    outdir = output_dir,
    num_reps = 1,
    abundance = abundance_matrix,
    fold_changes = fold_changes,
    paired = TRUE,
    readlen = 150,
    numreads = 10000000  # Specify the total number of reads to generate
  )
}

# Profile the simulation and save the results
prof <- profvis({
  simulate_with_profiling()
})

html_file <- file.path(output_dir, "profiling_report.html")
saveWidget(prof, html_file)

cat("Profiling report saved to: ", html_file, "\n")

# Track runtime and memory usage
start_time <- Sys.time()
simulate_with_profiling()
end_time <- Sys.time()

memory_usage <- gc()
runtime <- end_time - start_time

tracking_info <- paste("Runtime: ", runtime, "\nMemory Usage: \n", capture.output(print(memory_usage)), sep = "\n")
writeLines(tracking_info, file.path(output_dir, "simulation_tracking.txt"))
cat("Simulation tracking info saved to: ", file.path(output_dir, "simulation_tracking.txt"), "\n")
