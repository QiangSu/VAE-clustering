import argparse
import time
import multiprocessing
import os
import gzip
from collections import OrderedDict
import glob
from itertools import islice

def count_kmers(sequence, k, filter_set=None):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return kmers if filter_set is None else [kmer for kmer in kmers if kmer in filter_set]

def fastq_sequence_lines(file_path):
    with gzip.open(file_path, 'rt') as fastq_file:
        while True:
            identifier_line = fastq_file.readline()  # Skip the identifier line
            if not identifier_line:
                break  # EOF
            sequence_line = fastq_file.readline().strip()  # Read and strip the sequence line
            plus_line = fastq_file.readline()  # Skip the '+' line
            quality_line = fastq_file.readline()  # Skip the quality line
            yield sequence_line

def process_chunk(sequences, k, filter_set):
    chunk_kmers = []
    for sequence in sequences:
        chunk_kmers.extend(count_kmers(sequence, k, filter_set))
    return chunk_kmers

def main():
    parser = argparse.ArgumentParser(description='Count k-mer frequencies in a FASTQ file.')
    parser.add_argument('--k', type=int, help='Size of the k-mer.', required=True)
    parser.add_argument('--chunk_size', type=int, help='Number of records processed per chunk.', required=True)
    parser.add_argument('--fastq', type=str, help='Path to the FASTQ file.', required=True)
    parser.add_argument('--kmer_dir', type=str, help='Directory containing input CSV files with k-mer sequences.', required=True)
    parser.add_argument('--output', type=str, help='Output directory for storing CSV files.', required=True)
    parser.add_argument('--threads', type=int, help='Number of threads to use for processing.', default=multiprocessing.cpu_count())
    args = parser.parse_args()

    k = args.k
    chunk_size = args.chunk_size
    fastq_file_path = args.fastq
    kmer_dir = args.kmer_dir
    output_directory = args.output
    num_cores = args.threads

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Find all CSV files in the directory
    kmer_files = glob.glob(os.path.join(kmer_dir, "*_kmers.csv"))

    for kmer_csv_file_path in kmer_files:
        kmers_from_csv = OrderedDict()
        with open(kmer_csv_file_path, 'r') as csvfile:
            for line in csvfile:
                kmer = line.split(',')[0].strip()
                if not kmer.lower().startswith("kmer") and kmer:
                    kmers_from_csv[kmer] = 0

        csv_base_name = os.path.splitext(os.path.basename(kmer_csv_file_path))[0]
        output_file_path = os.path.join(output_directory, f"{csv_base_name}_counts.csv")

        start_time = time.time()

        with multiprocessing.Pool(processes=num_cores) as pool:
            chunk_results = []
            for chunk_index in range(0, chunk_size*num_cores, chunk_size):  # iterate over each chunk
                sequences = list(islice(fastq_sequence_lines(fastq_file_path), chunk_index, chunk_index + chunk_size))
                if not sequences:
                    break
                chunk_result = pool.apply_async(process_chunk, args=(sequences, k, set(kmers_from_csv.keys())))
                chunk_results.append(chunk_result)

            for chunk_result in chunk_results:
                chunk_kmers = chunk_result.get()
                for kmer in chunk_kmers:
                    kmers_from_csv[kmer] += 1

        with open(output_file_path, "w") as output_file:
            output_file.write("K-mer,Count\n")
            for kmer, count in kmers_from_csv.items():
                output_file.write(f"{kmer},{count}\n")

        end_time = time.time()
        execution_time = end_time - start_time
        print(f"k-mer count data for {csv_base_name} saved to {output_file_path}")
        print(f"Execution Time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    main()
