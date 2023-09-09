from collections import Counter
from Bio import SeqIO
import gzip
from multiprocessing import Pool
import numpy as np

# Function to convert a k-mer into a quaternary number
def kmer_to_quaternary(kmer):
    quaternary_number = np.zeros(k, dtype=int)
    for i, base in enumerate(kmer):
        if base == "A":
            quaternary_number[i] = 0
        elif base == "T":
            quaternary_number[i] = 1
        elif base == "C":
            quaternary_number[i] = 2
        elif base == "G":
            quaternary_number[i] = 3
    return ''.join(map(str, quaternary_number))

def process_record(record):
    sequence = str(record.seq)
    quaternary_numbers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        quaternary_number = kmer_to_quaternary(kmer)
        quaternary_numbers.append(quaternary_number)
    return quaternary_numbers

file_path1 = "ACTB_NM_001101_1812nt.fasta"
file_path2 = "1032245-N_cutadapt_trim_2P_1000000.fastq.gz"  # Update this to your actual gzipped Fastq file path

k = 50  # k-mer size

# Create a Pool of worker processes for parallel processing
num_processes = 15  # Adjust the number of processes as needed
pool = Pool(num_processes)

# Extract k-mers and convert to quaternary numbers for transcript1 using parallel processing
transcript1_records = list(SeqIO.parse(file_path1, "fasta"))
transcript1_quaternary_numbers = []
batch_size = 100  # Adjust the batch size as needed

for batch_start in range(0, len(transcript1_records), batch_size):
    batch_end = min(batch_start + batch_size, len(transcript1_records))
    batch_records = transcript1_records[batch_start:batch_end]
    batch_quaternary_numbers = pool.map(process_record, batch_records)
    transcript1_quaternary_numbers.extend([item for sublist in batch_quaternary_numbers for item in sublist])

# Extract k-mers and convert to quaternary numbers for transcript2 using parallel processing
transcript2_quaternary_numbers = []
with gzip.open(file_path2, 'rt') as fastq_file:
    transcript2_records = list(SeqIO.parse(fastq_file, "fastq"))

for batch_start in range(0, len(transcript2_records), batch_size):
    batch_end = min(batch_start + batch_size, len(transcript2_records))
    batch_records = transcript2_records[batch_start:batch_end]
    batch_quaternary_numbers = pool.map(process_record, batch_records)
    transcript2_quaternary_numbers.extend([item for sublist in batch_quaternary_numbers for item in sublist])

# Combine all quaternary numbers from transcript1 and transcript2
all_quaternary_numbers = transcript1_quaternary_numbers + transcript2_quaternary_numbers

# Split the combined quaternary numbers into smaller chunks for parallel counting
chunk_size = len(all_quaternary_numbers) // num_processes
chunks = [all_quaternary_numbers[i:i+chunk_size] for i in range(0, len(all_quaternary_numbers), chunk_size)]

# Create a Pool for parallel frequency counting
frequency_counts = pool.map(Counter, chunks)

# Combine frequency counts from all processes
total_frequency_counts = Counter()
for count in frequency_counts:
    total_frequency_counts += count

# Calculate the frequency of each quaternary number in transcript1
transcript1_quaternary_frequency = {}
for quaternary_number in transcript1_quaternary_numbers:
    frequency_in_transcript2 = total_frequency_counts.get(quaternary_number, 0)
    transcript1_quaternary_frequency[quaternary_number] = frequency_in_transcript2

# Save the frequency counts for transcript1 as a text file
output_file_path = "transcript1_frequency_data_1000000_ACTB_multiprocess_batch_chunk.txt"
with open(output_file_path, "w") as output_file:
    for quaternary_number, frequency in transcript1_quaternary_frequency.items():
        output_file.write(f"Quaternary Number: {quaternary_number}, Frequency in Transcript1: {frequency}\n")

print(f"Frequency data saved to {output_file_path}")
