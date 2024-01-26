#!/usr/bin/env python3

import sys
import subprocess
from Bio import SeqIO, AlignIO

def align_genomes(reference_genome_path, forward_reads_path, reverse_reads_path, output_file, threads:int=1):
    threads = str(threads)
    # Align paired-end reads to the reference genome using minimap2
    subprocess.run([
        "minimap2",
        "-ax", "sr",  # 'sr' is for short reads; change as per your data type
        "-t", threads,
        reference_genome_path,
        forward_reads_path,
        reverse_reads_path,
        "-o", output_file
    ])

def find_conserved_regions(alignment_file, threshold=0.9):
    alignment = AlignIO.read(alignment_file, "sam")
    conserved_regions = []

    for i in range(len(alignment[0])):  # Iterate over column positions
        column = alignment[:, i]  # Get a single column
        most_common_nuc = max(set(column), key = column.count)
        freq = column.count(most_common_nuc) / len(column)

        if freq >= threshold:  # Check if the most common nucleotide is above the threshold
            conserved_regions.append((i, most_common_nuc, freq))

    return conserved_regions

def main():
    reference_genome_path = "/Users/djfeistel/GitHub/SEWAGE/db/SARS-CoV-2.reference.fasta"
    forward_reads_path = sys.argv[1]
    reverse_reads_path = sys.argv[2]
    alignment_output = sys.argv[3]

    # Align the genomes
    align_genomes(reference_genome_path, forward_reads_path, reverse_reads_path, alignment_output)

    # # Identify conserved regions from the alignment
    conserved = find_conserved_regions(sys.argv[3])
    for position, nucleotide, frequency in conserved:
        print(f"Position: {position}, Nucleotide: {nucleotide}, Frequency: {frequency}")
    # # Further processing for primer design goes here...

if __name__ == "__main__":
    main()
