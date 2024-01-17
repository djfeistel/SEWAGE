import numpy as np
import random
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

illumina_qscore_dict = {
    0: '!',
    1: '"',
    2: '#',
    3: '$',
    4: '%',
    5: '&',
    6: "'",
    7: '(',
    8: ')',
    9: '*',
    10: '+',
    11: ',',
    12: '-',
    13: '.',
    14: '/',
    15: '0',
    16: '1',
    17: '2',
    18: '3',
    19: '4',
    20: '5',
    21: '6',
    22: '7',
    23: '8',
    24: '9',
    25: ':',
    26: ';',
    27: '<',
    28: '=',
    29: '>',
    30: '?',
    31: '@',
    32: 'A',
    33: 'B',
    34: 'C',
    35: 'D',
    36: 'E',
    37: 'F',
    38: 'G',
    39: 'H',
    40: 'I'
}

def main_2():
    # Parameters
    sequence_length = 5000000
    fragement_size = 200
    read_size = 150
    total_depth_coverage = 100
    random_seed = None

    fasta_dict = load_fasta_into_dict(sys.argv[1])
    cov_per_genome = genome_proporitons(fasta_dict, total_depth_coverage)

    for amplicon, full_sequence in fasta_dict.items():
        genome_name = amplicon.split("~~~")[1]
        if genome_name in cov_per_genome:
            genome_cov = cov_per_genome[genome_name]
        seq_len = len(full_sequence)
        ## list of indicies based on parameters above
        fragemnt_index_list = fragment_indices_list(
            seq_len=seq_len,
            frag_len=fragement_size,
            read_len=read_size,
            depth=genome_cov,
            seed=random_seed
        )


        for index in fragemnt_index_list:
            start, end = index
            fragemnt = full_sequence[start:end]
            f, r = fragemnt[:read_size], fragemnt[-read_size-1:-1]
            print(f"name: {amplicon}\nforward: {f}\nreverse: {r}\n")


def name_fastq_read(
        amplicon_defline: str, 
        sequence:str, 
        read_count:int, 
        reverse: bool=False
        ) -> list:
    '''
    add the name R1 or R2 to the defline of the read
        quality scores are (as of now) ignored and are given a Q40

    :param amplicon_defline: defline used in amplicon for tracking 
    :param sequence: Forward/R1 or reverse/R2 sequence
    :param read_count: Integer used in paired end reads
    :param reverse: set True if read is R2, else False is R1
    :return:
    '''
    if reverse:
        read_type = "R2"
    else:
        read_type = "R1"

    defline = f"@{amplicon_defline}:{read_type}:{read_count}"
    
    return [defline, sequence, "+", "I"*len(sequence)]

def genome_fasta_dictionary(fasta_file):
    """
    Reads a FASTA file and loads it into a dictionary.

    :param fasta_file: Path to the FASTA file.
    :return: A dictionary with headers as keys and sequences as values.
    """
    with open(fasta_file, 'r') as file:
        fasta_dict = {}
        header = None

        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                header = line[1:]  # Remove the '>' character
                if header in fasta_dict:
                    raise ValueError(f"Duplicate header found: {header}")
                fasta_dict[header] = ''
            else:
                if header is None:
                    raise ValueError("File format error: Sequence without a header")
                fasta_dict[header] += line

        return fasta_dict

def genome_proporitons(fasta_dict, coverage):
    genomes_keys = fasta_dict.keys()
    genomes_set = sorted(set([i.split('~~~')[1] for i in genomes_keys]))
    random_proportions = np.random.dirichlet(np.ones(len(genomes_set)), size=1)[0]
    return {i: j for i, j in zip(genomes_set, random_proportions*coverage)}

def get_truncated_normal(mean=0, sd=4, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def skewed_normal_distribution(size, mean:int=42, lower_limit:int=0, upper_limit:int=40, sd:int=10):
    # Create a truncated normal distribution
    tr_norm = get_truncated_normal(mean=mean, sd=sd, low=lower_limit, upp=upper_limit)

    # Generate values and round to get discrete numbers
    values = tr_norm.rvs(size)
    discrete_values = np.round(values)

    return discrete_values

def complement_dna(sequence):
    sequence = sequence
    compliment = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return ''.join(compliment[nucleotide] for nucleotide in sequence)

def fragment_indices_list(
        seq_len: int,
        frag_len: int=250,
        read_len: int=150,
        depth: int=30,
        seed=None
) -> list:
    '''
    creates a sorted list of random numbers which are used as the index for slicing an amplicon
        numbers are randomly picked from a uniform distribution
        upper limit is the length of the sequence minus 1 (bc of pythonic index) and minus the length of the fragment
        reverse list is frag_len + random number
    :param seq_len: length of sequence
    :param frag_len: length of fragment
    :param read_len: length of F/R read
    :param depth: depth of coverage
    :param seed: reproducibility
    :return: list of list with sorted values
    '''
    np.random.seed(seed)
    total_reads = int(depth * seq_len / read_len)
    upper_limit = (seq_len + 1) - frag_len

    forward_read_indicies = np.sort(np.random.randint(0, upper_limit, size=total_reads))
    reverse_read_indicies = forward_read_indicies + frag_len
    index_list = [[forward_read_indicies[x], reverse_read_indicies[x]] for x in range(len(forward_read_indicies))]
    return index_list


########
def main():
    num_reads = 10000  # Number of reads to generate
    read_length = 150  # Length of each read
    start_mean_quality = 40  # High quality at the start
    end_mean_quality = 30  # Lower quality towards the end
    std_dev_quality = 2  # Standard deviation for quality scores
    std_rate = 0.08
    # Generating reads
    # Generating paired-end reads
    paired_reads = [create_paired_end_reads(
        read_length,
        f"read_{i + 1}",
        start_mean=start_mean_quality,
        end_mean=end_mean_quality,
        std_dev=std_dev_quality,
        std_rate=std_rate
    )
        for i in range(num_reads)
    ]

    # Writing to separate FASTQ files for forward and reverse reads
    with open("forward_reads.fastq", "w") as fwd_output, open("reverse_reads.fastq", "w") as rev_output:
        for fwd_read, rev_read in paired_reads:
            SeqIO.write([fwd_read], fwd_output, "fastq")
            SeqIO.write([rev_read], rev_output, "fastq")

    return None

def random_dna_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def quality_scores_normal_distribution(
        length: int,
        start_mean: int,
        end_mean: int,
        std_dev: int
) -> list:
    '''
    quality
    :param length: int for length of sequence
    :param start_mean: start quality score between 40-0
    :param end_mean: end quality score between 40-0
    :param std_dev: STD of quality score
    :return: list of quality scores that are using for each read
    '''
    linear_scale = np.linspace(start_mean, end_mean, length)
    quality_scores = [int(np.random.normal(mean, std_dev)) for mean in linear_scale]
    # Clip the values to be within valid quality score range (e.g., 0 to 40)
    quality_scores = [min(max(0, score), 40) for score in quality_scores]
    return quality_scores

def quality_scores_normal_distribution_poor_quality_data(
        length=150,
        start_mean=40,
        end_mean=30,
        start_std_dev=2,
        std_dev_increase_rate=0.08
) -> list:
    '''
    create fastq reads that increase in SD towards the tail region
    :param length: int for length of sequence
    :param start_mean: start quality score between 40-0
    :param end_mean: end quality score between 40-0
    :param start_std_dev: STD of quality score start
    :param std_dev_increase_rate: rate at which to increase the SD
    :return:
    '''
    linear_scale = np.linspace(start_mean, end_mean, length)
    std_dev_scale = np.linspace(start_std_dev, start_std_dev + std_dev_increase_rate * length, length)
    quality_scores = [int(np.random.normal(mean, std_dev)) for mean, std_dev in zip(linear_scale, std_dev_scale)]
    # Clip the values to be within valid quality score range (e.g., 0 to 40)
    quality_scores = [min(max(0, score), 40) for score in quality_scores]
    return quality_scores

def create_paired_end_reads(seq_length, read_id, start_mean, end_mean, std_dev, std_rate):
    # Forward read
    fwd_sequence = random_dna_sequence(seq_length)
    #fwd_quality_scores = quality_scores_normal_distribution(seq_length, start_mean, end_mean, std_dev)
    fwd_quality_scores = quality_scores_normal_distribution_poor_quality_data(seq_length, start_mean, end_mean, std_dev, std_rate)
    fwd_record = SeqRecord(
        Seq(fwd_sequence),
        id=read_id,
        description="Forward read",
        letter_annotations={"phred_quality": fwd_quality_scores}
    )

    # Reverse read (reverse complement)
    start_mean_rev = start_mean - 3
    end_mean_rev = end_mean - 3
    rev_sequence = fwd_record.seq.reverse_complement()
    # Optionally, create a different quality score profile for the reverse read
    #rev_quality_scores = quality_scores_normal_distribution(seq_length, start_mean_rev, end_mean_rev, std_dev)
    rev_quality_scores = quality_scores_normal_distribution_poor_quality_data(seq_length, start_mean_rev, end_mean_rev, std_dev, std_rate)
    rev_record = SeqRecord(
        rev_sequence,
        id=read_id,
        description="Reverse read",
        letter_annotations={"phred_quality": rev_quality_scores}
    )

    return fwd_record, rev_record


if __name__ == "__main__":
    main_2()


