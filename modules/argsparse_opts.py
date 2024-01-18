import argparse
import os

def sewage_opts():
    parser = argparse.ArgumentParser(
        prog="SEWAGE",
        description=f"\nSynthetically Engineered Wastewater-like sequence data for Assessing Genomic and Environmental populations\n",
        epilog=f"Fast Start: SEWAGE -i <multi.fasta> -s <scheme>"
    )
    '''input and output'''
    required_pathway = parser.add_argument_group("Input and Scheme Parameters", "Required flags")
    required_pathway.add_argument(
        "-i", "--infasta", 
        help="Multifasta file or single column list of pathways to fasta files", 
        required=True,
        dest='infasta',
        type=str,
        default=None,
        metavar='STR'
    )
    required_pathway.add_argument(
        '-s', '--scheme', 
        help=f'Available primer scheme: \nArtic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"]\n \
        VarSkip: Long read = ["vsl1a"]; Short-read = ["vss1a", "vss2a", "vss2b"])', 
        required=True, 
        type=str,
        choices=["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vsl1a", "vss1a", "vss2a", "vss2b"],
        metavar='STR',
        dest='scheme',
        default=None
    )
    amplicon_pathway = parser.add_argument_group("Nameing Parameters", "")
    amplicon_pathway.add_argument(
        '-n', '--file_prefix_name', 
        help='File name prefix for all data generated [default="SEWAGE_"]', 
        required=False,
        default=None,
        type=str,
        metavar='STR',
        dest='file_prefix_name')
    amplicon_pathway.add_argument(
        '-sd', '--storage_dir', 
        help='Directory name for data storage [default="SEWAGE_[data_time]"]', 
        required=False,
        default=None,
        type=str,
        metavar='STR',
        dest='amplicon_storage_dir')

    # Create a group for the "Propirtions" section
    proportions_options = parser.add_argument_group("Proportion options", "")
    proportions_options.add_argument(
        "-p", "--proportion_model", 
        help="Generate equal (e), random (r), or dominate (d) variant of concern proportions of reads [default: d]", 
        required=False,
        default='r',
        type=str,
        choices=['r', 'e', 'd'],
        dest="proportion_model",
    )
    proportions_options.add_argument(
        "-dg", "--dVOC_genome",
        help='Name of dVOC. NOTE: name of dVOC must match the defline of the reference fasta file',
        required=False,
        default=None,
        type=str,
        metavar="STR",
        dest="dVOC_genome"
    )    
    proportions_options.add_argument(
        "-dp", "--dVOC_proporiton",
        help="Proportion of dDOV [default: 0.8]",
        required=False,
        default=0.8,
        metavar="FLOAT",
        type=float,
        dest="dVOC_proporiton",
    )    
    proportions_options.add_argument(
        "-ps", "--proportion_seed", 
        help="Random seed number for reproducing proportions [default: 13]", 
        metavar='INT',
        default=13,
        type=int,
        dest='proportion_seed'
    )
    read_options = parser.add_argument_group("Read generator options", "")
    read_options.add_argument(
        "-q", "--fastq_name", 
        help="Name of fastq files for F/R reads [default: SEWAGE_{R1/R2}.fastq]", 
        metavar='STR',
        default=None,
        type=str,
        dest='fastq_name'
    )
    read_options.add_argument(
        "-rl", "--read_length", 
        help="Read length in bp (value should not exceed amplicon length or will workflow fail) [default: 250]", 
        metavar='INT',
        default=250,
        type=int,
        dest='read_length'
    )
    read_options.add_argument(
    "-auto", "--auto_read_length_detection", 
    help="Automatically set the read length to the maximum possible based on amplicon length (i.e., half max amplicon + 50bp).  [default: False]",
    action='store_true',
    default=False,
    dest='auto_read_length_detection'
    )
    read_options.add_argument(
        "-cd", "--coverage_depth", 
        help="Total sequence depth coverage for each fastq file [default: 500]", 
        metavar='INT',
        default=500,
        type=int,
        dest='coverage_depth'
    )
    # self.max_reads
    read_options.add_argument(
        "-mr", "--max_reads", 
        help="Total number of reads for each fastq file [default: None]. NOTE: If set, --coverage_depth is ignored.", 
        metavar='INT',
        default=None,
        type=int,
        dest='max_reads'
    )
    read_options.add_argument(
        "-rs", "--read_seed", 
        help="Random seed number for reproducing reads [default: 13]",
        metavar='INT',
        default=13,
        type=int,
        dest='read_seed'
    )
    args = parser.parse_args()
    return args

def write_parameters_log(args, file_prefix_name, storage_dir):
    args_dict = vars(args)
    args_str = '\n'.join(f'{key}: {value}' for key, value in args_dict.items())
    if file_prefix_name is not None:
        log_file_path = os.path.join(storage_dir, f'{file_prefix_name}_parameters.txt')
    else:
        log_file_path = os.path.join(storage_dir, f'SEWAGE_parameters.txt')
    with open(log_file_path, 'w') as file:
        file.write(args_str + '\n')