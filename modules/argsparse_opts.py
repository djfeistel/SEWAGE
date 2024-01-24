import argparse
import os

def sewage_opts():
    parser = argparse.ArgumentParser(
        prog="SEWAGE",
        description=f"\nSynthetically Engineered Wastewater sequence data for Assessing Genomic Entities\n",
        epilog=f"Minimal Usage: SEWAGE -i <multi.fasta> -s <scheme>"
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
    # output options #
    output_pathway = parser.add_argument_group("Output naming", "")
    output_pathway.add_argument(
        '-n', '--file_prefix_name', 
        help='File name prefix for all generated data [default="SEWAGE"]', 
        required=False,
        default="SEWAGE",
        type=str,
        metavar='STR',
        dest='file_prefix_name')
    output_pathway.add_argument(
        '-o', '--storage_dir', 
        help='Directory name for data storage [default="SEWAGE_Workspace"]', 
        required=False,
        default="SEWAGE_Workspace",
        type=str,
        metavar='STR',
        dest='storage_dir')
    output_pathway.add_argument(
        "-t", "--time_stamp", 
        help="Append date and time stamp to the storage directory [default: False]",
        action='store_true',
        default=False,
        dest='time_stamp'
    )
    # Propirtion options #
    proportions_options = parser.add_argument_group("Proportion options", "")
    proportions_options.add_argument(
        "-p", "--proportion_model", 
        help="Generate equal (e), random (r), or dominate (d) variant of concern proportions of reads [default: r]", 
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
    # read options #
    read_options = parser.add_argument_group("Read generator options", "")
    read_options.add_argument(
        "-rl", "--read_length", 
        help="Read length in bp (value should not exceed amplicon length or will workflow fail) [default: 150]", 
        metavar='INT',
        default=150,
        type=int,
        dest='read_length'
    )
    read_options.add_argument(
        "-fl", "--frag_length", 
        help="Fragment length in bp (value should not exceed amplicon length or will workflow fail) [default: 300]", 
        metavar='INT',
        default=300,
        type=int,
        dest='frag_length'
    )
    read_options.add_argument(
        "-cd", "--coverage_depth", 
        help="Total sequence depth coverage for each fastq file [default: 500]", 
        metavar='INT',
        default=500,
        type=int,
        dest='coverage_depth'
    )
    read_options.add_argument(
        "-rs", "--read_seed", 
        help="Random seed number for reproducing reads [default: 13]",
        metavar='INT',
        default=13,
        type=int,
        dest='read_seed'
    ) #num_cpus
    # read_options.add_argument(
    #     "-cpus", "--num_cpus", 
    #     help="Number of CPUs for parallel processing (max CPUs allowed are (total_cpus - 1)) [default: 1]",
    #     metavar='INT',
    #     default=1,
    #     type=int,
    #     dest='num_cpus'
    # )
    # Quality Score Parameters
    q_score_options = parser.add_argument_group("Read generator options", "")
    q_score_options.add_argument(
        "-m", "--min_max_q", 
        help="Minimum and Maximum mean range quality scores [default: 30 40]",
        nargs=2,
        metavar=('INT', 'INT'),
        default=[30, 40],
        type=int,
        dest='min_max_q'
    )
    q_score_options.add_argument(
        "-sd", "--std_dev_q", 
        help="Standard deviation for quality scores [default: 3]",
        metavar='INT',
        default=3,
        type=int,
        dest='std_dev_q'
    )
    q_score_options.add_argument(
        "-bp", "--startup_effect_bp", 
        help="Number of bp's at the begining of a read that have reduced quality [default: 10]",
        metavar='INT',
        default=10,
        type=int,
        dest='startup_effect_bp'
    )
    q_score_options.add_argument(
        "-red", "--startup_effect_q_reduction", 
        help="Quality score reduction for '-bp' (maximum mean minus '-red') [default: 4]",
        metavar='INT',
        default=4,
        type=int,
        dest='startup_effect_q_reduction'
    )
    
    args = parser.parse_args()
    return args

def write_parameters_log(args, file_prefix_name, storage_pathway):
    args_dict = vars(args)
    args_str = '\n'.join(f'{key}: {value}' for key, value in args_dict.items())
    
    log_file_path = os.path.join(storage_pathway, f'{file_prefix_name}_parameters.txt')

    with open(log_file_path, 'w') as file:
        file.write(args_str + '\n')