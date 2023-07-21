import argparse

def opts():
    parser = argparse.ArgumentParser(
        prog="SEWAGE.py",
        description=f"\nSimulation of Environmental Wastewater sequence data for the Analysis of Genomics and Epidemiology\n",
        epilog=f"Minimal usage: ./SEWAGE -i <input>"
    )
    '''input and output'''
    required_pathway = parser.add_argument_group("Input and Output Parameters", "-i/--in is required")
    required_pathway.add_argument(
        "-i", "--in", 
        help="Pathway to directory with FASTA files or a text file with a list of pathways to FASTA files. \
                NOTE: FASTA files must end in .fasta, .fa, or .fsa when a pathway is specified.", 
        required=False,
        dest='fasta',
        type=str,
        default=None,
        metavar='PATHWAY or FILE'
    )
    required_pathway.add_argument(
        "-o", "--out", 
        help="Name of output directory for storage (if not specificed, default is 'SEWAGE_' + 10 random alphanumeric characters).", 
        required=False,
        dest='out',
        metavar='STR',
        type=str,
        default=None
    )
    required_pathway.add_argument(
        "-O", "--out_pathway", 
        help="Pathway to where output directory is stored [default='.']", 
        required=False,
        dest='out_pathway',
        metavar='DIR',
        type=str,
        default='.'
    )

    # Create a group for the "Propirtions" section
    proportions_options = parser.add_argument_group("Proportion options [default is -p r -rs 13]", "")
    proportions_options.add_argument(
        "-p", "--proportion", 
        help="Generate random (r) or equal (e) proportions of reads", 
        required=False,
        dest="proportion",
        default='r',
        type=str,
        choices=['r', 'e']
    )    
    proportions_options.add_argument(
        "-rs", "--rndSeed", 
        help="Random seed for generateing proportions [default=13]", 
        metavar='INT',
        default=13,
        type=int,
        dest='rndSeed'
    )
    '''art_illumina paramters''' 
    art_parameters = parser.add_argument_group(
        "ART parameters", 'Default parameters listed below are for \
        simulating "perfect" reads at 150bp and can be modified as needed. \
        All other parameters not listed here are in default setting or not used as defined by "art_illumia" \
        and cannot be access via SEWAGE. Please be familiar with how "art_illumina" \
        functins before modifying these parameters')
    art_parameters.add_argument(
        "-pf", "--pfold", 
        help="Value to be mutiplied by proportion for use with '--fcov' from 'art_illumina' [default=1000]\
            Example: proportion*pfold=fcov or fold coverage", 
        metavar='INT',
        default=1000,
        type=int,
        dest='pfold'
    )
    art_parameters.add_argument(
        "-ss", "--seqSys", 
        help="From 'art_illumina': 'The name of Illumina sequencing system of the built-in profile used for \
            simulation' [default=HS25]. Note: chosing a differnt Illumina sequecning system may require \
                modifying the 'art_illumina' paramters beforehand", 
        #metavar='STR',
        default='HS25',
        type=str,
        dest='ss',
        choices=['HS10', 'HS20', 'HS25', 
                 'HSXn', 'HSXt', 'MinS', 
                 'MSv1', 'MSv3', 'NS50', 
                 'GA1', 'GA2']
    )
    art_parameters.add_argument(
        "-l", "--len", 
        help="From 'art_illumina': the length of reads to be simulated [default=150]", 
        metavar='INT',
        default=150,
        type=int,
        dest='l'
    )
    art_parameters.add_argument(
        "-m", "--mflen", 
        help="From 'art_illumina': the mean size of DNA/RNA fragments for paired-end simulations [default=250]", 
        metavar='INT',
        default=250, #change from 200 to 250
        type=int,
        dest='m'
    )
    art_parameters.add_argument(
        "-s", "--sdev", 
        help="From 'art_illumina': the standard deviation of DNA/RNA fragment size for paired-end simulations [default=1]", 
        metavar='INT',
        default=1,
        type=int,
        dest='s'
    )
    art_parameters.add_argument(
        "-ir", "--insRate", 
        help="From 'art_illumina': the first-read insertion rate [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='ir'
    )
    art_parameters.add_argument(
        "-ir2", "--insRate2", 
        help="From 'art_illumina': the second-read insertion rate [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='ir2'
    )
    art_parameters.add_argument(
        "-dr", "--delRate", 
        help="From 'art_illumina': the first-read deletion rate [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='dr'
    )
    art_parameters.add_argument(
        "-dr2", "--delRate2", 
        help="From 'art_illumina': the second-read deletion rate [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='dr2'
    )
    art_parameters.add_argument(
        "-rsA", "--rndSeed_art_illumina", 
        help="From 'art_illumina': the seed for random number generator [default=13]", 
        metavar='INT',
        default=13,
        type=int,
        dest='rndSeed_art_illumina'
    )
    art_parameters.add_argument(
        "-k", "--maxIndel", 
        help="From 'art_illumina': the maximum total number of insertion and deletion per read [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='maxIndel'
    )
    art_parameters.add_argument(
        "-nf", "--maskN", 
        help="From 'art_illumina': the cutoff frequency of 'N' in a window size of the read length for masking genomic regions [default=0]", 
        metavar='INT',
        default=0,
        type=int,
        dest='nf'
    )
    art_parameters.add_argument(
        "-qL", "--minQ", 
        help="From 'art_illumina': the minimum base quality score [default=28]", 
        metavar='INT',
        default=28,
        type=int,
        dest='qL'
    )
    art_parameters.add_argument(
        "-qU", "--maxQ", 
        help="From 'art_illumina': the maximum base quality score [default=40]", 
        metavar='INT',
        default=40,
        type=int,
        dest='qU'
    )
    art_parameters.add_argument(
        "-qs", "--qShift", 
        help="From 'art_illumina': the amount to shift every first-read quality score by", 
        metavar='INT',
        type=int,
        dest='qs'
    )
    art_parameters.add_argument(
        "-qs2", "--qShift2", 
        help="From 'art_illumina': the amount to shift every second-read quality score by", 
        metavar='INT',
        type=int,
        dest='qs2'
    )
    tool_description = parser.add_argument_group("Tool Description", "Detailed information about the tool")
    tool_description.add_argument(
        "--details", 
        help="", 
        required=False,
        dest="details",
        action='store_true'
    )

    args = parser.parse_args()
    return args