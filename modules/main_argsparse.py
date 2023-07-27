import argparse
from enrichment_workflow import EnrichmentWorkflow
from amplicons_functions import generate_amplicons
from detailed_message import details

def main_sewage_menu():
    # Create the main parser
    parser = argparse.ArgumentParser(prog="SEWAGE",
                                     description=f"Synthetically Engineered Wastewater sequence data for Assessing Genomic Entities\n",
                                     epilog="")

    # Create subparsers
    subparsers = parser.add_subparsers(title='Options', 
                                       dest='subcommand', 
                                       description='', 
                                       required=True,
                                       metavar='Select one option')

    # Subparser for "amplicon" subcommand
    amplicon_parser = subparsers.add_parser('amplicon', help='Generate amplicons from a primer scheme')

    amplicon_parser.add_argument('-f', '--fasta', 
                                 help='Single or multi-fasta reference. Multi-fasta files should be unique genomes for each defline/sequence. \
                                    Output will produce "${defline}_amplicon.fasta" files for each genome', 
                                 required=True,
                                 type=str,
                                 dest='fasta')
    amplicon_parser.add_argument('-s', '--scheme', 
                                 help='Primer scheme: (Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"], \
                                    VarSkip = ["vsl1a", "vss1a", "vss2a", "vss2b"])', 
                                 required=True, 
                                 type=str,
                                 choices=["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vsl1a", "vss1a", "vss2a", "vss2b"],
                                 metavar='SCHEME',
                                 dest='scheme')
    amplicon_parser.add_argument('-o', '--output', 
                                 help='Output Prefix name for fasta file [default="SEWAGE_amplicons.fasta"]', 
                                 required=False,
                                 default="SEWAGE_amplicons",
                                 type=str,
                                 metavar='STR',
                                 dest='output')
    amplicon_parser.add_argument('-p', '--pathway', 
                                 help='Pathway to storgage directory [default="."]',
                                 required=False, 
                                 type=str,
                                 metavar='PATHWAY',
                                 dest='pathway',
                                 default='.')
    amplicon_parser.set_defaults(func=generate_amplicons)

    # Subparser for "sequence" subcommand
    enrich_parser = subparsers.add_parser('enrich', help='Generate enriched sequence data at equal or unequal proporitons from amplicons')
    required_pathway = enrich_parser.add_argument_group("Input and Output Parameters", "-i/--in is required")
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
        help="Name of output directory for storage [default='SEWAGE_enrich'].", 
        required=False,
        dest='out',
        metavar='STR',
        type=str,
        default="SEWAGE_enrich"
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
    proportions_options = enrich_parser.add_argument_group("Proportion options [default is -p v -V 0.8 -rs 13]", "")
    proportions_options.add_argument(
        "-p", "--proportion", 
        help="Generate a dVOC (v), random (r), or equal (e) proportion read set", 
        required=False,
        dest="proportion",
        default='v',
        type=str,
        choices=['v', 'r', 'e']
    )
    proportions_options.add_argument(
        "-V", "-dVOC",
        help="Proporiton of dVOC  [default=0.8]. NOTE: End proportion value might vary slighlty. \
            Valuse >= 1 are converted to 0.99",
        required=False,
        dest="dVOC",
        type=float,
        default=0.8,
        metavar="FLOAT"

    )
    proportions_options.add_argument(
        "-rs", "--rndSeed", 
        help="Random seed used for generating dVOC and random proportions [default=13]", 
        metavar='INT',
        default=13,
        type=int,
        dest='rndSeed'
    )
    '''art_illumina paramters''' 
    art_parameters = enrich_parser.add_argument_group(
        "ART parameters", 'Default parameters listed below are for \
        simulating "perfect" reads at 150bp and can be modified as needed. \
        All other parameters not listed below are either not in use or in default \
        setting as defined by "art_illumia" and cannot be access via SEWAGE. \
        Please be familiar with how "art_illumina" functins before modifying these parameters')
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
    # details about tool
    enrich_parser.set_defaults(func=EnrichmentWorkflow)

    details_parser = subparsers.add_parser('details', help='Details about this tool and other information')
    details_parser.set_defaults(func=details)
    
    args = parser.parse_args()

    return args
