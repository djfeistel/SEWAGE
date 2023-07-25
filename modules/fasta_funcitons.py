import glob
import os 
import sys

def globFastAfiles(pathway)->list:
    '''glob fasta files in pathway, directory must only harbor fasta files'''
    pathway = os.path.realpath(pathway)
    file_pathway_list = glob.glob(pathway + "/*")
    valid_extensions = ('.fasta', '.fa', '.fsa')
    fasta_pathway_list = [file for file in file_pathway_list if os.path.splitext(file)[1].lower() in valid_extensions]
    if len(fasta_pathway_list) == 0:
        print(f"\nThere are not fasta files located in\n\t{pathway}\n", file=sys.stderr)
        print(f"Check to make sure that your fasta files end with .fasta, .fa, or .fsa", file=sys.stderr)
        sys.exit(0)
    #print(f"Total FASTA files detected in pathway are: {len(fasta_pathway_list)}", file=sys.stderr)
    return fasta_pathway_list

def listOfFastaPathways(fasta_in):
    with open(fasta_in, 'r') as fh:
        #pathways_list = [line.strip().rsplit() for line in file]
        pathways_list = fh.read().splitlines()
    return pathways_list