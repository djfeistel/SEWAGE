import glob
import os 
'''
gather and load reference genomes and ictionary with defline as key and sequence as value
note: increasing the number of genomes or geneoms with longer sequeces will make this load slower
'''

def check_fasta_format(fasta:str)->bool:
    with open(fasta, 'r') as fh:
        if fh.readline().startswith(">"):
            return True
        else:
            return False

def glob_fasta_files_in_cwd(extention:str) -> list:
    '''
    used for getting fasta files in directory
    note that fasta files being used need to have the same extention!
    '''
    return glob.glob(os.path.realpath('.') + f"/*{extention}")

def fasta_pathway_list(pathway_list:str) -> list:
    '''used for a fasta pathway list'''
    with open(pathway_list, 'r') as fh:
        fasta_pathways_list = fh.read().splitlines()
    return fasta_pathways_list

def concatenate_fasta_files(fasta_paths:list):
    '''creates a variable with fasta files as a multi fasta for glob and pathway list'''
    for fasta_path in fasta_paths:
        with open(fasta_path, 'r') as file:
            concatenated_fasta += file.read() + "\n"  # Adding a newline character between files
    return concatenated_fasta

def read_fasta_into_dict(fasta, fasta_dict:dict) -> dict:
    '''use this with mulitfasta, or concatenated_fasta'''
    current_defline = ""
    current_sequence = ""
    with open(fasta, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):  # Defline
                if current_defline: # if true
                    fasta_dict[current_defline] = current_sequence # seqeuence is done, add to dict
                current_defline = line[1:] # remove > 
                current_sequence = "" #restart next sequence 
            else:
                current_sequence += line #assumes that sequence is broken up on new lines
        # Add the last sequence after reaching the end of file
        if current_defline and current_sequence:
            fasta_dict[current_defline] = current_sequence
    return fasta_dict

