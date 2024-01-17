import glob
import os 
import sys
'''
gather and load reference genomes and create dictionary with defline as key and sequence as value
note: increasing the number of genomes or geneoms with longer sequeces may make this load slower
'''

class LoadFastaClass():
    '''
    class function for loading fasta genome data into dictionary
        key: defline
        value: sequence
    '''
    def __init__(
            self,
            fasta:str,
    ):
        self.fasta = fasta
        
    
    def check_input_for_fasta_or_pathways(self):
        '''open file (fasta or pathway list) and test if first line if fasta defline or pathway file'''
        with open(self.fasta, 'r') as fh:
            if fh.readline().strip().startswith(">"):
                return True
            elif os.path.isfile(fh.readline().strip()):
                return False

    def glob_fasta_files_in_cwd(self) -> list:
        '''
        used for getting fasta files in directory
        note that fasta files being used need to have the same extention!
        '''
        return glob.glob(os.path.realpath('.') + f"/*{self.extention}")

    def fasta_pathway_list(self) -> list:
        '''used for a fasta pathway list'''
        with open(self.fasta, 'r') as fh:
            fasta_pathways_list = fh.read().splitlines()
        return fasta_pathways_list

    def create_list_of_pathways(self)->list:
        pathway_list = []
        with open(self.fasta, 'r') as fasta_pathways:
            for fasta_pathway in fasta_pathways:
                fasta_pathway = fasta_pathway.strip()
                pathway_list.append(fasta_pathway)
        return pathway_list

    def read_multifasta_into_dict(self) -> dict:
        '''use this with mulitfasta, or concatenated_fasta'''
        fasta_dict = {}
        current_defline = ""
        current_sequence = ""
        with open(self.fasta, 'r') as fh:
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
    
    def read_fasta_pathways_into_dict(self, fasta_pathway_list) -> dict:
            '''use this with fasta pathway list'''
            fasta_dict = {}
            current_defline = ""
            current_sequence = ""
            for fasta_pathway in fasta_pathway_list:    
                with open(fasta_pathway, 'r') as fh:
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

def load_fasta_workflow(infasta:str)->dict:
    '''
    workflow for creating dictionary with genomes
        key: genome
        value: sequence
    '''

    fasta_data = LoadFastaClass(fasta=infasta)

    if fasta_data.check_input_for_fasta_or_pathways():
        fasta_dict = fasta_data.read_multifasta_into_dict()
        print(fasta_dict)
    else:
        pathway_list = fasta_data.create_list_of_pathways()
        for i in pathway_list:
            print(i)
        fasta_dict = fasta_data.read_fasta_pathways_into_dict(pathway_list)
        print(fasta_dict)

if __name__ == "__main__":
    # use this to run work flow and test data 
    load_fasta_workflow(infasta=sys.argv[1])