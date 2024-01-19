import os 
import textwrap
'''
gather and load reference genomes and create dictionary with defline as key and sequence as value
note: increasing the number of genomes or geneoms with longer sequeces may make this load slower
'''

class FastaUtils():
    '''
    class function for loading fasta genome data into dictionary
        key: defline
        value: sequence
    '''
    def __init__(
            self,
            fasta:str,
            file_prefix_name:str="SEWAGE",
            
    ):
        self.fasta = fasta
        self.file_prefix_name = file_prefix_name
        self.fasta_dict = None

    def get_fasta_dict(self):
        return self.fasta_dict
        
    def check_input_for_fasta_or_pathways(self):
        '''open file (fasta or pathway list) and test if first line if fasta defline or pathway file'''
        with open(self.fasta, 'r') as fh:
            if fh.readline().strip().startswith(">"):
                return True
            elif os.path.isfile(fh.readline().strip()):
                return False
            
    def fasta_pathway_list(self) -> list:
        '''used for a fasta pathway list'''
        with open(self.fasta, 'r') as fh:
            self.fasta = fh.read().splitlines()

    def check_pathway_list_for_fasta(self):
        for path in self.fasta:
            try:
                with open(path, 'r') as file:
                    first_line = file.readline().strip()
                    if not first_line.startswith(">"):
                        return False
            except Exception:
                return False
        return True

    def create_list_of_pathways(self)->list:
        pathway_list = []
        with open(self.fasta, 'r') as fasta_pathways:
            for fasta_pathway in fasta_pathways:
                fasta_pathway = fasta_pathway.strip()
                pathway_list.append(fasta_pathway)
        return pathway_list

    def read_multifasta_into_dict(self) -> dict:
        '''use this with mulitfasta, or concatenated_fasta'''
        self.fasta_dict = {}
        current_defline = ""
        current_sequence = ""
        with open(self.fasta, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                if line.startswith(">"):  # Defline
                    if current_defline: # if true
                        self.fasta_dict[current_defline] = current_sequence # seqeuence is done, add to dict
                    current_defline = line[1:] # remove > 
                    current_sequence = "" #restart next sequence 
                else:
                    current_sequence += line #assumes that sequence is broken up on new lines
            # Add the last sequence after reaching the end of file
            if current_defline and current_sequence:
                self.fasta_dict[current_defline] = current_sequence
          
    def read_fasta_pathways_into_dict(self, fasta_pathway_list) -> dict:
            '''use this with fasta pathway list'''
            self.fasta_dict = {}
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
                                self.fasta_dict[current_defline] = current_sequence # seqeuence is done, add to dict
                            current_defline = line[1:] # remove > 
                            current_sequence = "" #restart next sequence 
                        else:
                            current_sequence += line #assumes that sequence is broken up on new lines
                    # Add the last sequence after reaching the end of file
                    if current_defline and current_sequence:
                        self.fasta_dict[current_defline] = current_sequence
            
    def write_reference_fasta(self, storage_pathway):
        file_output = f"{self.file_prefix_name}_Reference_genomes.fasta"
        with open(os.path.join(storage_pathway, file_output), "w") as wf:
            for genome, sequence in self.fasta_dict.items():
                wf.write(f">{genome}\n")
                wrapped_sequence = textwrap.wrap(sequence, width=60)
                for line in wrapped_sequence:
                    wf.write(f"{line}\n")
