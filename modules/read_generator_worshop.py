import numpy as np
import pandas as pd
import random
import sys
from fasta_funcitons_class_workshop import load_fasta_workflow_test
from amplicons_functions_class_worshop import generate_amplicon_workflow_test
from proportion_functions_workshop import GenomeProporitons


class ReadGenerator():
    def __init__(
            self,
            df_amplicon:pd.DataFrame,
            proportion_dict:dict
    ):
        self.df_amplicon = df_amplicon
        self.proportion_dict = proportion_dict
    
    def add_proportion_column(self):
        self.df_amplicon['proportion'] = self.df_amplicon['reference_defline'].map(self.proportion_dict)
    
    def read_indices_og(
                self,
                read_length: int=250,
                coverage_depth: int=60,
                seed=None
        ) -> list:
            '''
            creates a sorted list of random numbers which are used 
                as the indicies for slicing an amplicon into reads
            one amplicion at a time
            '''
            np.random.seed(seed)

            # calulate total reads need per amplicon
            self.df_amplicon['total_reads'] = coverage_depth * self.df_amplicon['length_bp'] // read_length
            
            # Initialize columns for forward and reverse read indices
            self.df_amplicon['read_indicies'] = None
            
            # Iterate over DataFrame rows
            for idx, row in self.df_amplicon.iterrows():
                total_reads = int(row['total_reads'])
                upper_limit = int(row['length_bp'] + 1 - read_length)

                # Generate indices
                forward_indices = np.sort(np.random.randint(0, upper_limit, size=total_reads))
                reverse_indices = forward_indices + read_length
                read_indicies = list(zip(forward_indices, reverse_indices))
                
                # Store indices in DataFrame
                self.df_amplicon.at[idx, 'read_indicies'] = read_indicies
            self.df_amplicon.to_csv("temp.tsv", sep='\t', index=False)
            print(self.df_amplicon)

    def create_reads_from_amplicons(
                    self,
                    read_length: int=250,
                    coverage_depth: int=60
                    ):
                '''
                Create reads from amplicons
                '''

                # calulate total reads need per amplicon for desired sequencing depth
                self.df_amplicon['total_reads'] = coverage_depth * self.df_amplicon['length_bp'] // read_length
                
                self.df_amplicon['R1_read'] = None
                self.df_amplicon['R2_read'] = None

                # Iterate over DataFrame rows to get reads
                for idx, row in self.df_amplicon.iterrows():

                    R1_read = row['amplicon_sequence'][:read_length]
                    R2_index_start = row['length_bp'] - read_length
                    R2_read = row['amplicon_sequence'][R2_index_start:]
                    self.df_amplicon.at[idx, 'R1_read'] = R1_read
                    self.df_amplicon.at[idx, 'R2_read'] = R2_read
                
    def write_reads(self):
        for idx, row in self.df_amplicon.iterrows():
            amplicon_defline = ':'.join([
                row['reference_defline'], 
                row['primer_scheme'], 
                row['primer_name']
            ])

            R1 = row['R1_read']
            R2 = row['R2_read']

            for read_number in range(1,row['total_reads']+1):
                #forward reads
                fastq_line_R1 = ReadGenerator.name_fastq_read(
                    amplicon_defline=amplicon_defline,
                    read_sequence=R1,
                    read_number=read_number
                    )
                fastq_line_R2 = ReadGenerator.name_fastq_read(
                    amplicon_defline=amplicon_defline,
                    read_sequence=R2,
                    read_number=read_number,
                    reverse=True
                    )
                ReadGenerator.write_fastq(R1=fastq_line_R1, R2=fastq_line_R2, fastq_name="temp")
    
    @staticmethod
    def name_fastq_read(
            amplicon_defline: str,
            read_sequence:str, 
            read_number:int, 
            reverse: bool=False
            ) -> list:
        
        if reverse:
            read_type = "R2"
        else:
            read_type = "R1"

        defline = f"@{amplicon_defline}:{read_type}_{read_number}"
        #currently there is no quality score implimentation 
        return [defline, read_sequence, "+", "I"*len(read_sequence)]
    
    @staticmethod
    def write_fastq(R1, R2, fastq_name:str):
        with open(fastq_name + "_R1.fastq", "a") as wf_r1, open(fastq_name + "_R2.fastq", "a") as wf_r2:
            R1_output = f"{R1[0]}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n"
            R2_output = f"{R2[0]}\n{R2[1]}\n{R2[2]}\n{R2[3]}\n"
            wf_r1.write(R1_output)
            wf_r2.write(R2_output)

    # def read_slicing()




########
def read_generator_workflow_test():
    reference_fasta_dict = load_fasta_workflow_test(sys.argv[1])
    
    df_amplicon = generate_amplicon_workflow_test(sys.argv[1])
    proportion_dict = GenomeProporitons(reference_fasta_dict=reference_fasta_dict).randomProportions()
    
    rg = ReadGenerator(df_amplicon=df_amplicon, proportion_dict=proportion_dict)
    # combine propriton values on defline columns
    rg.add_proportion_column()
    # get list of indicies for reads
    rg.create_reads_from_amplicons()
    # create reads
    rg.write_reads()


if __name__ == "__main__":
    read_generator_workflow_test()


