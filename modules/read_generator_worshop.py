import numpy as np
import pandas as pd
import random
import sys
from fasta_funcitons_class_workshop import load_fasta_workflow_test
from amplicons_functions_class_worshop import generate_amplicon_workflow_test
from proportion_functions_workshop import GenomeProporitons

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
    
    def fragment_indices_list_df(
                self,
                frag_length: int=250,
                read_length: int=150,
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
                upper_limit = int(row['length_bp'] + 1 - frag_length)

                # Generate indices
                forward_indices = np.sort(np.random.randint(0, upper_limit, size=total_reads))
                reverse_indices = forward_indices + frag_length
                read_indicies = list(zip(forward_indices, reverse_indices))
                
                # Store indices in DataFrame
                self.df_amplicon.at[idx, 'read_indicies'] = read_indicies
            print(self.df_amplicon)
    @staticmethod
    def name_fastq_read(
            amplicon_defline: str,
            read_sequence:str, 
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

        defline = f"@{amplicon_defline}:{read_type}_{read_count}"
        
        yield [defline, read_sequence, "+", "I"*len(read_sequence)]
    
    @staticmethod
    def write_fastq(amplicon_reads_list:dict, fastq_name:str):
        with open(fastq_name, "w") as wf, open(amplicon_reads_list, 'r') as reads:
            for read in reads:
                read_out = '\n'.join(read)
                #read_out = f"{read[0]}\n{read[1]}\n{read[2]}\n{read[3]}"
                wf.write(read_out)


    # def read_slicing()




########
def read_generator_workflow_test():
    reference_fasta_dict = load_fasta_workflow_test(sys.argv[1])
    
    df_amplicon = generate_amplicon_workflow_test(sys.argv[1])
    proportion_dict = GenomeProporitons(reference_fasta_dict=reference_fasta_dict).randomProportions()
    
    rg = ReadGenerator(df_amplicon=df_amplicon, proportion_dict=proportion_dict)
    rg.add_proportion_column()
    rg.fragment_indices_list_df()
    print(df_amplicon['total_reads'].sum())


if __name__ == "__main__":
    read_generator_workflow_test()


