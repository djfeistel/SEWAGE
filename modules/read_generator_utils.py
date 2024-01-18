import numpy as np
import pandas as pd
import os
import sys
# from fasta_funcitons_class_workshop import load_fasta_workflow_test
# from amplicons_functions_class_worshop import generate_amplicon_workflow_test
# from proportion_functions_workshop import GenomeProporitons


class ReadGenerator():
    def __init__(
            self,
            df_amplicon:pd.DataFrame,
            proportion_dict:dict,
            storage_dir:str,
            fastq_name:str=None,
            read_length:int=250,
            coverage_depth:int=60,
            seed:int=13
    ):
        self.df_amplicon = df_amplicon
        self.proportion_dict = proportion_dict
        self.storage_dir = storage_dir
        self.fastq_name = fastq_name
        self.read_length = read_length
        self.coverage_depth = coverage_depth
        self.seed = seed

    def add_proportion_column(self):
        self.df_amplicon['proportion'] = self.df_amplicon['reference_defline'].map(self.proportion_dict)

    def create_reads_from_amplicons(self):
                '''
                Create reads from amplicons
                '''

                # calulate total reads need per amplicon for desired sequencing depth
                self.df_amplicon['total_reads'] = self.coverage_depth * self.df_amplicon['length_bp'] // self.read_length
                
                self.df_amplicon['R1_read'] = None
                self.df_amplicon['R2_read'] = None

                # Iterate over DataFrame rows to get reads
                for idx, row in self.df_amplicon.iterrows():

                    R1_read = row['amplicon_sequence'][:self.read_length]
                    R2_index_start = int(row['length_bp'] - self.read_length)

                    R2_read = row['amplicon_sequence'][R2_index_start:]
                    self.df_amplicon.at[idx, 'R1_read'] = R1_read
                    self.df_amplicon.at[idx, 'R2_read'] = R2_read
    
    def write_fastq(self, R1, R2):
        
        if self.fastq_name is None:
             fastq_name_out = "SEWAGE"
        else:
             fastq_name_out = self.fastq_name

        fastq_r1 = os.path.join(self.storage_dir, fastq_name_out + "_R1.fastq")
        fastq_r2 = os.path.join(self.storage_dir, fastq_name_out + "_R2.fastq")
        with open(fastq_r1, "w") as wf_r1, open(fastq_r2, "w") as wf_r2:
            for i in range(len(R1)):
                R1_output = f"{R1[i][0]}\n{R1[i][1]}\n{R1[i][2]}\n{R1[i][3]}\n"
                R2_output = f"{R2[i][0]}\n{R2[i][1]}\n{R2[i][2]}\n{R2[i][3]}\n"
                wf_r1.write(R1_output)
                wf_r2.write(R2_output)

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
       
    def write_reads(self):
        R1_list = []
        R2_list = []        
        for idx, row in self.df_amplicon.iterrows():
            amplicon_defline = ':'.join([
                row['reference_defline'], 
                row['primer_scheme'], 
                row['primer_name']
            ])

            R1 = row['R1_read']
            R2 = row['R2_read']

            for read_number in range(1, row['total_reads']+1):
                #forward reads
                fastq_line_R1 = ReadGenerator.name_fastq_read(
                    amplicon_defline=amplicon_defline,
                    read_sequence=R1,
                    read_number=read_number
                    )
                R1_list.append(fastq_line_R1)
                fastq_line_R2 = ReadGenerator.name_fastq_read(
                    amplicon_defline=amplicon_defline,
                    read_sequence=R2,
                    read_number=read_number,
                    reverse=True
                    )
                R2_list.append(fastq_line_R2)
        self.write_fastq(R1=R1_list, R2=R2_list)
    
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


