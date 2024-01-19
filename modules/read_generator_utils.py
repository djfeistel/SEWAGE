import numpy as np
import pandas as pd
import os
import sys

class ReadGenerator():
    def __init__(
            self,
            df_amplicon:pd.DataFrame,
            proportion_dict:dict,
            fastq_name:str="SEWAGE",
            propotion_file_name:str="SEWAGE",
            read_length:int=150,
            frag_length:int=300,
            coverage_depth:int=500,
            seed:int=13,
    ):
        self.df_amplicon = df_amplicon
        self.df_reads = df_amplicon.copy()
        self.df_reads = self.df_reads.dropna(subset=['start', 'end'], how='all')
        self.proportion_dict = proportion_dict
        self.fastq_name = fastq_name
        self.propotion_file_name = propotion_file_name
        self.read_length = read_length
        self.frag_length = frag_length
        self.coverage_depth = coverage_depth
        self.seed = seed
        
    def get_df_reads(self):
        return self.df_reads

    def map_proportions_to_df(self):
        self.df_reads['proportion'] = self.df_reads['reference_defline'].map(self.proportion_dict)
         
    def calculate_total_frags_to_cover_amplicon(self): 
        total_frags = (self.coverage_depth * self.df_reads['amplicon_length_bp'] // self.frag_length)
        self.df_reads['total_frags'] = (self.df_reads['proportion'] * total_frags).astype(int)
    
    def generate_random_fragment_indicies(self, row):
        total_frags = row['total_frags']
        amplicon_length_bp = row['amplicon_length_bp']
        start_inx = np.random.randint(low=0, high=(amplicon_length_bp - self.frag_length), size=total_frags)
        end_inx = start_inx + self.frag_length
        return sorted(list(zip(start_inx, end_inx)), key=lambda x: x[0])

    def slice_fragments_into_reads(self, row):
        #return list of reads
        return [[row['amplicon_sequence'][start:end][:self.read_length], row['amplicon_sequence'][start:end][self.read_length:]] for start, end in row['fragment_indices']]

    @staticmethod    
    def reverse_compliment(reverse_strand:str):
        compliment = {"A": "T", "T": "A", "G": "C", "C": "G"}
        complement_strand = [compliment[base] if base in compliment else 'N' for base in reverse_strand]
        return complement_strand[::-1]
    
    def create_reads_workflow(self):

        self.map_proportions_to_df()
        self.calculate_total_frags_to_cover_amplicon()
        self.df_reads['fragment_indices'] = self.df_reads.apply(self.generate_random_fragment_indicies, axis=1)
        self.df_reads['read_lists'] = self.df_reads.apply(self.slice_fragments_into_reads, axis=1)

    def save_read_df(self, storage_pathway):
        cols_removed = ['fragment_indices', 'read_lists']
        self.df_reads.drop(cols_removed, axis=1, inplace=True)

        keys = ['reference_defline', 'primer_name', 'primer_scheme', 'amplicon_sequence']
        self.df_amplicon = self.df_amplicon[keys]
        df_merge = pd.merge(self.df_reads, self.df_amplicon, on=keys, how='outer')

        cols = [col for col in df_merge if col != 'amplicon_sequence']
        cols.append('amplicon_sequence')
        df_merge = df_merge[cols]

        df_merge.sort_values(by=['reference_defline', 'start'], inplace=True)
        df_merge.to_csv(os.path.join(storage_pathway, f"{self.propotion_file_name}_metadata.tsv"), index=False, sep='\t')

    def write_fastq_files(self, storage_pathway):
        fastq_r1 = os.path.join(storage_pathway, self.fastq_name + "_R1.fastq")
        fastq_r2 = os.path.join(storage_pathway, self.fastq_name + "_R2.fastq")
        #read_count_total = self.df_reads['total_frags'].sum()
        with open(fastq_r1, "w") as R1, open(fastq_r2, "w") as R2:
            read_count = 1
            #reference_defline	primer_scheme	primer_name
            for indx, row in self.df_reads.iterrows():
                defline = f"@{row['reference_defline']}:{row['primer_scheme']}:{row['primer_name']}"
                
                for forward, reverse in row['read_lists']:
                    
                    reverse = self.reverse_compliment(reverse)
                    defline_forward = defline + f"_R1/{read_count}"
                    defline_reverse = defline + f"_R2/{read_count}"
                    
                    quality_score_R1 = len(forward) * "I"
                    quality_score_R2 = len(reverse) * "I"

                    R1.write(f"{defline_forward}\n")
                    R1.write(f"{forward}\n")
                    R1.write(f"+\n")
                    R1.write(f"{quality_score_R1}\n")

                    R2.write(f"{defline_reverse}\n")
                    R2.write(f"{reverse}\n")
                    R2.write(f"+\n")
                    R2.write(f"{quality_score_R2}\n")
                    read_count +=1