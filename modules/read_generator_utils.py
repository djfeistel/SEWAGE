import numpy as np
import pandas as pd
import os
import sys
import time
#from modules.quality_scores_util_workshop import quailty_scores_generator_array

class ReadGenerator():
    def __init__(
            self,
            df_amplicon:pd.DataFrame,
            fastq_name:str="SEWAGE",
            file_prefix_name:str="SEWAGE",
            read_length:int=150,
            frag_length:int=300,
            coverage_depth:int=500,
            seed:int=13,
    ):
        self.df_amplicon = df_amplicon.reset_index(drop=True)
        self.fastq_name = fastq_name
        self.file_prefix_name = file_prefix_name
        self.read_length = np.int64(read_length)
        self.frag_length = np.int64(frag_length)
        self.coverage_depth = np.int64(coverage_depth)
        self.seed = seed
        
    def get_df_reads(self):
        return self.df_reads

    def get_df_amplicons(self):
        return self.df_amplicon
    
    def calculate_total_frags_to_cover_amplicon(self):
        total_frags = (self.coverage_depth * self.df_amplicon['amplicon_length_bp'] // self.frag_length)
        self.df_amplicon['total_frags'] = (self.df_amplicon['proportion'] * total_frags).astype("int64")

    def add_amplicon_number(self):
        self.df_amplicon['amplicon_number'] = range(1, len(self.df_amplicon) + 1)    

    # def create_df_reads(self):
    #     '''create df_reads'''
    #     repeats = self.df_amplicon['total_frags'].values
    #     self.df_reads = pd.DataFrame(np.repeat(self.df_amplicon.values, repeats, axis=0), columns=self.df_amplicon.columns)

    def create_df_reads(self):
        '''create df_reads in a more efficient way.'''
        # Instead of repeating rows, create a new DataFrame with the necessary number of rows
        # and use a more efficient method to populate it.
        total_rows = self.df_amplicon['total_frags'].sum()
        self.df_reads = pd.DataFrame(index=range(total_rows), columns=self.df_amplicon.columns)

        # Populate the new DataFrame efficiently.
        start_idx = 0
        for _, row in self.df_amplicon.iterrows():
            end_idx = start_idx + row['total_frags']
            self.df_reads.iloc[start_idx:end_idx] = row
            start_idx = end_idx


    def drop_irrelevant_columns_from_df_reads(self):
        drop_columns = ['forward_primer', 'reverse_primer', 'start', 'end', 'proportion']
        self.df_reads.drop(columns=drop_columns, inplace=True) 
    
    def add_read_number(self):
        self.df_reads['read_number'] = range(1, len(self.df_reads) + 1)    

    def generate_random_fragment_indicies_new(self, row):
        np.random.seed(None)
        amplicon_length_bp = row['amplicon_length_bp']
        start_inx = np.random.randint(low=0, high=(amplicon_length_bp - self.frag_length), size=1)
        end_inx = start_inx + self.frag_length
        return int(start_inx), int(end_inx)
    
    def slice_fragments_into_reads(self, row):
        start = row['fragment_start']
        end = row['fragment_end']
        amplicon_sequence = row['amplicon_sequence']
        
        # Adjusted to handle a list containing one tuple
        R1 = amplicon_sequence[start:end][:self.read_length]
        R2 = amplicon_sequence[start:end][self.read_length:]
    
        return R1, R2
    
    def quailty_scores_generator_array(
        self,
        max_q:int = 40,
        min_q:int = 20,
        std_dev=3,
        begining_q:int=4,
        begining_bp:int=10
    ):
        '''workshoping this one
        currently creates a linear decrease in qvalues'''
        
        def map_values_to_chars(array_2d):
            illumina_qscore_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"
            #max_index = 40  # max quality score
            lookup_array = np.array(list(illumina_qscore_string))

            # Vectorized mapping
            # element in array acts as index for lookup
            # more effienct
            mapped_array = lookup_array[array_2d.astype(int)]
            joined_arrays = [''.join(map(str, sub_array)) for sub_array in mapped_array]
            return joined_arrays
        
        total_reads = int(len(self.df_reads))*2
        half_reads = int(total_reads/2)
        array_2d = np.tile(np.linspace(max_q ,min_q, self.read_length), (total_reads, 1))
        # create SD for each position
        flux_2d = np.random.normal(0, std_dev, (total_reads, self.read_length))
        # add the SD to original array
        final_array = array_2d + flux_2d
        # clip values outside teh rang eof quality scores
        final_array = np.clip(final_array, 0, 40) # qscores are between 0-40
        # round values
        final_array = np.round(final_array)
        # mimic illumina data where first ~10bp are slighly less quality 
        final_array[:, :begining_bp] -= begining_q
        #final_array[num_reads:] -= 1
        final_array = map_values_to_chars(final_array)
        
        R1_q = final_array[:half_reads]
        R2_q = final_array[half_reads:]

        return R1_q, R2_q

    @staticmethod    
    def reverse_compliment(reverse_strand:str):
        compliment = {"A": "T", "T": "A", "G": "C", "C": "G"}
        complement_strand_list = [compliment[base] if base in compliment else 'N' for base in reverse_strand]
        complement_strand = ''.join(complement_strand_list)
        return complement_strand[::-1]
    
    def create_reads_workflow(self):
        start_time = time.time()
        self.calculate_total_frags_to_cover_amplicon()
        self.df_amplicon.dropna(how='any', inplace=True)
        self.add_amplicon_number()
        self.create_df_reads()
        self.drop_irrelevant_columns_from_df_reads()
        self.add_read_number() # add in amplicon section later
        self.df_reads[['fragment_start', 'fragment_end']] = self.df_reads.apply(self.generate_random_fragment_indicies_new, axis=1, result_type='expand')
        self.df_reads[['R1_read', 'R2_read']] = self.df_reads.apply(self.slice_fragments_into_reads, axis=1, result_type='expand')
        R1_q, R2_q = self.quailty_scores_generator_array()
        self.df_reads['R1_q'] = R1_q
        self.df_reads['R2_q'] = R2_q

        end_time = time.time() - start_time
        print(end_time)
        
        
    def save_read_df(self, storage_pathway):
        self.df_reads.to_csv(os.path.join(storage_pathway, f"{self.file_prefix_name}_reads.tsv"), index=False, sep='\t')

    
    def write_fastq_files(self, storage_pathway):
        fastq_r1 = os.path.join(storage_pathway, self.fastq_name + "_R1.fastq")
        fastq_r2 = os.path.join(storage_pathway, self.fastq_name + "_R2.fastq")
        
        with open(fastq_r1, "w") as R1, open(fastq_r2, "w") as R2:
            read_count = 1
            #reference_defline	primer_scheme	primer_name
            for indx, row in self.df_reads.iterrows():
                defline_R1 = f"@{row['reference_defline']}:{row['primer_scheme']}:{row['primer_name']}:A{row['amplicon_number']}:R1/{row['read_number']}"
                defline_R2 = f"@{row['reference_defline']}:{row['primer_scheme']}:{row['primer_name']}:A{row['amplicon_number']}:R2/{row['read_number']}"
                R1_read = row['R1_read']
                R2_read = self.reverse_compliment(row['R2_read'])
                R1_q = row['R1_q']
                R2_q = row['R2_q']
                
                R1.write(f"{defline_R1}\n")
                R1.write(f"{R1_read}\n")
                R1.write(f"+\n")
                R1.write(f"{R1_q}\n")

                R2.write(f"{defline_R2}\n")
                R2.write(f"{R2_read}\n")
                R2.write(f"+\n")
                R2.write(f"{R2_q}\n")
                