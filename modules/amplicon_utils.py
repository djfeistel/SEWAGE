import os
import numpy as np
import pandas as pd
import textwrap

'''
Generate amplicions from a fasta dicitonary created in fasta_funcitons_class.py
'''
class GenerateAmplicons:
    
    def __init__(self,
                 reference_fasta_dict: dict,
                 scheme: str,
                 file_prefix_name: str="SEWAGE"
                 ):
        #
        self.reference_fasta_dict = reference_fasta_dict
        self.scheme = scheme
        self.file_prefix_name = file_prefix_name
        #
        self.df_amplicion = None
        self.primer_dict = None
        self.short_read_schemes = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vss1a", "vss2a", "vss2b"]
        self.long_read_schemes = ["vsl1a"]
        self.scheme_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "schemes")
        
    def get_df_amplicon(self):
        return self.df_amplicion.reset_index(drop=True)
        
    @staticmethod
    def reverse_complimentary_sequence(sequence: str) -> str:
        '''
        return reverse complimentary sequence for the reverse primer
        only used with ARTIC and VarSkip short reads primers is not 
        needed for VarSkip long reads
        '''
        compliment_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
        }
        complementary_sequence_list = [compliment_dict.get(nucl, nucl) for nucl in sequence]
        return "".join(complementary_sequence_list[::-1])
    
    def scheme_primer_dictionary(self) -> dict:
        '''return dictionary with primer name as key and F/R primer tuple as value'''
        scheme_pathway = os.path.join(self.scheme_dir, self.scheme + ".tsv")
        
        if "vsl1a" in scheme_pathway:
            # do not take reverse complimentary for varskip long reads scheme vsl1a
            df = pd.read_csv(scheme_pathway, sep='\t', header=0)
            df.set_index('primer_name', inplace=True)
            self.primer_dict = df[['forward_primer', 'reverse_primer']].apply(tuple, axis=1).to_dict()
        else:
            df = pd.read_csv(scheme_pathway, sep='\t', header=0)
            df['reverse_complimentary'] = df['reverse_primer'].apply(self.reverse_complimentary_sequence)
            df.set_index('primer_name', inplace=True)
            self.primer_dict = df[['forward_primer', 'reverse_complimentary']].apply(tuple, axis=1).to_dict()

    @staticmethod    
    def primer_short_read_index(reference_sequence: str, forward_primer: str, reverse_primer: str) -> list:
            '''
            find exact match to forward and reverse primer sequences
            if at least one primer sequence is not found, return False
            else determe the starting index for the forward primer and
            the ending index of the reverse primer
            return list of each index as an element
            '''
            forward_index_start = reference_sequence.find(forward_primer)
            reverse_index_end = reference_sequence.find(reverse_primer)
            if any(x == -1 for x in [forward_index_start, reverse_index_end]):
                return False
            # need to add length of primer for end of primer index
            reverse_index_end = reverse_index_end + len(reverse_primer)
            return [forward_index_start, reverse_index_end]
    
    @staticmethod
    def primer_long_read_index(reference_sequence: str, forward_primer: str,reverse_primer: str) -> list:
            '''
            R2 primer for Varskip does not need to be reverse complimentary
            varskip long-read reverse primer is sometimes upstream of forward
            if this occurs, return False
            '''
            forward_index_start = reference_sequence.find(forward_primer)
            reverse_index_end = reference_sequence.find(reverse_primer)
            if any(x == -1 for x in [forward_index_start, reverse_primer]):
                return False
            elif forward_index_start > reverse_index_end:
                return False
            else:
                reverse_index_end = reverse_index_end + len(reverse_primer)
                return [forward_index_start, reverse_index_end]

    def create_amplicon_dataframe(self) -> pd.DataFrame:
        '''find amplicons and return a pandas DataFrame'''
        if self.scheme in self.short_read_schemes:
            primer_indices_function = self.primer_short_read_index
        elif self.scheme in self.long_read_schemes:
            primer_indices_function = self.primer_long_read_index
        else:
            raise ValueError("Unknown scheme")
        
        # Preallocate DataFrame structure
        amplicon_data = []

        self.scheme_primer_dictionary()
        for reference_defline, reference_sequence in self.reference_fasta_dict.items():
            for primer_name, primer_sequences in self.primer_dict.items():
                forward_primer, reverse_primer = primer_sequences

                index_amplicon_list = primer_indices_function(reference_sequence=reference_sequence, 
                                                    forward_primer=forward_primer, 
                                                    reverse_primer=reverse_primer)
                if not index_amplicon_list:
                    start, end = pd.NA, pd.NA 
                    amplicon = "No Amplification"
                    amplicon_length = 0
                else:
                    start, end = [int(x) for x in index_amplicon_list]
                    amplicon = reference_sequence[start:end]
                    amplicon_length = int(len(amplicon))
                    # remove the primers from the ends of the amplicon 
                    # to avoid higher coverages at the end 
                    # perahps make that a flag?
                
                # Collect data in a tuple
                amplicon_data.append((reference_defline, self.scheme, primer_name, forward_primer, reverse_primer, start, end, amplicon_length, amplicon))

        # Create DataFrame from collected data
        column_names = ["reference_defline", "primer_scheme", "primer_name", "forward_primer", "reverse_primer", "start", "end", "amplicon_length_bp", "amplicon_sequence"]
        self.df_amplicion = pd.DataFrame(amplicon_data, columns=column_names).sort_values(by=['reference_defline', 'start'])

        # Optimize data types
        convert_cols_to_int = ['amplicon_length_bp']
        self.df_amplicion[convert_cols_to_int] = self.df_amplicion[convert_cols_to_int].astype('int64')
        
    def save_amplicon_meta_data(self, storage_pathway):
        metadata_file_name = f"{self.file_prefix_name}_amplicons.tsv"
        self.df_amplicion.to_csv(os.path.join(storage_pathway, metadata_file_name), sep='\t', index=False)
    
    def save_amplicon_fasta(self, storage_pathway):
        fasta_file_name = f"{self.file_prefix_name}_amplicons.fasta"

        # write fasta file
        with open(os.path.join(storage_pathway, fasta_file_name), 'w') as f:
            for idx, row in self.df_amplicion.iterrows():
                if row['amplicon_sequence'] == "No Amplification":
                    continue
                defline = f">{row['reference_defline']}~~~{row['primer_scheme']}~~~{row['primer_name']}~~~{row['start']}~~~{row['end']}~~~{row['amplicon_length_bp']}bp"
                f.write(defline + '\n')
                # Wrap the sequence
                sequence = row['amplicon_sequence']
                wrapped_sequence = textwrap.wrap(sequence, width=60)
                for line in wrapped_sequence:
                    f.write(line + '\n')
