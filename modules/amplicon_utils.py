import os
import numpy as np
import pandas as pd
from datetime import datetime
#remove after finsiehd 
#from fasta_funcitons_class_workshop import load_fasta_workflow_test
'''
Generate amplicions from a fasta dicitonary created in fasta_funcitons_class.py
'''
class GenerateAmplicons:
    
    def __init__(self,
                 fasta_dict: dict,
                 scheme: str,
                 amplicon_fasta_name=None,
                 amplicon_storage_dir:str=None
                 ):
        
        self.fasta_dict = fasta_dict
        self.scheme = scheme
        self.scheme_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "schemes")
        self.amplicon_fasta_name = amplicon_fasta_name
        self.amplicon_storage_dir = amplicon_storage_dir
        self.date_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.short_read_schemes = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vss1a", "vss2a", "vss2b"]
        self.long_read_schemes = ["vsl1a"]
        
    def check_amplicon_storage_dir(self):
        # Set directory name based on whether self.amplicon_storage_dir is None
        current_datetime = self.date_time
        dir_name = self.amplicon_storage_dir if self.amplicon_storage_dir is not None else f"SEWAGE_{current_datetime}"
        
        try:
            os.mkdir(dir_name)
            #print(f"Directory '{dir_name}' for amplicon stroage was created successfully.")
        except FileExistsError:
            raise FileExistsError(f"Directory '{dir_name}' for amplicon stroage already exists.")
        except FileNotFoundError:
            raise FileNotFoundError()
        except Exception as e:
            print(f"An error occurred: {e}")
            raise e
        self.amplicon_storage_dir = os.path.realpath(dir_name)
        return self.amplicon_storage_dir
        
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
        
        scheme_pathway = os.path.join(self.scheme_dir, self.scheme + ".tsv")
        
        if "vsl1a" not in scheme_pathway:
            df = pd.read_csv(scheme_pathway, sep='\t', header=0)
            df['reverse_complimentary'] = df['reverse_primer'].apply(GenerateAmplicons.reverse_complimentary_sequence)
            df.set_index('primer_name', inplace=True)
            primer_scheme_dict = df[['forward_primer', 'reverse_complimentary']].apply(tuple, axis=1).to_dict()
            # primer_scheme_dict = Key: primer name; value: (forward primer, reverse primer)
            return primer_scheme_dict
        else:
            # do not take reverse complimentary for varskip long reads
            df = pd.read_csv(scheme_pathway, sep='\t', header=0)
            df.set_index('primer_name', inplace=True)
            primer_scheme_dict = df[['forward_primer', 'reverse_primer']].apply(tuple, axis=1).to_dict()
            return primer_scheme_dict

    @staticmethod    
    def primer_short_read_index(
            reference_sequence: str, 
            forward_primer: str, 
            reverse_primer: str
            ) -> list:
            '''
            find exact match to forward and reverse primer sequences
            if at least one primer sequence is not found, return False
            else determing the starting index for the forward primer and
                the ending index of the reverse primer
            return a list of each index as an element
            '''
            forward_index_start = reference_sequence.find(forward_primer)
            reverse_index_end = reference_sequence.find(reverse_primer)
            if any(x == -1 for x in [forward_index_start, reverse_index_end]):
                return False
            # need to add length of primer for end of primer index
            reverse_index_end = reverse_index_end + len(reverse_primer)
            return [forward_index_start, reverse_index_end]
    
    @staticmethod
    def primer_long_read_index(
            reference_sequence: str, 
            forward_primer: str,
            reverse_primer: str

            ) -> list:
            forward_index_start = reference_sequence.find(forward_primer)
            # R2 primer for Varskip does not need to be reverse complimentary
            reverse_index_end = reference_sequence.find(reverse_primer)
            if any(x == -1 for x in [forward_index_start, reverse_primer]):
                return False
            # for some reason, the reverse_index_end is upstream (i.e., smaller index) 
            # of the forward_index_start when searching for long reads.
            # for now thiselogic removes those when this happens 
            # and calls them as No Amplification
            elif forward_index_start > reverse_index_end:
                return False
            else:
                # need to add length of primer for end of primer index
                reverse_index_end = reverse_index_end + len(reverse_primer)
                return [forward_index_start, reverse_index_end]

    def amplicon_dataframe(
            self, 
            genome_sequence_dict: dict, 
            primer_scheme_dict: dict
    ) -> pd.DataFrame:
        '''find amplicions and return a pandas df'''
        if self.scheme in self.short_read_schemes:
            primer_indices_function = GenerateAmplicons.primer_short_read_index
        elif self.scheme in self.long_read_schemes:
            primer_indices_function = GenerateAmplicons.primer_long_read_index
        else:
            raise ValueError("Unknown scheme")
        
        amplicon_dict = {
            "reference_defline": [],
            "primer_scheme": [],
            "primer_name": [],
            "start": [],
            "end": [],
            "length_bp": [],
            "amplicon_sequence": []
        }
            
        for reference_defline, reference_sequence in genome_sequence_dict.items():
            for primer_name, primer_sequences in primer_scheme_dict.items():
                forward_primer = primer_sequences[0]
                reverse_primer = primer_sequences[1]

                index_amplicon_list = primer_indices_function(reference_sequence=reference_sequence, 
                                                      forward_primer=forward_primer, 
                                                      reverse_primer=reverse_primer)
                if not index_amplicon_list:
                    start, end = np.nan, np.nan 
                    amplicion = "No Amplification"
                    amplicion_length = np.nan
                else:
                    start, end = [int(x) for x in index_amplicon_list]
                    amplicion = reference_sequence[start:end]
                    amplicion_length = int(len(amplicion))
                
                amplicon_dict["reference_defline"].append(reference_defline)
                amplicon_dict["primer_scheme"].append(self.scheme)
                amplicon_dict["primer_name"].append(primer_name)
                amplicon_dict["start"].append(start)
                amplicon_dict["end"].append(end)
                amplicon_dict["length_bp"].append(amplicion_length)
                amplicon_dict["amplicon_sequence"].append(amplicion)
        
        df = pd.DataFrame(amplicon_dict).sort_values(by=['reference_defline', 'start'])
        convert_cols_to_int = ['start', 'end', 'length_bp']
        df[convert_cols_to_int] = df[convert_cols_to_int].apply(lambda x: x.astype('Int64'))
        return df
    
    def write_fasta(self, df):
        '''write fasta file with amplicions'''

        # chcek if none and assign name
        if self.amplicon_fasta_name is None:
            fasta_file_name = f"SEWAGE_amplicons.fasta"
            metadata_file_name = f"SEWAGE_amplicons_metadata.tsv"
        else:
            fasta_file_name = self.amplicon_fasta_name + ".fasta"
            metadata_file_name = self.amplicon_fasta_name + "_metadata.tsv"

        #save original amplicon file
        df.to_csv(os.path.join(self.amplicon_storage_dir, metadata_file_name), sep='\t', index=False)

        # Create the defline using vectorized string concatenation
        df = df.dropna(subset=['start', 'end', 'length_bp'], how='all')
        deflines = ('>' + df['reference_defline'].astype(str) + '~~~' +
                    df['primer_scheme'].astype(str) + '~~~' +
                    df['primer_name'].astype(str) + '~~~' +
                    df['start'].astype(str) + '~~~' +
                    df['end'].astype(str) + '~~~' +
                    df['length_bp'].astype(str) + 'bp')

        # Concatenate deflines and sequences with newline characters
        fasta_strings = deflines + '\n' + df['amplicon_sequence'] + '\n'

        # Write all lines to the file at once
        with open(os.path.join(self.amplicon_storage_dir, fasta_file_name), 'w') as f:
            f.write('\n'.join(fasta_strings))
        return df