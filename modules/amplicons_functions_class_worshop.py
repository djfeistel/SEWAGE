import os
import sys
import argparse
import pandas as pd
from datetime import datetime
#remove after finsiehd 
from fasta_funcitons_class_workshop import load_fasta_workflow_test
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
        self.amplicon_fasta_name = amplicon_fasta_name
        self.amplicon_storage_dir = amplicon_storage_dir
        self.scheme_dir = "../schemes"
        self.date_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.short_read_schemes = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vss1a", "vss2a", "vss2b"]
        self.long_read_schemes = ["vsl1a"]
        # if not os.path.exists(self.fasta_pathway):
        #     raise FileNotFoundError(f"FASTA pathway '{self.fasta_pathway}' does not exist.")
        # if not os.path.isdir(self.fasta_pathway):
        #     raise FileNotFoundError(f"FASTA pathway '{self.fasta_pathway}' is not a directory.")
        # fasta_files = [f for f in os.listdir(self.fasta_pathway) if os.path.isfile(os.path.join(self.fasta_pathway, f)) and f.lower().endswith(('.fasta', '.fsa', '.fa'))]
        # if not fasta_files:
        #     raise FileNotFoundError(f"No fasta files with allowed extensions [.fasta, .fsa, .fa] found in the directory '{self.fasta_pathway}'.")
        # if self.scheme is None and self.user is None:
        #     raise TypeError(f"Must choose at least one: --scheme or --user")
        # if not os.path.exists(self.sewage_dir):
        #     raise FileNotFoundError(f"Sewage directory '{self.sewage_dir}' does not exist.")
        # if not os.path.isdir(self.sewage_dir):
        #     raise ValueError(f"Sewage directory '{self.sewage_dir}' is not a directory.")
        # if not os.path.exists(self.amplicon_storage_pathway):
        #     raise FileNotFoundError(f"Amplicon storage pathway '{self.amplicon_storage_pathway}' does not exist.")
        # if not os.path.isdir(self.amplicon_storage_pathway):
        #     raise ValueError(f"Amplicon storage pathway '{self.amplicon_storage_pathway}' is not a directory.")
    

    def check_amplicon_storage_dir(self):
        # Set directory name based on whether self.amplicon_storage_dir is None
        current_datetime = self.date_time
        dir_name = self.amplicon_storage_dir if self.amplicon_storage_dir is not None else f"SEWAGE_amplicons_{current_datetime}"
        
        try:
            os.mkdir(dir_name)
            print(f"Directory '{dir_name}' for amplicon stroage was created successfully.")
        except FileExistsError:
            raise FileExistsError(f"Directory '{dir_name}' for amplicon stroage already exists.")
        except FileNotFoundError:
            raise FileNotFoundError()
        except Exception as e:
            print(f"An error occurred: {e}")
            raise e
        self.amplicon_storage_dir = os.path.realpath(dir_name)
        
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
            df = pd.read_csv(scheme_pathway, sep='\t', header=None)
            df[3] = df[2].apply(GenerateAmplicons.reverse_complimentary_sequence)
            df.set_index(0, inplace=True)
            primer_scheme_dict = df[[1, 3]].apply(tuple, axis=1).to_dict()
            # primer_scheme_dict = Key: primer name; value: (forward primer, reverse primer)
            return primer_scheme_dict
        else:
            df = pd.read_csv(scheme_pathway, sep='\t', header=None)
            df.set_index(0, inplace=True)
            primer_scheme_dict = df[1].to_dict()
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
                    start, end = None, None
                    amplicion = "No Amplification"
                else:
                    start, end = index_amplicon_list
                    amplicion = reference_sequence[start:end]
                
                amplicon_dict["reference_defline"].append(reference_defline)
                amplicon_dict["primer_scheme"].append(self.scheme)
                amplicon_dict["primer_name"].append(primer_name)
                amplicon_dict["start"].append(start)
                amplicon_dict["end"].append(end)
                amplicon_dict["length_bp"].append(len(amplicion))
                amplicon_dict["amplicon_sequence"].append(amplicion)
        
        return pd.DataFrame(amplicon_dict)
    
    def write_fasta(self, df):
        '''write fasta file with amplicions'''
        if self.amplicon_fasta_name is None:
            file_name = f"SEWAGE_amplicon_{self.date_time}.fasta"
        else:
            file_name = self.amplicon_fasta_name
        # Create the defline using vectorized string concatenation
        df = df[df['amplicon_sequence'] != "No Amplification"]
        deflines = ('>' + df['reference_defline'].astype(str) + '~~~' +
                    df['primer_scheme'].astype(str) + '~~~' +
                    df['primer_name'].astype(str) + '~~~' +
                    df['start'].astype(str) + '~~~' +
                    df['end'].astype(str) + '~~~' +
                    df['length_bp'].astype(str) + 'bp')

        # Concatenate deflines and sequences with newline characters
        fasta_strings = deflines + '\n' + df['amplicon_sequence'] + '\n'

        # Write all lines to the file at once
        with open(os.path.join(self.amplicon_storage_dir, file_name), 'w') as f:
            f.write('\n'.join(fasta_strings))

def opts():
    # Create the main parser
    amplicon_parser = argparse.ArgumentParser(prog="SEWAGE amplicon",
                                     description=f"SEWAGE NEEDS A NAME\n",
                                     add_help=True,
                                     epilog="minimal usage: SEWAGE amplicon -f <input> -s <scheme>")


    amplicon_parser.add_argument('-f', '--fasta', 
                                 help='Pathway to directory with single reference genome fasta files [.fasta, .fsa, .fa]', 
                                 required=True,
                                 type=str,
                                 dest='fasta')
    amplicon_parser.add_argument('-s', '--scheme', 
                                 help='Primer scheme: (Artic = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2"], \
                                    VarSkip = ["vsl1a", "vss1a", "vss2a", "vss2b"]) [default=V5.3.2]', 
                                 required=False, 
                                 type=str,
                                 choices=["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vsl1a", "vss1a", "vss2a", "vss2b"],
                                 metavar='SCHEME',
                                 dest='scheme',
                                 default=None) #make sure that code reflect when None is reached it does not work
    amplicon_parser.add_argument('-u', '--user_scheme', 
                                 help='User defined primer scheme TSV file with three columns with: Forward_seq, Reverse_seq, Primer_Name', 
                                 required=False, 
                                 type=str,
                                 metavar='FILE',
                                 dest='user',
                                 default=None) #change this later
    amplicon_parser.add_argument('-o', '--output', 
                                 help='Output directory name for amplicon [default="SEWAGE_amplicons"]', 
                                 required=False,
                                 default=None,
                                 type=str,
                                 metavar='STR',
                                 dest='output')
    amplicon_parser.add_argument('-p', '--pathway', 
                                 help='Pathway to storgage directory [default="."]',
                                 required=False, 
                                 type=str,
                                 metavar='PATHWAY',
                                 dest='pathway',
                                 default='.')
    args = amplicon_parser.parse_args()    
    return args

def generate_amplicon_workflow_test(fasta_input):

    fasta_dict = load_fasta_workflow_test(fasta_input)
    amplicon_generator = GenerateAmplicons(fasta_dict=fasta_dict, scheme="V5.3.2")
    #amplicon_generator.check_amplicon_storage_dir()
    primer_scheme_dict = amplicon_generator.scheme_primer_dictionary()
    df_amplicon = amplicon_generator.amplicon_dataframe(
        genome_sequence_dict=fasta_dict,
        primer_scheme_dict=primer_scheme_dict
        )
    #amplicon_generator.write_fasta(df=df_amplicon)
    #for whatever reason i need to filter No Amplified amplicons here
    df_amplicon = df_amplicon[df_amplicon['amplicon_sequence'] != "No Amplification"]
    return df_amplicon

if __name__ == "__main__":
    generate_amplicon_workflow_test()