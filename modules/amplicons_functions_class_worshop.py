import os
import sys
import argparse
import panda as pd

class GenerateAmplicons:
    
    def __init__(self,
                 fasta_dict: dict,
                 scheme: str,
                 amplicon_storage_dir:str=None
                 ):
        
        self.fasta_dict = fasta_dict
        self.scheme = scheme
        self.amplicon_storage_dir = amplicon_storage_dir
        self.scheme_dir = "../schemes"
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
        if self.amplicon_storage_dir is not None:
            if not os.path.isdir(self.amplicon_storage_dir):
                try:
                    os.mkdir(self.amplicon_storage_dir)
                except Exception as e:
                    raise e
            return self.amplicon_storage_dir
        else:
            amplicon_dir = "SEWAGE_amplicons"
            if not os.path.isdir(amplicon_dir):
                try:
                    os.mkdir(amplicon_dir)
                except Exception as e:
                    raise e
            return amplicon_dir

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
    
    def scheme_dictionary(self) -> dict:

        scheme_pathway = os.path.join(self.scheme_dir, self.scheme + ".tsv")
        if "vsl1a" not in scheme_pathway:
            df = pd.read_csv(scheme_pathway, sep='\t', header=None)
            df[3] = df[2].apply(reverse_complimentary_sequence)
            df.set_index(0, inplace=True)
            primer_scheme_dict = df[[1, 3]].apply(tuple, axis=1).to_dict()
            # primer_scheme_dict = Key: primer name; value: (forward primer, reverse primer)
            return primer_scheme_dict
        else:
            df = pd.read_csv(scheme_pathway, sep='\t', header=None)
            df.set_index(0, inplace=True)
            primer_scheme_dict = df[1].to_dict()
            return primer_scheme_dict
        
    def primer_paired_end_indices(
            self,
            reference_sequence: str, 
            forward_primer: str, 
            reverse_compliment_primer: str
            ) -> list:
            '''
            find exact match to forward and reverse primer sequences
            if at least one primer sequence is not found, return False
            else determing the starting index for the forward primer and
                the ending index of the reverse primer
            return a list of each index as an element
            '''
            forward_index = reference_sequence.find(forward_primer)
            reverse_index_start = reference_sequence.find(reverse_compliment_primer)
            if any(x == -1 for x in [forward_index, reverse_index_start]):
                return False
            reverse_index_end = reverse_index_start + len(reverse_compliment_primer)
            return [forward_index, reverse_index_end]
    
    def primer_long_read_indices(
            self,
            reference_sequence: str, 
            forward_primer: str
            ) -> list:
            start = reference_sequence.find(forward_primer)
            if start == -1:
                return False
            else:
                end = start + len(start)
                return [start, end]
    
    def create_amplicon_dict_for_fasta_files(self,
                                            ref_seq_dict: dict,
                                            primer_scheme_dict: dict) -> dict:
        '''
        requires a for loop for calling reference_sequence dict key : values from read_referecne_fasta_list()
        '''
        short_read_schemes = [
            "V1", "V2", "V3", "V4", "V4.1", 
            "V5.3.2", "vss1a", "vss2a", "vss2b"
            ]
        long_read_schemes = ["vsl1a"]
        amplicon_dict = {}
        #key==primer name; value==amplicon or None
        for reference_sequence_defline, reference_sequence_values in ref_seq_dict.items():
            for primer_name, primer_seqs in primer_scheme_dict.items():
                    
                if type(primer_seqs) == tuple:
                    forward_primer = primer_seqs[0]
                    reverse_primer = primer_seqs[1]
                    index_amplicon_list = primer_paired_end_indices(reference_sequence=reference_sequence_values, forward_primer=forward_primer, reverse_primer=reverse_primer
                    )
                else:
                    index_amplicon_list = primer_long_read_indices(reference_sequence=primer_seqs, forward_primer=primer_seqs)

'''             STOPED HERE             '''
                
                '''amplicon dict is 
                amplicon_dict[reference_sequence_name] = {defline: sequence}'''

                if not index_amplicon_list:
                    defline = "~~~".join([primer_name, reference_sequence_defline, scheme_profile])
                    if reference_fasta_basename in amplicon_dict:
                        amplicon_dict[reference_fasta_basename][defline] = None
                    else:
                        amplicon_dict[reference_fasta_basename] = {defline: None}
                else:
                    amplicion_seq = str(reference_sequence_seq[index_amplicon_list[0]:index_amplicon_list[1]])
                    amplicon_length = str(len(amplicion_seq)) + "bp"
                    defline = "~~~".join([primer_name, reference_sequence_defline, scheme_profile, amplicon_length])
                    if reference_fasta_basename in amplicon_dict:
                        amplicon_dict[reference_fasta_basename][defline] = [amplicion_seq, index_amplicon_list]
                    else:
                        amplicon_dict[reference_fasta_basename] = {defline: [amplicion_seq, index_amplicon_list]}
        return amplicon_dict
    
    def write_amplicion_fasta_files(self, amplicon_storage_dir_pathway: str, amplicon_dict: dict) -> None:

        if self.scheme and not self.user:
            scheme_profile = self.scheme
        elif self.user and not self.scheme:
            scheme_profile = self.user

        for ref_name, seq_dict in amplicon_dict.items():

            output_fasta = os.path.join(amplicon_storage_dir_pathway, ref_name[0] + "_amplicons.fasta")
            output_log = os.path.join(amplicon_storage_dir_pathway, ref_name[0] + "_amplicons.log")
            
            with open(output_fasta, 'w') as wf_amps, open(output_log, 'w') as wf_logs:
            
                for defline, amplicon_items_list in seq_dict.items():
                    primer_name = defline.split('~~~')[0]

                    if amplicon_items_list is None:
                        wf_logs.write(f"{primer_name}\t{scheme_profile}\t{ref_name[0]}\tNo Amplification\n")

                    elif amplicon_items_list is not None:
                        amplicon_sequence, position_list = amplicon_items_list
                        start_pos, end_pos = position_list

                        line = '\t'.join([primer_name, scheme_profile, ref_name[0], str(start_pos), str(end_pos), str(len(amplicon_sequence))])
                        wf_logs.write(f"{line}\n")

                        wf_amps.write(f">{defline}\n")
                        while amplicon_sequence:
                            chunk = amplicon_sequence[:100]
                            amplicon_sequence = amplicon_sequence[100:]
                            wf_amps.write(f"{chunk}\n")
        return None

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

def GenerateAmplicon_MainWorkflow(fasta, scheme, user, sewage_home_dir, output, pathway):
    '''Main workflow for creating amplicons'''
    '''removed the opts as the main_argsparser is used'''

    sewage_amplicon = GenerateAmplicons(fasta, 
                                        scheme,
                                        user, 
                                        sewage_home_dir,
                                        output, 
                                        pathway)
    
    fasta_pathway_list = sewage_amplicon.fasta_reference_pathway_list()
    reference_sequence_dict = sewage_amplicon.read_referecne_pathway_fasta_list(fasta_pathway_list)
    primer_scheme_dict = sewage_amplicon.read_primer_scheme_file()
    amplicon_dict = sewage_amplicon.create_amplicon_dict_for_fasta_files(ref_seq_dict=reference_sequence_dict,
                                                                         primer_scheme_dict=primer_scheme_dict)
    amplicon_storage_dir_pathway = sewage_amplicon.make_amplicon_storage_dir()
    sewage_amplicon.write_amplicion_fasta_files(amplicon_storage_dir_pathway=amplicon_storage_dir_pathway,
                                                amplicon_dict=amplicon_dict)
    
    print(f"\nSEWAGE has finished creating amplicon fasta files.", file=sys.stderr)
    print(f"Data is stored in {amplicon_storage_dir_pathway}.\n", file=sys.stderr)
    print(f"Use 'SEWAGE enrich' command with any set of amplicon fasta files to generate mix proporitons.\n")

    return