import os
import sys
import glob
import argparse

'''some Class functions are not use but should be updated'''
class GenerateAmplicons:
    
    def __init__(self: str,
                 fasta_pathway: str,
                 scheme: str,
                 sewage_dir: str,
                 amplicon_storage_dir: str,
                 amplicon_storage_pathway: str
                 ):
        self.fasta_pathway = os.path.realpath(fasta_pathway)
        self.scheme = scheme
        self.sewage_dir = os.path.realpath(sewage_dir)
        self.amplicon_storage_dir = amplicon_storage_dir
        self.amplicon_storage_pathway = os.path.realpath(amplicon_storage_pathway)

        if not os.path.exists(self.fasta_pathway):
            raise FileNotFoundError(f"FASTA pathway '{self.fasta_pathway}' does not exist.")
        if not os.path.isdir(self.fasta_pathway):
            raise FileNotFoundError(f"FASTA pathway '{self.fasta_pathway}' is not a directory.")
        fasta_files = [f for f in os.listdir(self.fasta_pathway) if os.path.isfile(os.path.join(self.fasta_pathway, f)) and f.lower().endswith(('.fasta', '.fsa', '.fa'))]
        if not fasta_files:
            raise FileNotFoundError(f"No fasta files with allowed extensions [.fasta, .fsa, .fa] found in the directory '{self.fasta_pathway}'.")
        if not os.path.exists(self.sewage_dir):
            raise FileNotFoundError(f"Sewage directory '{self.sewage_dir}' does not exist.")
        if not os.path.isdir(self.sewage_dir):
            raise ValueError(f"Sewage directory '{self.sewage_dir}' is not a directory.")
        if not os.path.exists(self.amplicon_storage_pathway):
            raise FileNotFoundError(f"Amplicon storage pathway '{self.amplicon_storage_pathway}' does not exist.")
        if not os.path.isdir(self.amplicon_storage_pathway):
            raise ValueError(f"Amplicon storage pathway '{self.amplicon_storage_pathway}' is not a directory.")
    
    def make_amplicon_storage_dir(self):
        if self.amplicon_storage_dir is None:
            pathway = os.path.join(self.amplicon_storage_pathway, "SEWAGE_amplicons")
        else:
            pathway = os.path.join(self.amplicon_storage_pathway, self.amplicon_storage_dir)
        if not os.path.isdir(pathway):
            os.mkdir(pathway)
        return pathway
    
    def Check_fasta_reference_input(self) -> bool:
        '''is pathway: True, if file, False, else None'''
        if os.path.isdir(self.fasta_pathway):
            return True
        elif os.path.isfile(self.fasta_pathway):
            return False
        else:
            return None
    
    def fasta_reference_pathway_list(self) -> list:
        '''
        return list of tuples
        [(fasta pathway, fasta file), (fasta pathway, fasta file), ...]
        '''
        file_pathway_list = glob.glob(self.fasta_pathway + "/*")
        valid_extensions = ('.fasta', '.fa', '.fsa')
        fasta_file_pathway_list = [(file, os.path.basename(file)) for file in file_pathway_list if os.path.splitext(file)[1].lower() in valid_extensions]
        if len(fasta_file_pathway_list) == 0:
            print(f"\nPathway to fasta directory did not find fasta files.", file=sys.stderr)
            print(f"Check pathway and that valid fasta files are accepted [.fasta, .fsa, .fa]", file=sys.stderr)
            sys.exit(0)
        else:
            return fasta_file_pathway_list
        
    def fasta_reference_file_to_list(self) -> list:
        '''
        return tuple as [(fasta pathway)]
        '''
        return [(self.fasta_pathway, os.path.basename(self.fasta_pathway))]
    
    def read_referecne_pathway_fasta_list(self, fasta_pathway_list: list) -> dict:
        '''
        fasta_pathway_list is a list with at least one pathway to a fasta file
        each fasta pathway can be wither single of multifasta 
        '''
        sequences = {}
        current_defline = ""  # Initialize with an empty string
        current_sequence = ""
        for fasta_pathway, fasta_file in fasta_pathway_list:
            with open(fasta_pathway, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue  # Skip empty lines
                    if line.startswith(">"):  # Def line
                        if current_defline:
                            sequences[current_defline] = [current_sequence, fasta_file.split('.')[0]]
                        current_defline = line[1:] #remove > 
                        current_sequence = ""
                    else:
                        current_sequence += line
                # Add the last sequence after reaching the end of the file
                if current_defline and current_sequence:
                    sequences[current_defline] = [current_sequence, fasta_file.split('.')[0]]
            current_defline = ""
            current_sequence = ""
        return sequences

    def read_primer_scheme_file(self) -> dict:
        '''
        open custom primer file in scheme dir
        return dict with primer name as key and values as [forward, reverse] primer seqs
        '''
        scheme_pathway = os.path.join(os.path.join(self.sewage_dir, "schemes"), self.scheme + ".tsv")
        primer_scheme_dict = {}
        with open(scheme_pathway, 'r') as fh:
            for line in fh:
                # lines starting with # indicate information about the scheme
                if line.startswith("#"): continue
                line = line.strip().split('\t')
                primer_scheme_dict[line[0]] = line[1:]
        return primer_scheme_dict

    def create_amplicon_dict_for_fasta_files(self,
                                            ref_seq_dict: dict,
                                            primer_scheme_dict: dict) -> dict:
        '''
        requires a for loop for calling reference_sequence dict key : values from read_referecne_fasta_list()
        '''
        def reverse_complimentary_sequence(sequence: str) -> str:
            '''
            return reverse complimentary sequence for the reverse primer
            only used with ARTIC and VarSkip short reads primers
            is not needed for VarSkip long reads
            '''
            compliment_dict = {
                "A": "T",
                "T": "A",
                "G": "C",
                "C": "G",
            }
            complementary_sequence_list = [compliment_dict.get(nucl, nucl) for nucl in sequence]
            return "".join(complementary_sequence_list[::-1])

        def primer_index_for_amplicon_seq(reference_sequence: str, 
                                        forward_primer: str, 
                                        reverse_primer: str) -> list:
            '''
            find exact match to forward and reverse primer sequences
            if at least one primer sequence is not found, return False
            else determing the starting index for the forward primer and
                the ending index of the reverse primer
            return a list of each index as an element
            '''
            forward_index = reference_sequence.find(forward_primer)
            reverse_index_start = reference_sequence.find(reverse_primer)
            if any(x == -1 for x in [forward_index, reverse_index_start]):
                return False
            reverse_index_end = reverse_index_start + len(reverse_primer)
            return [forward_index, reverse_index_end]
        
        amplicon_dict = {} 
        #key==primer name; value==amplicon or None
        for reference_sequence_defline, reference_sequence_values in ref_seq_dict.items():
            reference_sequence_seq, reference_fasta_basename = reference_sequence_values

            for primer_name, primer_seqs in primer_scheme_dict.items():

                short_read_schemes = ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vss1a", "vss2a", "vss2b"]
                long_read_schemes = ["vsl1a"]
                
                if any(True for x in short_read_schemes if self.scheme == x):
                    forward_primer = primer_seqs[0]
                    reverse_primer = primer_seqs[1]
                    reverse_primer = reverse_complimentary_sequence(reverse_primer)
                elif any(True for x in long_read_schemes if self.scheme == x):
                    forward_primer = primer_seqs[0]
                    reverse_primer = primer_seqs[1]

                index_amplicon_list = primer_index_for_amplicon_seq(
                    reference_sequence=reference_sequence_seq,
                    forward_primer=forward_primer,
                    reverse_primer=reverse_primer
                    )
                '''amplicon dict is 
                amplicon_dict[reference_sequence_name] = {defline: sequence}'''

                if not index_amplicon_list:
                    defline = "~~~".join([primer_name, reference_sequence_defline, self.scheme, ])
                    if reference_fasta_basename in amplicon_dict:
                        amplicon_dict[reference_fasta_basename][defline] = None
                    else:
                        amplicon_dict[reference_fasta_basename] = {defline: None}
                else:
                    amplicion_seq = str(reference_sequence_seq[index_amplicon_list[0]:index_amplicon_list[1]])
                    amplicon_length = str(len(amplicion_seq)) + "bp"
                    defline = "~~~".join([primer_name, reference_sequence_defline, self.scheme, amplicon_length])
                    if reference_fasta_basename in amplicon_dict:
                        amplicon_dict[reference_fasta_basename][defline] = [amplicion_seq, index_amplicon_list]
                    else:
                        amplicon_dict[reference_fasta_basename] = {defline: [amplicion_seq, index_amplicon_list]}
        return amplicon_dict
    
    def write_amplicion_fasta_files(self, amplicon_storage_dir_pathway: str, amplicon_dict: dict) -> None:

        for ref_name, seq_dict in amplicon_dict.items():

            output_fasta = os.path.join(amplicon_storage_dir_pathway, ref_name + "_amplicons.fasta")
            output_log = os.path.join(amplicon_storage_dir_pathway, ref_name + "_amplicons.log")
            
            with open(output_fasta, 'w') as wf_amps, open(output_log, 'w') as wf_logs:
            
                for defline, amplicon_items_list in seq_dict.items():
                    primer_name = defline.split('~~~')[0]

                    if amplicon_items_list is None:
                        wf_logs.write(f"{primer_name}\tNo Amplification\n")

                    elif amplicon_items_list is not None:
                        amplicon_sequence, position_list = amplicon_items_list
                        start_pos, end_pos = position_list

                        line = '\t'.join([primer_name, ref_name, str(start_pos), str(end_pos), str(len(amplicon_sequence))])
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
                                    VarSkip = ["vsl1a", "vss1a", "vss2a", "vss2b"])', 
                                 required=True, 
                                 type=str,
                                 choices=["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vsl1a", "vss1a", "vss2a", "vss2b"],
                                 metavar='SCHEME',
                                 dest='scheme',
                                 default='vsl1a') #change this later
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

def GenerateAmplicon_MainWorkflow(fasta, scheme, sewage_home_dir, output, pathway):
    '''Main workflow for creating amplicons'''
    '''removed the opts as the main_argsparser is used'''

    sewage_amplicon = GenerateAmplicons(fasta, 
                                        scheme, 
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