import os
import sys
import glob

def generate_amplicons(
        referecne_fasta_raw: str, 
        scheme: str,
        sewage_dir: str,
        output: str,
        output_pathway: str):
    '''
    iterate through a dictionary with seq_name as key and sequence as values
    then find forward, reverse primers in sequence
    return dictionary with amplicion_name as ksy and amplicon sequence as value
    '''
    def pathway_or_fastfile(file_path):
        file_path = os.path.realpath(file_path)
        if os.path.isdir(file_path):
            file_pathway_list = glob.glob(file_path + "/*")
            valid_extensions = ('.fasta', '.fa', '.fsa')
            fasta_pathway_list = [file for file in file_pathway_list if os.path.splitext(file)[1].lower() in valid_extensions]
        elif os.path.isfile(file_path):
            fasta_pathway_list = [file_path]
        return fasta_pathway_list
    
    def read_referecne_fasta(fasta_pathway_list):
        '''
        open single or mulit-reference fasta file, or a pathway to fasta files that 
         end in .fasta .fa or .fsa and return single string sequence for each reference
        '''
        sequences = {}
        current_defline = ""  # Initialize with an empty string
        current_sequence = ""
        for pathway in fasta_pathway_list:
            with open(pathway, 'r') as fasta_file:
                for line in fasta_file:
                    line = line.strip()

                    if not line:
                        continue  # Skip empty lines

                    if line.startswith(">"):  # Def line
                        if current_defline:
                            sequences[current_defline] = current_sequence
                        current_defline = line[1:]
                        current_sequence = ""
                    else:
                        current_sequence += line

                # Add the last sequence after reaching the end of the file
                if current_defline and current_sequence:
                    sequences[current_defline] = current_sequence

        return sequences

    def read_primer_file(scheme: str, sewage_dir: str) -> dict:
        '''open custom primer file and return dict'''
        scheme_pathway = os.path.join(os.path.join(sewage_dir, "schemes"), scheme + ".tsv")
        primer_dict = {}
        with open(scheme_pathway, 'r') as fh:
            for line in fh:
                if line.startswith("#"): continue
                line = line.strip().split('\t')
                primer_dict[line[0]] = line[1:]
        return primer_dict

    def reverse_complimentary_sequence(sequence: str) -> str:
        '''
        return reverse complimentary sequence for the reverse primer
        only used with ARTIC and VarSkip short reads primers
        '''
        compliment_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "N": "N"
        }
        complementary_sequence_list = [compliment_dict.get(nucl, nucl) for nucl in sequence]
        return "".join(complementary_sequence_list[::-1])

    def find_primers_index(ref_seq: str, forward_primer: str, reverse_primer_revcomp: str) -> list:
        '''
        find exact match to forward and reverse primer sequences
        if at least one primer sequence is not found, return False
        else determing the starting index for the forward primer and
            the ending index of the reverse primer
        return a list of each index as an element
        '''
        forward_index = ref_seq.find(forward_primer)
        reverse_index_start = ref_seq.find(reverse_primer_revcomp)
        if any(x == -1 for x in [forward_index, reverse_index_start]):
            return False
        reverse_index_end = reverse_index_start + len(reverse_primer_revcomp)
        return [forward_index, reverse_index_end]

    def write_amplicions(amplicon_dict: dict, output_pathway: str, ref_fasta_name: str) -> None:
        output_fasta = os.path.join(output_pathway, ref_fasta_name + ".fasta")
        with open(output_fasta, 'w') as wf:
            for defline, sequence in amplicon_dict.items():
                wf.write(f">{defline}\n")
                while sequence:
                    chunk = sequence[:100]
                    sequence = sequence[100:]
                    wf.write(f"{chunk}\n")
        return

    def makeDir(output_pathway: str):
        output_pathway = os.path.realpath(output_pathway)
        if not os.path.exists(output_pathway):
            print(f"\nPathway {output_pathway} does not exist\n", file=sys.stderr)
            sys.exit(1)
        pathway = os.path.join(output_pathway, "SEWAGE_amplicons")
        if not os.path.isdir(pathway):
            os.mkdir(pathway)
        return pathway
    
    amplicon_pathway = makeDir(output_pathway=output_pathway)

    fasta_pathway_list = pathway_or_fastfile(referecne_fasta_raw)
    ref_seq_dict = read_referecne_fasta(fasta_pathway_list)
    primer_dict = read_primer_file(scheme=scheme, sewage_dir=sewage_dir)

    for ref_name, ref_sequence in ref_seq_dict.items():
        amplicon_dict = {}
        for primer_name, primer_seqs in primer_dict.items():

            if any(True for x in ["V1", "V2", "V3", "V4", "V4.1", "V5.3.2", "vss1a", "vss2a", "vss2b"] if x in scheme):
                forward_primer = primer_seqs[0]
                reverse_primer = reverse_complimentary_sequence(sequence=primer_seqs[1])
            elif "vsl1a" in scheme:
                forward_primer = primer_seqs[0]
                reverse_primer = primer_seqs[1]

            index_amplicon_list = find_primers_index(
                ref_seq=ref_sequence,
                forward_primer=forward_primer,
                reverse_primer_revcomp=reverse_primer
            )

            if not index_amplicon_list:
                defline = "__".join([primer_name, scheme, ref_name])
                amplicon_dict[defline] = None
            else:
                amplicion_seq = str(ref_sequence[index_amplicon_list[0]:index_amplicon_list[1]])
                amplicon_length = str(len(amplicion_seq)) + "bp"
                defline = "__".join([primer_name, scheme, ref_name, amplicon_length])
                amplicon_dict[defline] = amplicion_seq
    
        write_amplicions(amplicon_dict=amplicon_dict, 
                         output_pathway=amplicon_pathway,
                         ref_fasta_name=ref_name)

    return
    
# if __name__ == "__main__":
#     generate_amplicons(sys.argv[1], sys.argv[2], sys.argv[3])
'''
still need to connect this script to the amplicon subcommand using argsparse
still need to add in all the available primers lists
'''