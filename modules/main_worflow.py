import os
import sys
from modules import fasta_utils
from modules import amplicon_utils
from modules import proportion_utils
from modules import read_generator_utils
from modules import argsparse_opts
from modules import storage_utils

def main():
    
    args = argsparse_opts.sewage_opts()
    args.infasta = os.path.realpath(args.infasta)
    
    infasta = args.infasta
    scheme = args.scheme
    file_prefix_name = args.file_prefix_name
    storage_dir = args.storage_dir
    time_stamp = args.time_stamp
    proportion_model = args.proportion_model
    dVOC_genome = args.dVOC_genome
    dVOC_proporiton = args.dVOC_proporiton
    proportion_seed = args.proportion_seed
    read_length = args.read_length
    frag_length = args.frag_length
    coverage_depth = args.coverage_depth
    read_seed = args.read_seed

    ### fasta utilities ###
    fastaUtilsInstance = fasta_utils.FastaUtils(
        fasta=infasta,
        file_prefix_name=file_prefix_name
    )

    if fastaUtilsInstance.check_input_for_fasta_or_pathways():
        reference_fasta_dict = fastaUtilsInstance.read_multifasta_into_dict()
    else:
        fastaUtilsInstance.create_list_of_pathways()
        try:
            if fastaUtilsInstance.check_pathway_list_for_fasta():
                reference_fasta_dict = fastaUtilsInstance.read_fasta_pathways_into_dict()
            else:
                print("Error: The pathway list does not contain valid FASTA files.")
                sys.exit(1)  # Terminate the program with a non-zero exit code
        except Exception as e:
            print(f"An error occurred: {e}")
            sys.exit(1)

    reference_fasta_dict = fastaUtilsInstance.get_fasta_dict()

    if not reference_fasta_dict:
        raise ValueError(f"{infasta} file was not imported correctly")
    
    ### amplicon utilities ###
        
    ampliconUtilsInstance = amplicon_utils.GenerateAmplicons(
        reference_fasta_dict=reference_fasta_dict, 
        scheme=scheme,
        file_prefix_name=file_prefix_name
    )

    ampliconUtilsInstance.create_amplicon_dataframe()
    df_amplicon = ampliconUtilsInstance.get_df_amplicon()

    if df_amplicon.empty:
        raise ValueError("No amplicons were generated. \
                         Check reference genomes for compatability with {scheme} primer scheme.")
    
    ### proporitons utilities ###

    proportionUtilsInstance = proportion_utils.GenomeProporitons(
        reference_fasta_dict=reference_fasta_dict,
        dVOC_genome=dVOC_genome,
        dVOC_proporiton=dVOC_proporiton,
        random_seed=proportion_seed
        )
    
    if proportion_model == "r":
        proportionUtilsInstance.random_proportions()
    elif proportion_model == "e":
        proportionUtilsInstance.equal_proportions()
    elif proportion_model == "d":
        if dVOC_genome is not None:
            proportionUtilsInstance.dvoc_proportions_genome_assignment()
        else:
            proportionUtilsInstance.dvoc_proportions_random_genome_assignment()

    proportions_dict = proportionUtilsInstance.get_proportions_dict()

    if not proportions_dict:
        raise ValueError("No primers found")
    
    ### generate reads utilities ###

    readGeneratorUtilsInstance = read_generator_utils.ReadGenerator(
        df_amplicon=df_amplicon, 
        proportion_dict=proportions_dict,
        fastq_name=file_prefix_name,
        propotion_file_name=file_prefix_name,
        read_length=read_length,
        frag_length=frag_length,
        coverage_depth=coverage_depth,
        seed=read_seed
        )
    
    readGeneratorUtilsInstance.create_reads_workflow()

    df_reads = readGeneratorUtilsInstance.get_df_reads()
    
    ### save results ###
    
    storageUtilsInstance = storage_utils.StorageUtils(
        storage_dir=storage_dir, 
        time_stamp=time_stamp
    )

    storageUtilsInstance.add_time_stamp()
    storageUtilsInstance.create_storage_dir()
    storage_pathway = storageUtilsInstance.get_storage_pathway()

    # # save fasta
    argsparse_opts.write_parameters_log(args=args, file_prefix_name=file_prefix_name, storage_pathway=storage_pathway)
    fastaUtilsInstance.write_reference_fasta(storage_pathway=storage_pathway)
    ampliconUtilsInstance.save_amplicon_fasta(storage_pathway=storage_pathway)
    readGeneratorUtilsInstance.write_fastq_files(storage_pathway=storage_pathway)
    readGeneratorUtilsInstance.save_read_df(storage_pathway=storage_pathway)
    





    