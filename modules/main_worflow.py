import numpy as np
import pandas as pd
import random
import argparse
import os
import fasta_funcitons_class_workshop
import amplicons_functions_class_worshop
import proportion_functions_workshop
import read_generator_worshop
from argsparse_opts_workshop import sewage_opts

if __name__ == "__main__":

    args = sewage_opts()

    infasta = args.infasta
    scheme = args.scheme
    amplicon_fasta_name = args.amplicon_fasta_name
    amplicon_storage_dir = args.amplicon_storage_dir
    proportion_model = args.proportion_model
    dVOC_genome = args.dVOC_genome
    dVOC_proporiton = args.dVOC_proporiton
    proportion_seed = args.proportion_seed
    fastq_name = args.fastq_name
    read_length = args.read_length
    coverage_depth = args.coverage_depth
    read_seed = args.read_seed

    # load fasta data
    load_fasta = fasta_funcitons_class_workshop.LoadFasta(fasta=infasta)
    if load_fasta.check_input_for_fasta_or_pathways():
        reference_fasta_dict = load_fasta.read_multifasta_into_dict()
    else:
        pathway_list = load_fasta.create_list_of_pathways()
        reference_fasta_dict = load_fasta.read_fasta_pathways_into_dict(pathway_list)
    
    # amplicon generation 
        
    amplicon_generator = amplicons_functions_class_worshop.GenerateAmplicons(
        fasta_dict=reference_fasta_dict, 
        scheme=scheme,
        amplicon_fasta_name=amplicon_fasta_name,
        amplicon_storage_dir=amplicon_storage_dir
    )
    #use storage_dir as input to read_generator
    storage_dir = amplicon_generator.check_amplicon_storage_dir()
    primer_scheme_dict = amplicon_generator.scheme_primer_dictionary()
    df_amplicon = amplicon_generator.amplicon_dataframe(
        genome_sequence_dict=reference_fasta_dict,
        primer_scheme_dict=primer_scheme_dict
        )
    
    
    amplicon_generator.write_fasta(df=df_amplicon)
    #for whatever reason i need to filter No Amplified amplicons here
    df_amplicon = df_amplicon[df_amplicon['amplicon_sequence'] != "No Amplification"]

    #generate proportions
    genome_props = proportion_functions_workshop.GenomeProporitons(
        reference_fasta_dict=reference_fasta_dict,
        dVOC_genome=dVOC_genome,
        dVOC_proporiton=dVOC_proporiton,
        random_seed=proportion_seed
        )
    
    if proportion_model == "e":
        proportion_dict = genome_props.equal_proportions()
    elif proportion_model == "r":
        proportion_dict = genome_props.random_proportions()
    elif proportion_model == "d":
        if dVOC_genome is not None:
            proportion_dict = genome_props.dvoc_proportions_genome_assignment()
        else:
            proportion_dict = genome_props.dvoc_proportions_random_genome_assignment()

    #generate reads
    read_generator = read_generator_worshop.ReadGenerator(
        df_amplicon=df_amplicon, 
        proportion_dict=proportion_dict,
        storage_dir=storage_dir,
        fastq_name=fastq_name,
        read_length=read_length,
        coverage_depth=coverage_depth,
        seed=read_seed
        )
    read_generator.add_proportion_column()
    read_generator.create_reads_from_amplicons()
    read_generator.write_reads()