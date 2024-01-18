import time
import os
import glob
import random
import string
import shutil

def countdown():
    print(f"\nSimulating Sewage reads starts in...\n")
    for i in range(3, 0, -1):
        print(f"{i} seconds\n")
        time.sleep(1)
    print(f"Blast off!\n")
    return

def writeProportionFile(fasta_pathway_list:str, parent_dir:str)->None:
    with open(os.path.join(parent_dir, 'proportions_list.txt'), 'w') as wf:
        for file, prop in fasta_pathway_list:
            wf.write(f"{file}\t{prop}\n")

def parentDir_pathway(dir_pathway:str, dir_name:str):
    return os.path.realpath(os.path.join(dir_pathway, dir_name))


def concatenateFastQfiles(rawreads_dir, parent_dir, output_name):
    '''Concatenate fastq files into a single file'''
    forward_reads = glob.glob(rawreads_dir + '/*_1.fq')
    reverse_reads = glob.glob(rawreads_dir + '/*_2.fq')

    def writeProportionFastQ(parent_dir, output_name, list_of_fastq_files, extention):
        output_file = output_name + extention
        pathway = os.path.join(parent_dir, output_file)
        with open(pathway, 'wb') as wf:
            for fastq_file in list_of_fastq_files:
                with open(fastq_file, 'rb') as rf:
                    shutil.copyfileobj(rf, wf)
        return
    
    writeProportionFastQ(parent_dir, output_name, forward_reads, "_1.fastq")
    writeProportionFastQ(parent_dir, output_name, reverse_reads, "_2.fastq")

    return