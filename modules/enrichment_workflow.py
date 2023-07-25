import sys
import os
import glob
import logging
import random
import string
import shutil
import subprocess
import numpy as np
import time
import logging
from tqdm import tqdm

modules_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "modules")
sys.path.append(modules_dir)
# internal functions
from log_functions import logFile, logFile_initials
from fasta_funcitons import globFastAfiles, listOfFastaPathways
from proportion_functions import equalProportions, randomProportions
from art_functions import check_art_illumina_command, artIlluminaSubprocess
from write_functions import (
    countdown,
    writeProportionFile,
    parentDir_pathway,
    concatenateFastQfiles,
)


def EnrichmentWorkflow(
    args,
    fasta_in: str,
    dir_name: str,
    out_pathway: str,
    proportion: int,
    rndSeed: int,
    pfold: int,
    ss: str,
    l: int,
    m: int,
    s: int,
    ir: int,
    ir2: int,
    dr: int,
    dr2: int,
    rdmSd_art: int,
    maxIndel: int,
    nf: int,
    qL: int,
    qU: int,
    qs: int,
    qs2: int,
):
    
    
    #args = args

    """check if art_illumina is in pathway before starting any work"""
    if not check_art_illumina_command()[0]:
        print(
            f"\n'art_illumina' was not found in your pathway {check_art_illumina_command()[1]}\n",
            file=sys.stderr,
        )
        sys.exit(1)

    out_pathway = os.path.realpath(out_pathway)

    """check if fasta_in is a pathway or file exists or None"""
    fasta_pathway_list = None
    fasta_type = None
    if fasta_in is None:
        print(
            f"\nPlease specify a pathway or file name for the -i flag\n",
            file=sys.stderr,
        )
        sys.exit(0)
    elif os.path.isdir(fasta_in):
        fasta_pathway_list = globFastAfiles(fasta_in)
        fasta_type = "PATHWAY"
    elif os.path.isfile(fasta_in):
        fasta_pathway_list = listOfFastaPathways(fasta_in)
        fasta_type = "FILE"
    else:
        print(f"\n{fasta_in} is neither a file nor a pathway.\n", file=sys.stderr)
        sys.exit(1)

    """check if out_pathway and dir_name exits"""
    if not os.path.isdir(out_pathway):
        print(f"\n{out_pathway} is not a valid pathway\n", file=sys.stderr)
        sys.exit(1)

    """make storage dir"""
    parent_dir = parentDir_pathway(dir_pathway=out_pathway, dir_name=dir_name)
    try:
        os.mkdir(parent_dir)
        rawreads_dir = os.path.join(parent_dir, "00.RAWREADS")
        logs_dir = os.path.join(parent_dir, "99.LOGS")  # logs for art commands
        os.mkdir(rawreads_dir)
        os.mkdir(logs_dir)
    except FileExistsError:
        print(
            f"\nDirectory '{os.path.basename(parent_dir)}' already exists in pathway. Choose a different name.\n",
            file=sys.stderr,
        )
        sys.exit(0)
    except OSError as e:
        print(f"Failed to create directory: {e}")
        sys.exit(0)

    """
    if you have made it this far, there is no turning back...
    your journy starts now...
    """
    """logging starts here"""
    logFile(parent_dir=parent_dir)  # initiat logging file
    logFile_initials(
        args=args,
        fasta_in=fasta_in,
        fasta_type=fasta_type,
        parent_dir=parent_dir,
        fasta_pathway_list=fasta_pathway_list,
    )

    """which proportion do you want?"""

    if proportion == "e":
        fasta_proportions_list = equalProportions(fasta_pathway_list=fasta_pathway_list)
    else:
        fasta_proportions_list = randomProportions(
            fasta_pathway_list=fasta_pathway_list, random_seed=rndSeed
        )

    """change to rawreads directory to add reads there"""
    try:
        os.chdir(rawreads_dir)
        # print(f"Data will be saved in: {os.getcwd()}\n", file=sys.stderr)
        logging.info(f"Data saved in: {rawreads_dir}")
        logging.info(f"Simulated information saved in: {logs_dir}")
    except Exception as e:
        print(f"Error occurred while changing to the directory: {e}")

    """just for fun :)"""
    countdown()

    for item in tqdm(
        fasta_proportions_list, desc="Processing", unit="tuple", dynamic_ncols=False
    ):
        # for fasta_file, fasta_proportion in fasta_proportions_list:
        fasta_file, fasta_proportion = item
        logging.info(f"Simulation reads for {fasta_file}")
        logging.info(f"Proporiton is set at {fasta_proportion}")

        artIlluminaSubprocess(
            fastaFile=fasta_file,
            log_dir=logs_dir,
            proportion=fasta_proportion,
            pfold=pfold,
            ss=ss,
            l=l,
            m=m,
            s=s,
            ir=ir,
            ir2=ir2,
            dr=dr,
            dr2=dr2,
            rndSeed=rdmSd_art,
            maxIndel=maxIndel,
            nf=nf,
            qL=qL,
            qU=qU,
            qs=qs,
            qs2=qs2,
        )

    """concatinate fastq files into one file with basename of parent_dir as name"""
    concatenateFastQfiles(
        parent_dir=parent_dir,
        rawreads_dir=rawreads_dir,
        output_name=os.path.basename(parent_dir),
    )
    """write a file with the pathways and proporitons"""
    writeProportionFile(
        fasta_pathway_list=fasta_proportions_list, parent_dir=parent_dir
    )

    print(f"\nFinished creating 'Sewage' samples :)\n", file=sys.stderr)
    logging.info(f"Finished creating 'Sewage' samples :)")
