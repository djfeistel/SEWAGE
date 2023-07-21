import logging

def logFile(parent_dir):
    dir_name = os.path.basename(parent_dir)
    log_filename = os.path.join(parent_dir, f"{dir_name}_logfile.log")
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    return

def logFile_initials(args, fasta_in:str, fasta_type:str, parent_dir:str, fasta_pathway_list:list):
    args_dict = vars(args)
    args_list = list(args_dict.items())
    args_list = [(str(i), str(j)) for i, j in args_list]
    logging.info("SEWAGE initiated")
    logging.info("Parameters used in SEWAGE:")
    for item in args_list:
        logging.info(f"{item[0]}\t{item[1]}")
    logging.info("Additional information")
    logging.info(f"\nFastA files found in {fasta_in}")
    logging.info(f"FastA input type: {fasta_type}")
    logging.info(f"Total FASTA files detected: {len(fasta_pathway_list)}")
    logging.info(f"Output storage: {parent_dir}")
    return