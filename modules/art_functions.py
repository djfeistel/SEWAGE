import os
import subprocess

def check_art_illumina_command():
    try:
        subprocess.run(["art_illumina"], 
                       check=True, 
                       stdout=subprocess.DEVNULL, 
                       stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        return [False, FileNotFoundError]
    except subprocess.CalledProcessError:
        return [True, subprocess.CalledProcessError]


def artIlluminaSubprocess(fastaFile:str,
                          log_dir:str,
                          proportion:float,
                          pfold:int=1000,
                          ss:str='HS25',
                          l:int=150,
                          m:int=200,
                          s:int=1,
                          ir:int=0,
                          ir2:int=0,
                          dr:int=0,
                          dr2:int=0,
                          rndSeed:int=13,
                          maxIndel:int=0,
                          nf:int=0,
                          qL:int=28,
                          qU:int=40,
                          qs:int=0,
                          qs2:int=0):
    '''simulate reads'''
    file_name = os.path.basename(fastaFile) #removed pathway for as fasta file only
    id = file_name[:file_name.rfind('.')] # remove extention to add to reads defline
    fcov = float(proportion)*pfold

    command = [
        "art_illumina",
        "-ss", ss,
        "-i", fastaFile,
        "--paired",
        "-l", str(l),
        "-f", str(fcov),
        "-m", str(m),
        "-s", str(s),
        "-o", f"{id}_",
        "-d", f"FASTAid_{id}_proportion_{round(proportion, 3)}_",
        "-ir", str(ir),
        "-ir2", str(ir2),
        "-dr", str(dr),
        "-dr2", str(dr2),
        "-rs", str(rndSeed),
        "-k", str(maxIndel),
        "-na",
        "-nf", str(nf),
        "-qL", str(qL),
        "-qU", str(qU),
        "-qs", str(qs),
        "-qs2", str(qs2)
    ]

    stderr_file = f"{id}_log.txt"
    log_pathway = os.path.join(log_dir, stderr_file)
    with open(log_pathway, 'w') as stderr_output:
        process = subprocess.Popen(' '.join(command), shell=True, stdout=stderr_output, stderr=stderr_output)
        process.wait()
    return