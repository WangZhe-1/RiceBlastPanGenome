'''
Author: your name
Date: 2020-10-18 11:06:40
LastEditTime: 2020-10-18 11:29:27
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/mafft.py
'''
import subprocess
from joblib import Parallel,delayed
def mafft_popen(Single_Copy_Orthologue_file_path,out_dir_path,out_stderr_dir_path):
    ma_popen=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/miniconda3/bin/mafft",
            "--genafpair",
            "--maxiterate",
            "16",
            "--reorder",
            Single_Copy_Orthologue_file_path.resolve()
        ],
        stdout=(out_dir_path/Single_Copy_Orthologue_file_path.stem).open("w"),
        stderr=(out_stderr_dir_path/(Single_Copy_Orthologue_file_path.stem+"_stderr.txt")).open("w")
    )
    ma_popen.wait()
def mafft_run(Single_Copy_Orthologue_dir_path,mafft_out_dir_path,mafft_stderr_dir_path):
    '''
    input 1: Single_Copy_Orthologue_dir_path
    output 1: mafft_out_dir_path
    outpur 2: mafft_stderr_dir_path
    '''
    Parallel(n_jobs=1)(
        delayed(mafft_popen)(
            Single_Copy_Orthologue_file_path,
            mafft_out_dir_path,
            mafft_stderr_dir_path
            ) 
            for Single_Copy_Orthologue_file_path in Single_Copy_Orthologue_dir_path.iterdir()
        )