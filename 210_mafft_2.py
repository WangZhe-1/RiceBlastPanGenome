'''
Author: your name
Date: 2020-10-27 08:56:03
LastEditTime: 2020-10-28 21:01:19
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/210_mafft.py
'''
from dask_jobqueue import PBSCluster
from dask import delayed
import dask
cluster = PBSCluster(
    job_extra=["-l nodes=1:ppn=24","-l mem=5000MB"],
    header_skip=["select"],
    processes=24,
    walltime='25:00:00'
    )
cluster.scale(jobs=10)
from dask.distributed import Client
client = Client(cluster)

print(cluster.job_script())
import subprocess
from pathlib import Path


def mafft_popen(gene_file_path,out_dir_path,out_stderr_dir_path):
    out_file=out_dir_path/gene_file_path.stem
    if out_file.exists() is True: return
    ma_popen=subprocess.Popen(
        [
            "mafft",
            "--genafpair",
            "--maxiterate",
            "16",
            "--reorder",
            gene_file_path.resolve()
        ],
        stdout=out_file.open("w"),
        stderr=(out_stderr_dir_path/(gene_file_path.stem+"_stderr.txt")).open("w")
    )
    ma_popen.wait()

gene_file_dir_path=Path("/gpfshome/home/Baojd/wangzhe/Single_Copy_Orthologue_gene")
calls=[]
for gene_file_path in gene_file_dir_path.iterdir():
    calls.append(
        delayed(mafft_popen)(
            gene_file_path,
            Path("/gpfshome/home/Baojd/wangzhe/mafft_out/"),
            Path("/gpfshome/home/Baojd/wangzhe/mafft_stderr/")
        )
        )
dask.compute(*calls)