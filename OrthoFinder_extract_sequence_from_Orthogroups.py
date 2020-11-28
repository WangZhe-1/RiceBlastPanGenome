#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

from Bio import SeqIO
from pathlib import Path

#input 1 orthogroup
# Orthogroup_Sequences_path=Path("../../wangzhe2/OrthoFinder_result/Results_Feb11/Orthogroup_Sequences/test/")
Orthogroup_Sequences_path=Path("../../wangzhe2/OrthoFinder_result/Results_Feb11/Orthogroup_Sequences/")
#output 1:pan genome sequence
out_pan=open ("../Pan_genome_data/pan_protein.fasta","w+")
#output 2:pan genome og_id protein_id
out_id=open("../Pan_genome_data/pan_id.txt","w+")
def sort_head(std_in,flag,group_id):
    global out_pan
    # 1表示有70-15的基因，以grep结果排序
    if flag==1:
        sort_cmd=subprocess.Popen(
            [
                "seqkit",
                "sort",
                "-l",
                "-r",
                "--quiet"
            ],
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            encoding='utf-8'
        )
    else:
        # 0表示没有70-15的基因，用原文件排序
        sort_cmd=subprocess.Popen(
            [
                "seqkit",
                "sort",
                "-l",
                "-r",
                "--quiet",
                std_in
            ],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            encoding='utf-8'
        )
    head_cmd=subprocess.Popen(
        [
            "seqkit",
            "head",
            "-n",
            "1"
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        encoding='utf-8'
    )
    (sort_out,sort_err)=sort_cmd.communicate(input=std_in)    
    if sort_cmd.returncode!=0:
        for line_err in sort_err.split('\n'):
            print (line_err)
        return
    (head_out,head_err) = head_cmd.communicate(input=sort_out)
    if head_cmd.returncode!=0:
        for line_err in head_err.split('\n'):
            print (line_err)
        return
    out_pan.write(head_out)
    protein_id=head_out.split('\n')[0].strip('>').strip()
    out_id.write("{}\t{}\n".format(group_id,protein_id))
    
for Orthogroup_Sequences_file in Orthogroup_Sequences_path.iterdir():
    if Orthogroup_Sequences_file.is_dir(): continue
    seqkit_cmd=subprocess.Popen(
        [
            "seqkit",
            "grep",
            "-r",
            "-p",
            "MGG_",
            str(Orthogroup_Sequences_file)
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        encoding='utf-8'
    )
    (std_out,std_err)=seqkit_cmd.communicate()
    if std_err.__len__()==0:
        if std_out.__len__()==0:
            sort_head(str(Orthogroup_Sequences_file),0,Orthogroup_Sequences_file.stem)
        else:
            sort_head(std_out,1,Orthogroup_Sequences_file.stem)
    else:
        for line_err in std_err.split('\n'):
            print (line_err)

out_pan.close()
out_id.close()

        