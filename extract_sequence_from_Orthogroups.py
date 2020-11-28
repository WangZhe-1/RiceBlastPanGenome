#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

from Bio import SeqIO
from pathlib import Path

out_pan=None
out_id=None
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
                "--quiet"
            ],
            stdin=subprocess.PIPE,
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

def extract_sequence_from_Orthogroups(Orthogroup_Sequences_path_name,out_pan_name,out_id_name):
    '''
    input 1 orthogroup_sequence_path_name
    output 1:pan genome sequence
    output 2:pan genome og_id protein_id
    '''
    global out_pan,out_id
    out_pan=open(out_pan_name,'w+')
    out_id=open(out_id_name,'w+')
    Orthogroup_Sequences_path=Path(Orthogroup_Sequences_path_name)
    for Orthogroup_Sequences_file in Orthogroup_Sequences_path.iterdir():
        if Orthogroup_Sequences_file.is_dir(): continue
        seqkit_trim_x_cmd=subprocess.Popen(
            [
                "seqkit",
                "grep",
                "-s",
                "-i",
                "-v",
                "-p",
                "x",
                str(Orthogroup_Sequences_file)
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            encoding='utf-8'
        )
        (trim_x_out,trim_x_err)=seqkit_trim_x_cmd.communicate()
        if trim_x_err.__len__()==0:
            if trim_x_out.__len__()==0:
                continue
        else:
            for line_err in trim_x_err.split('\n'):
                print (line_err)
        seqkit_cmd=subprocess.Popen(
            [
                "seqkit",
                "grep",
                "-r",
                "-p",
                "MGG_"
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            encoding='utf-8'
        )
        (std_out,std_err)=seqkit_cmd.communicate(input=trim_x_out)
        if std_err.__len__()==0:
            if std_out.__len__()==0:
                sort_head(trim_x_out,0,Orthogroup_Sequences_file.stem)
            else:
                sort_head(std_out,1,Orthogroup_Sequences_file.stem)
        else:
            for line_err in std_err.split('\n'):
                print (line_err)
    out_pan.close()
    out_id.close()



        