'''
@Author: your name
@Date: 2020-08-02 15:46:50
@LastEditTime: 2020-08-02 21:16:10
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/copy_for_protein_2.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
from extract_strain_id import extract_strain_id

# input 1
# R code contain busco result


def Copy_contig(input_path,MGG_in,ina168_in,output_dir_path,std_out_err_path):
    '''
    input1: copy file
    input 2: 70-15 file path
    input 3: ina168 file path
    output1: copy to dir path
    output 2: copy std out err
    '''
    strain_95_list=extract_strain_id("../Pan_genome_data/read_busco.txt")
    for query_id in strain_95_list:
        file_path=input_path/(query_id.strip('\n')+'.fasta')
        gene_cp=subprocess.Popen(
        ["cp",str(file_path),str(output_dir_path)],
        stdout=(std_out_err_path/'cp_{}_stdout.txt'.format(file_path.stem)).open('w'),
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        encoding='UTF-8',
        )
        gene_cp.wait()

    gene_cp=subprocess.Popen(
    ["cp",MGG_in,str(output_dir_path)],
    stdout=(std_out_err_path/'cp_{}_stdout.txt'.format("70-15")).open('w'),
    stderr=subprocess.STDOUT,
    universal_newlines=True,
    encoding='UTF-8',
    )

    gene_cp=subprocess.Popen(
    ["cp",ina168_in,str(output_dir_path)],
    stdout=(std_out_err_path/'cp_{}_stdout.txt'.format("ina168")).open('w'),
    stderr=subprocess.STDOUT,
    universal_newlines=True,
    encoding='UTF-8',
    )
