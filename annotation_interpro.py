'''
Author: your name
Date: 2020-09-09 13:51:37
LastEditTime: 2020-09-09 14:19:00
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/annotation_interpro.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
pan_protein=Path('../Pan_genome_data_2/pan/pan_protein_no_shorter_20.fasta')
out_file='../Pan_genome_data_2/annotation/'


interpro_cmd=subprocess.Popen(
    [
        "/mnt/d/zhes_learning_space/software_in_ubuntu/interproscan-5.46-81.0/interproscan.sh",
        "-i",
        pan_protein,
        "-f",
        "GFF3",
        "-d",
        out_file,
        "-appl",
        "Pfam",
        "-goterms",
        "-pa",
        "-iprlookup",
        "-cpu",
        "12"
    ],
    stdout=open('../Pan_genome_data_2/annotation/interpro_out','w+'),
    stderr=open('../Pan_genome_data_2/annotation/interpro_err',"w+"),
    universal_newlines=True,
    encoding='UTF-8',
)
interpro_cmd.wait()
