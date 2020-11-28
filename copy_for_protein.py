#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import re
from pathlib import Path
from extract_strain_id import extract_strain_id

# input 1
# R code contain busco result
strain_95_list=extract_strain_id("../Pan_genome_data/a_qc_contig_busco_read_busco.txt")
# 输入
base_path=Path('../../GFF/')

for query_id in strain_95_list:
    for file in base_path.glob(query_id.strip('\n')+'*_protein.fasta'):
        gene_cp=subprocess.Popen(
        ["cp",str(file),'../Pan_genome_data/b_orthofinder_input_protein/'],
        stdout=open('../Pan_genome_data/out/cp_{}_stdout.txt'.format(file.stem),'w+'),
        stderr=open('../Pan_genome_data/out/cp_{}_stderr.txt'.format(file.stem),'w+'),
        universal_newlines=True,
        encoding='UTF-8',
        )
        gene_cp.wait()

gene_cp=subprocess.Popen(
["cp","../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_proteins_T0.fasta",'../Pan_genome_data/b_orthofinder_input_protein/'],
stdout=open('../Pan_genome_data/out/cp_{}_stdout.txt'.format("70-15"),'w+'),
stderr=open('../Pan_genome_data/out/cp_{}_stderr.txt'.format("70-15"),'w+'),
universal_newlines=True,
encoding='UTF-8',
)

gene_cp=subprocess.Popen(
["cp","../New_add_ina168/out.aa",'../Pan_genome_data/b_orthofinder_input_protein/ina168_protein.fasta'],
stdout=open('../Pan_genome_data/out/cp_{}_stdout.txt'.format("ina168"),'w+'),
stderr=open('../Pan_genome_data/out/cp_{}_stderr.txt'.format("ina168"),'w+'),
universal_newlines=True,
encoding='UTF-8',
)
