'''
@Author: your name
@Date: 2020-07-31 11:17:33
@LastEditTime: 2020-07-31 17:26:35
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/fgenesh_result_to_fasta.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
from re import search

rgx_protein_start=re.compile('^>FGENESH:\s')
rgx_protein_end_1=re.compile('^>FGENESH:\[mRNA\]\s*')
rgx_protein_end_2=re.compile('^//')
rgx_gene_end=re.compile("\>FGENESH\:\[exon\]")
def fgenesh_result_to_fasta(in_file_path,out_gene_file,out_protein_file):
    with open(out_protein_file,'w+') as out_protein_fl:
        with open(out_gene_file,'w+') as out_gene_fl:
            for file_path in in_file_path.rglob('result.txt'):
                if file_path.is_file():
                    with file_path.open() as in_fl:
                        for line in in_fl:
                            search_protein_result=rgx_protein_start.search(line)
                            if search_protein_result is not None:
                                out_protein_fl.write(line)
                                line=next(in_fl)
                                while((rgx_protein_end_1.search(line) is None) and (rgx_protein_end_2.search(line) is None)):
                                    out_protein_fl.write(line)
                                    line=next(in_fl)
                        in_fl.seek(0)
                        for line in in_fl:
                            search_gene_result=rgx_protein_end_1.search(line)
                            if search_gene_result is not None:
                                out_gene_fl.write(line)
                                line=next(in_fl)
                                while(rgx_gene_end.search(line) is None):
                                    out_gene_fl.write(line)
                                    line=next(in_fl)

if __name__ == "__main__":
    fgenesh_result_path=Path('/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000030/')
    out_protein_file='../Pan_genome_data/TW-6-2-2-B-1_protein.fasta'
    out_gene_file='../Pan_genome_data/TW-6-2-2-B-1_gene.fasta'
    fgenesh_result_to_fasta(fgenesh_result_path,out_gene_file,out_protein_file)