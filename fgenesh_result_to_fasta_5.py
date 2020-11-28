'''
@Author: your name
@Date: 2020-08-01 23:20:38
@LastEditTime: 2020-08-01 23:58:24
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/fgenesh_result_to_fasta_5.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re

rgx_protein_start=re.compile('^>FGENESH:\s')
rgx_protein_end_1=re.compile('^>FGENESH:\[mRNA\]\s*')
rgx_protein_end_2=re.compile('^//')
rgx_contig_start=re.compile("Seq\s+name:(.+)Magnaporthe.+isolate(.+?)\s")
rgx_gene_end=re.compile("\>FGENESH\:\[exon\]")
def fgenesh_result_to_fasta(in_file_path,out_gene_file,out_protein_file):
    with open(out_protein_file,'w+') as out_protein_fl:
        with open(out_gene_file,'w+') as out_gene_fl:
            for file_path in in_file_path.rglob('result.txt'):
                if file_path.is_file():
                    with file_path.open() as in_fl:
                        # 撮箕一有两个部分，一个部分在撮箕二前面，另一部分在撮箕二后面。
                        for line in in_fl:
                            # 判断撮箕起点，撮箕一
                            contig_name_search=rgx_contig_start.search(line)
                            if contig_name_search is not None:
                                count_for_contig_id=1
                                contig_name=contig_name_search.group(1)
                                strain_name=contig_name_search.group(2)
                                for line in in_fl:
                                    search_protein_result=rgx_protein_start.search(line)
                                    search_gene_result=rgx_protein_end_1.search(line)
                                    if search_gene_result is not None:
                                        out_gene_fl.write(">{}_{}_gene_{}\n".format(strain_name.strip(),contig_name.strip(),count_for_contig_id))
                                        line=next(in_fl)
                                        # 判断撮箕一终点
                                        while(rgx_gene_end.search(line) is None):
                                            out_gene_fl.write(line)
                                            line=next(in_fl)
                                    if search_protein_result is not None:
                                        out_protein_fl.write(">{}_{}_gene_{}\n".format(strain_name.strip(),contig_name.strip(),count_for_contig_id))
                                        line=next(in_fl)
                                        while((rgx_protein_end_1.search(line) is None) and (rgx_protein_end_2.search(line) is None)):
                                            out_protein_fl.write(line)
                                            line=next(in_fl)
                                        count_for_contig_id=count_for_contig_id+1
                                        # 撮箕一起点和撮箕二终点重合，把撮箕一再拿来用，
                                        # 当前这一行是撮箕二终点行，正好是撮箕一的起点，一旦不加撮箕一，for循环的next把这行跳过，则会丢失撮箕一的内容
                                        search_gene_result=rgx_protein_end_1.search(line)
                                        if search_gene_result is not None:
                                            out_gene_fl.write(">{}_{}_gene_{}\n".format(strain_name.strip(),contig_name.strip(),count_for_contig_id))
                                            line=next(in_fl)
                                            while(rgx_gene_end.search(line) is None):
                                                out_gene_fl.write(line)
                                                line=next(in_fl)
                                        else:
                                            break

if __name__ == "__main__":
    fgenesh_result_path=Path('/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000030/')
    out_protein_file='../Pan_genome_data/TW-6-2-2-B-1_protein_5.fasta'
    out_gene_file='../Pan_genome_data/TW-6-2-2-B-1_gene_5.fasta'
    fgenesh_result_to_fasta(fgenesh_result_path,out_gene_file,out_protein_file)
    