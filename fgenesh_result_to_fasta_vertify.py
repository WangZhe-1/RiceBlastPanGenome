'''
@Author: your name
@Date: 2020-08-02 10:05:58
@LastEditTime: 2020-08-04 20:26:05
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/fgenesh_result_to_fasta_vertify.py
'''
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
import gffutils

rgx_protein_start=re.compile('^>FGENESH:\s')
rgx_protein_end_1=re.compile('^>FGENESH:\[mRNA\]\s*')
# 也是contig end，contig end 1
rgx_protein_end_2=re.compile('^//')
# contig end 2
rgx_contig_end=re.compile("no reliable predictions")
rgx_contig_start=re.compile("Seq\s+name:(?:(?P<contig>.+)Magnaporthe oryzae\s(?:isolate|strain)*(?P<strain>.+?)\s)|(?P<contig_ina168>Ina168.+?\s)")
rgx_gene_end=re.compile("\>FGENESH\:\[exon\]")
def fgenesh_result_to_fasta_verify(in_file_path,out_gene_file,out_protein_file,out_verify):
    with open(out_protein_file,'w+') as out_protein_fl:
        with open(out_gene_file,'w+') as out_gene_fl:
            with open(out_verify,'w') as out_verify_fl:
                for file_path in in_file_path.rglob('result.txt'):
                    if file_path.is_file():
                        with file_path.open() as in_fl:
                            # 撮箕一有两个部分，一个部分在撮箕二前面，另一部分在撮箕二后面。
                            for line in in_fl:
                                # 判断撮箕起点，撮箕一
                                contig_name_search=rgx_contig_start.search(line)
                                if contig_name_search is not None:
                                    count_for_contig_id=1
                                    contig_name=contig_name_search.group("contig")
                                    strain_name=contig_name_search.group("strain")
                                    if strain_name is None:
                                        contig_name=contig_name_search.group("contig_ina168")
                                        strain_name="ina168"
                                    elif strain_name=="genome":
                                        strain_name="FR13"
                                    contig_name=contig_name.strip()
                                    strain_name=strain_name.strip()
                                    for line in in_fl:
                                        search_protein_result=rgx_protein_start.search(line)
                                        search_gene_result=rgx_protein_end_1.search(line)
                                        if search_gene_result is not None:
                                            gene_id=">{}_{}_gene_{}\n".format(strain_name.strip(),contig_name.strip(),count_for_contig_id)
                                            out_gene_fl.write(gene_id)
                                            line_split=line.split()
                                            out_verify_fl.write("{}\t{}\t{}\n".format("{}_gene_{}".format(contig_name.strip(),count_for_contig_id),line_split[5],line_split[7]))
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
                                                gene_id=">{}_{}_gene_{}\n".format(strain_name.strip(),contig_name.strip(),count_for_contig_id)
                                                out_gene_fl.write(gene_id)
                                                line_split=line.split()
                                                out_verify_fl.write("{}\t{}\t{}\n".format("{}_gene_{}".format(contig_name.strip(),count_for_contig_id),line_split[5],line_split[7]))
                                                line=next(in_fl)
                                                while(rgx_gene_end.search(line) is None):
                                                    out_gene_fl.write(line)
                                                    line=next(in_fl)
                                            else:
                                                break
                                        if rgx_contig_end.search(line) is not None:
                                            break
def gff_mRNA(gff_file_path,out_file):
    db_file=Path("../Pan_genome_data/"+gff_file_path.stem+".db")
    if db_file.is_file() is False:
        gffutils.create_db(
        str(gff_file_path),
        dbfn=str(db_file), 
        force=True, 
        keep_order=True,
        id_spec=None
        )
    db=gffutils.FeatureDB(db_file, keep_order=True)
    with open(out_file,'w') as out_fl:
        for mRNA in db.features_of_type("mRNA"):
            out_fl.write("{}\t{}\t{}\n".format(
                list(db.parents(mRNA, featuretype='gene'))[0].id,
                mRNA.start,
                mRNA.end
            )
            )
if __name__ == "__main__":
    fgenesh_result_path=Path('/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000001/000082')
    out_protein_file='../Pan_genome_data/ina168_protein_verify.fasta'
    out_gene_file='../Pan_genome_data/ina168_gene_verify.fasta'
    out_verify_file_fasta="../Pan_genome_data/ina168_fasta_verify.txt"
    out_verify_file_GFF="../Pan_genome_data/ina168_GFF_verify.txt"
    fgenesh_result_to_fasta_verify(fgenesh_result_path,out_gene_file,out_protein_file,out_verify_file_fasta)
    gff_path=Path('/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000002/000082/result.txt')
    # gff_mRNA(gff_path,out_verify_file_GFF)
    
