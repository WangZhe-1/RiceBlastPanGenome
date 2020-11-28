'''
Author: your name
Date: 2020-09-11 17:22:38
LastEditTime: 2020-09-29 09:36:06
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/extract_nucleic_2.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 把1891序列读进去，拿出核酸id（摘T之前的），核酸读到index结构，用id找对应的核酸
from joblib import Parallel,delayed
import re
import subprocess
from pathlib import Path
from Bio import SeqIO
import gffutils
from extract_strain_id import extract_strain_id
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet
def merge_to_one(GFF_path_name,list_file,mRNA_base,cat_err_file):
    '''
    input 1: GFF path
    output 1: list file
    output 2: gene base
    output 3: cat_err_file
    '''
    fasta_list=[
        'cat'
        ]
    GFF_path=Path(GFF_path_name)
    with open (list_file,'w+') as list_out:
        for fasta_file in GFF_path.glob("*.fasta"):
            list_out.write('{}\n'.format(fasta_file))
            fasta_list.append(str(fasta_file.resolve()))
        
    fasta_cat=subprocess.Popen(
        fasta_list,
        stdout=open(mRNA_base,'w+'),
        stderr=open(cat_err_file,'w+'),
        universal_newlines=True,
    )
    fasta_cat.wait()

def extract_mRNA(protein_file,mRNA_base_name,pan_mRNA_file_name):
    '''
    input 1: pan protein fasta file
    input 2: merged gene file
    output 1:pan_mRNA_file_name
    '''
    mRNA_base=SeqIO.index(
        mRNA_base_name,
        "fasta",
        )
    with open(pan_mRNA_file_name,'w+') as out_fl:
        for protein_sequence in SeqIO.parse(protein_file,"fasta"):
            mRNA_sequence=mRNA_base[protein_sequence.id]
            SeqIO.write(mRNA_sequence,out_fl,"fasta")
def make_gffutils_db(strain_gff,db_out_file):
    try:
        gffutils.create_db(
        str(strain_gff),
        dbfn=str(db_out_file/(strain_gff.stem+".db"))
        )
    except:
        pass
def contig_gff3_gffutils_db(strain_gff_dir_path,db_out_dir):
    Parallel(n_jobs=6)(delayed(make_gffutils_db)(i,db_out_dir) for i in strain_gff_dir_path.glob("*.gff"))
def extract_gene_gff(protein_file_path,strain_db_dir_path,contig_path,out_gene_file):
    global_names=globals()
    for db_file in strain_db_dir_path.iterdir():
        global_names[db_file.stem.strip()+"_db"]=gffutils.FeatureDB(db_file)
    MGG_db=SeqIO.index("../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta","fasta")
    with out_gene_file.open("w") as out_fl:
        for protein_sequence in SeqIO.parse(protein_file_path,"fasta"):
            protein_id_list=protein_sequence.id.split("_",1)
            strain_id=protein_id_list[0]
            gff_protein_id=protein_id_list[1]
            if strain_id=="MGG":
                SeqIO.write(MGG_db[protein_sequence.id[0:9]],out_fl,"fasta")
            else:
                if strain_id=="WD-3-1":
                    strain_id=protein_sequence.id[0:8]
                    gff_protein_id=protein_sequence.id[9:]
                gene_plain_sequences=globals().get(strain_id+"_db")[gff_protein_id].sequence(
                    str(contig_path/(strain_id+".fasta")),
                    use_strand=True
                    )
                record = SeqRecord(
                    Seq(gene_plain_sequences,Bio.Alphabet.IUPAC.unambiguous_dna),
                    id=protein_sequence.id,
                    description=""
                    )
                SeqIO.write(record, out_fl,"fasta")