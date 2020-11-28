#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 把1891序列读进去，拿出核酸id（摘T之前的），核酸读到index结构，用id找对应的核酸

import re
import subprocess
from pathlib import Path
from Bio import SeqIO
import uniformize_gene_to_protein
from extract_strain_id import extract_strain_id
def merge_to_one(GFF_path_name,id_file,list_file,gene_base,cat_err_file):
    '''
    input 1: GFF path
    input 2: strain list
    output 1: list file
    output 2: gene base
    output 3: cat_err_file
    '''
    fasta_list=[
        'cat',
        '/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta',
        "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/New_add_ina168/ina168_gene.fasta"
        ]
    strain_95_list=extract_strain_id(id_file)
    GFF_path=Path(GFF_path_name)
    with open (list_file,'w+') as list_out:
        for species_id in strain_95_list:
            for fasta_file in GFF_path.glob(species_id.strip('\n')+'*_gene.fasta'):
                list_out.write('{}\n'.format(fasta_file))
                fasta_list.append(str(fasta_file.resolve()))
        
    fasta_cat=subprocess.Popen(
        fasta_list,
        stdout=open(gene_base,'w+'),
        stderr=open(cat_err_file,'w+'),
        universal_newlines=True,
    )
    fasta_cat.wait()

def extract_gene(protein_file,gene_base_name,pan_gene_file_name,gene_protein_mapping_table):
    '''
    input 1: pan protein fasta file
    input 2: merged gene file
    input 3: gene_protein_mapping_table
    output 1:pan_gene_file_name
    don't froget remove 
    '''
    gene_base=SeqIO.index(
        gene_base_name,
        "fasta",
        key_function=uniformize_gene_to_protein.get_id_gene
        )
    with open(gene_protein_mapping_table,'w+') as gene_protein_fl:
        with open(pan_gene_file_name,'w+') as out_fl:
            for protein_sequence in SeqIO.parse(protein_file,"fasta"):
                # if re.search("\.t1",protein_sequence.id) is not None: continue
                # if re.search("^70-15",protein_sequence.id) is not None: continue
                gene_id=uniformize_gene_to_protein.get_id_protein(protein_sequence.id)
                gene_sequence=gene_base[gene_id]
                SeqIO.write(gene_sequence,out_fl,"fasta")
                gene_protein_fl.write("{}\t{}\n".format(gene_sequence.id,protein_sequence.id))