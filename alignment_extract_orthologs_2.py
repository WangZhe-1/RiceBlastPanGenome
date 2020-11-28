#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import rpy2.robjects as robjects
from pathlib import Path
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.packages import importr
import re
import rpy2
from uniformize_gene_to_protein import get_id_gene
from uniformize_gene_to_protein import get_id_protein
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastformatterCommandline
from Bio.Blast import NCBIXML
from Directory_creater import directory_creater
import gffutils
import json
from extract_strain_id import extract_strain_id
from gffutils.helpers import asinterval
import pybedtools
import re

R_code_read_ortholog_df='''
read_ortholog_df_to_list=function(ortholog_file_name,strain_name,i,ortholog_df_list){
    ortholog_df=read.table(ortholog_file_name,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    ortholog_df_list[[i]]=ortholog_df
    names(ortholog_df_list)[i]=strain_name
    print(colnames(ortholog_df))
  }
'''
def read_file(ortholog_path_name,joined_df_file_name):
    '''
    input 1: 70-15 ortholog_path_name
    output 1: joined_df_file_name 
    '''
    rgx_strain_name=re.compile("v__(.+)$")    
    base=importr("base")
    R_ortholog_df_list=base.vector("list",2)
    i=1
    utils=importr("utils")
    for ortholog_file in Path(ortholog_path_name).iterdir():
        strain_name=rgx_strain_name.search(ortholog_file.stem).group(1)
        R_ortholog_df_list.rx2[i]=utils.read_table(
            str(ortholog_file),
            header = True,
            sep = "\t",
            **{'stringsAsFactors': False},
            **{'check.names': False}
            )
        i=i+1
    plyr=importr("plyr")
    ortholog_joined_df=plyr.join_all(R_ortholog_df_list,**{'by': "70-15"},**{'type': "full"},**{'match': "all"})
    utils.write_table(ortholog_joined_df,**{'file': joined_df_file_name},**{'append': False},**{'quote': False},**{'sep': "\t"},**{'row.names': False},**{'col.names': True})
gene_base=None
def extract_gene(protein_id,out_fl):
    gene_id=get_id_protein(protein_id)
    gene_sequence=gene_base[gene_id]
    SeqIO.write(gene_sequence,out_fl,"fasta")

def extract_gene_two_head(protein_id):
    gene_id=get_id_protein(protein_id)
    gene_sequence=gene_base[gene_id]
    return(gene_sequence)
def blastdb (in_file,db_file):
    make_db_cmd=NcbimakeblastdbCommandline(
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.1+/bin/makeblastdb',
        dbtype='nucl',
        input_file=in_file,
        out=db_file
    )
    make_db_cmd()
def blast(db_file,query_file):
    blast_asn_cmd=NcbiblastnCommandline(
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.1+/bin/blastn',
        query=query_file,
        db=db_file,
        outfmt=11,
        out=blast_out_asn_file
    )
    blast_asn_cmd()
    blast_xml_cmd=NcbiblastformatterCommandline(
        archive=blast_out_asn_file,
        outfmt=5,
        out=blast_out_xml_file,
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.1+/bin/blast_formatter'
    )
    blast_xml_cmd()
    blast_txt_cmd=NcbiblastformatterCommandline(
        archive=blast_out_asn_file,
        outfmt=7,
        out=blast_out_txt_file,
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.1+/bin/blast_formatter'
    )
    blast_txt_cmd()

def ortholog_blast(head_id,body_id):
    blast_head_sequence_file=blast_head_sequence_path / (head_id+"_blast_head_sequence.fasta")
    blast_db_file=blast_db_path/(head_id+"_blast_head_db")
    if blast_head_sequence_file.is_file() is False:
        with blast_head_sequence_file.open('w+') as blast_head_sequence_fl:
            extract_gene(head_id,blast_head_sequence_fl)
        blastdb(blast_head_sequence_file,blast_db_file)
    blast_body_sequence_file=blast_body_sequence_path / (head_id+"_"+body_id+"_blast_body_sequence.fasta")
    global blast_out_asn_file,blast_out_xml_file,blast_out_txt_file
    blast_out_xml_file=blast_out_xml_path/(head_id+"_"+body_id+"_blast_out.xml")
    blast_out_asn_file=blast_out_asn_path/(head_id+"_"+body_id+"_blast_out.asn")
    blast_out_txt_file=blast_out_txt_path/(head_id+"_"+body_id+"_blast_out.txt")
    if blast_body_sequence_file.is_file() is False:
        with blast_body_sequence_file.open('w+') as blast_body_sequence_fl:
            extract_gene(body_id,blast_body_sequence_fl)
        blast(blast_db_file,blast_body_sequence_file)
    
    with open(blast_out_xml_file) as xml_fl:
        for record in NCBIXML.parse(xml_fl):
            gene_name=record.query.split()[0]
            if record.alignments:
                max_flag=-1
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if max_flag == -1:
                            identity_discriminant_for_length=hsp.align_length/record.query_length
                            identity_discriminant_for_identity_perscent=hsp.identities/hsp.align_length
                            max_flag=max_flag+2
                            return identity_discriminant_for_identity_perscent,identity_discriminant_for_length
def generate_bed(gff_feature_item):
    yield asinterval(gff_feature_item)

def slop_list2sequence(head_id,protein_id,out_fl,protein_id_vs_chrom_fl):
    if protein_id[0:3]=="MGG":
        bedtools_sequence=pybedtools.BedTool(generate_bed(MGG_db[protein_id[0:9]]))
        extended_genes_slop = bedtools_sequence.slop(b=1000,genome=globals().get("70-15_json"),s=True)
        sequence=extended_genes_slop.sequence(fi=str(contig_dir_path_name/"70-15.fasta"),s=True)
        slop_sequence=SeqIO.read(sequence.seqfn,"fasta")
        SeqIO.write(slop_sequence,out_fl,"fasta")
        protein_id_vs_chrom_fl.write("{}\t{}\n".format(protein_id,slop_sequence.id))
    else:
        strain_id,protein_ordinal=strain_protein_id_pattern.search(protein_id).group(1,2)
        bedtools_sequence=pybedtools.BedTool(generate_bed(globals().get(strain_id+"_db")["gene_"+protein_ordinal]))
        extended_genes_slop = bedtools_sequence.slop(b=1000,genome=globals().get(strain_id+"_json"),s=True)
        sequence=extended_genes_slop.sequence(fi=str(contig_dir_path_name/(strain_id+".fasta")),s=True)
        slop_sequence=SeqIO.read(sequence.seqfn,"fasta")
        SeqIO.write(slop_sequence,out_fl,"fasta")
        protein_id_vs_chrom_fl.write("{}\t{}\n".format(protein_id,slop_sequence.id))

def slop_list2gff():
    global_names=globals()
    for protein_id in slop_list:
        if protein_id[0:3]=="MGG":
            head_id=protein_id
            yield asinterval(MGG_db[protein_id[0:9]])
        else:
            strain_id,protein_ordinal=strain_protein_id_pattern.search(protein_id).group(1,2)
            yield asinterval(global_names.get(strain_id+"_db")["gene_"+protein_ordinal])

def two_head(head_list,df_row):
    list_fst_sequence_two_head=[]
    for R_cell in df_row:
        if (R_cell[0]=="NA") or (str(R_cell[0])=="NA"):
            continue
        body_list = R_cell[0].split()
        for body_id in body_list:
            body_id=body_id.strip(",")
            # if len(body_list) >1:
            #     body_variable_list_fl.write("{}\n".format(body_id))
            #     blast_general_same_strain_list_fl.write("{}\n".format(body_id))
            list_fst_sequence_two_head.append(extract_gene_two_head(body_id))
    for head_id in head_list:
        fst_sequence_two_head_file=fst_sequence_two_head / (head_id[0:11]+"_fst_sequence.fasta")
        with fst_sequence_two_head_file.open('w+') as fst_sequence_two_head_fl:
            for sequence in list_fst_sequence_two_head:
                SeqIO.write(sequence,fst_sequence_two_head_fl,"fasta")
            SeqIO.write(extract_gene_two_head(head_id[0:11]),fst_sequence_two_head_fl,"fasta")
def one2one(head_id,df_row):
    fst_sequence_file=fst_sequence / (head_id+"_fst_sequence.fasta")
    body_variable_list_file_name=blast_per_MGG_gene / (head_id+"_body_variable_list.txt")
    global body_blast_value_list_fl
    body_blast_value_list_file_name=blast_per_MGG_gene/(head_id+"_blast_value_list.txt")
    global slop_list
    slop_list=[]
    fst_sequence_slop_1K_file=fst_sequence_slop_1K_dir/(head_id+"_slop_1K.fasta")
    protein_id_vs_chrom_file=protein_id_vs_chrom_dir/(head_id+"_protein_id_vs_chrom.txt")
    with protein_id_vs_chrom_file.open("w+") as protein_id_vs_chrom_fi:
        with fst_sequence_slop_1K_file.open("w+") as fst_sequence_slop_1K_fl:
            with body_variable_list_file_name.open('w+') as body_variable_list_fl:
                with fst_sequence_file.open('w+') as fst_sequence_fl:
                    with body_blast_value_list_file_name.open('w+') as body_blast_value_list_fl:
                        extract_gene(head_id,fst_sequence_fl)
                        slop_list.append(head_id)
                        slop_list2sequence(head_id,head_id,fst_sequence_slop_1K_fl,protein_id_vs_chrom_fi)
                        body_blast_value_dic={}
                        for R_cell in df_row:
                            if R_cell[0]=='NA' or (str(R_cell[0])=="NA"):
                                continue
                            body_list=R_cell[0].split()
                            for body_id in body_list:
                                body_id=body_id.strip(",")
                                if len(body_list)>1:
                                    body_variable_list_fl.write("{}\n".format(body_id))
                                    blast_general_same_strain_list_fl.write("{}\n".format(body_id))
                                body_blast_value_dic.setdefault(body_id,ortholog_blast(head_id,body_id))
                            if (None in body_blast_value_dic.values()):
                                for item in body_blast_value_dic:
                                    item_value=body_blast_value_dic[item]
                                    if item_value==None:
                                        body_blast_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item+"*",0,0))
                                        blast_general_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item+"*",0,0))
                                    else:
                                        body_blast_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item+"*",item_value[0],item_value[1]))
                                        blast_general_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item+"*",item_value[0],item_value[1]))
                            else:
                                for item in body_blast_value_dic:
                                    item_value=body_blast_value_dic[item]
                                    body_blast_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item,item_value[0],item_value[1]))
                                    blast_general_value_list_fl.write("{}\t{}\t{}\t{}\n".format(head_id,item,item_value[0],item_value[1]))
                                max_blast_value_protein_id=max(body_blast_value_dic.items(), key=lambda item: item[1][0])[0]
                                extract_gene(max_blast_value_protein_id,fst_sequence_fl)
                                slop_list.append(max_blast_value_protein_id)
                                slop_list2sequence(head_id,max_blast_value_protein_id,fst_sequence_slop_1K_fl,protein_id_vs_chrom_fi)
                            body_blast_value_dic.clear()

def extract_ortholog_gene(gene_base_name,joined_df_file_name,id_file,ortholog_blast_path_name):
    '''
    input 1: gene_base_name
    input 2: joined_df_file_name
    input 3: id_file
    output 1: ortholog_blast_path_name
    '''
    global gene_base
    gene_base=SeqIO.index(
        gene_base_name,
        "fasta",
        key_function=get_id_gene
        )
    blast_path=Path(ortholog_blast_path_name)
    global blast_general,blast_per_MGG_gene,blast_body_sequence_path,blast_head_sequence_path,blast_db_path,blast_out_path,blast_out_xml_path,blast_out_txt_path,blast_sequence,fst_sequence,blast_result,blast_out_asn_path
    sequence_out_dir=directory_creater(blast_path.parent/"sequence_out")
    blast_general=blast_path / "blast_general"
    if blast_general.exists() is False:
        blast_general.mkdir()
    blast_per_MGG_gene=blast_path/"blast_per_MGG_gene"
    if blast_per_MGG_gene.exists() is False:
        blast_per_MGG_gene.mkdir()
    blast_body_sequence_path=blast_path/"blast_body_sequence"
    if blast_body_sequence_path.exists() is False:
        blast_body_sequence_path.mkdir()
    blast_head_sequence_path=blast_path/"blast_head_sequence"
    if blast_head_sequence_path.exists() is False:
        blast_head_sequence_path.mkdir()
    blast_db_path=blast_path/"blast_db"
    if blast_db_path.exists() is False:
        blast_db_path.mkdir()
    blast_out_path=blast_path/"blast_out"
    if blast_out_path.exists() is False:
        blast_out_path.mkdir()
    blast_out_xml_path=blast_out_path/"blast_out_xml"
    if blast_out_xml_path.exists() is False:
        blast_out_xml_path.mkdir()
    blast_out_asn_path=blast_out_path/"blast_out_asn"
    if blast_out_asn_path.exists() is False:
        blast_out_asn_path.mkdir()
    blast_out_txt_path=blast_out_path/"blast_out_txt"
    if blast_out_txt_path.exists() is False:
        blast_out_txt_path.mkdir()
    
    blast_sequence=blast_path/"blast_sequence"
    if blast_sequence.exists() is False:
        blast_sequence.mkdir()
    fst_sequence=sequence_out_dir/"fst_sequence"
    if fst_sequence.exists() is False:
        fst_sequence.mkdir()
    global fst_sequence_two_head
    fst_sequence_two_head=sequence_out_dir/"fst_sequence_two_head"
    if fst_sequence_two_head.exists() is False:
        fst_sequence_two_head.mkdir()
    blast_result=blast_path/"blast_result"
    if blast_result.exists() is False:
        blast_result.mkdir()
    global blast_exception_path_name
    blast_exception_path_name=directory_creater(blast_path/"blast_exception")
    global gffutils_db_dir_path_name
    gffutils_db_dir_path_name=directory_creater(blast_path.parent/"gffutils_db")
    global json_dir_path_name
    json_dir_path_name=Path("../Pan_genome_data/contig_length_json/")
    global contig_dir_path_name
    contig_dir_path_name=Path("../../contig/")
    global gff_path_name
    gff_path_name=Path("../../GFF/")
    
    
    gff_out_path=directory_creater(sequence_out_dir/"gene_slop_1K_gff")
    global fst_sequence_slop_1K_dir
    fst_sequence_slop_1K_dir=directory_creater(sequence_out_dir/"fst_sequence_slop_1K")
    global protein_id_vs_chrom_dir
    protein_id_vs_chrom_dir=directory_creater(sequence_out_dir/"protein_id_vs_chrom")

    species_95_list=extract_strain_id(id_file)
    species_95_list.append("ina168")
    species_95_list.remove("magnaporthe_oryzae_70-15_8_proteins_T0")
    global_names=globals()
    for json_file in json_dir_path_name.iterdir():
        with json_file.open() as json_fl:
            global_names[json_file.stem+"_json"]=json.load(json_fl)
    for strain_id in species_95_list:
        for gff_file in gff_path_name.glob(strain_id+".gff"):
            gff_db_name=gffutils_db_dir_path_name/(strain_id+".db")
            if gff_db_name.is_file() is False:
                gffutils.create_db(str(gff_file),str(gff_db_name),force=True,id_spec=None)
            global_names[strain_id+"_db"]=gffutils.FeatureDB(gff_db_name)


    MGG_db_file_path=gffutils_db_dir_path_name/"MGGdb.db"
    if MGG_db_file_path.is_file() is False:
        gffutils.create_db("../../70-15_refference_genome/70-15_Gff/magnaporthe_oryzae_70-15_8_genome_summary_per_gene_amend.txt",str(MGG_db_file_path),id_spec=':source:')
    global MGG_db
    MGG_db=gffutils.FeatureDB(MGG_db_file_path)
    global strain_protein_id_pattern
    strain_protein_id_pattern=re.compile("(.+)_protein_(.+)_")

    base=importr("base")
    utils=importr("utils")
    ortholog_joined_df=utils.read_table(
        joined_df_file_name,
        sep = "\t",
        header = True,
        **{'stringsAsFactors': False},
        **{'check.names': False}
        )
    blast_general_same_strain_list=blast_general/"same_strain_list.txt"
    blast_general_value_list=blast_general/"value_list.txt"
    global blast_general_same_strain_list_fl,blast_general_value_list_fl
    ortholog_joined_df_sub=ortholog_joined_df.rx(True,-1)
    with blast_general_same_strain_list.open('w+') as blast_general_same_strain_list_fl:
        with blast_general_value_list.open('w+') as blast_general_value_list_fl:
            for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)):
                df_row=ortholog_joined_df_sub.rx(i, True)
                df_row_iter=iter(df_row)
                head_list=next(df_row_iter)[0].split()
                if len(head_list)==1:
                    one2one(head_list[0],df_row_iter)
                    bedtool_file = pybedtools.BedTool(slop_list2gff()).saveas(gff_out_path/(head_list[0]+".gff"))
                else:
                    two_head(head_list,df_row_iter)

            
        





