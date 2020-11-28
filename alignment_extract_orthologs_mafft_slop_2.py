#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sys import path
from joblib import Parallel,delayed
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio
from rpy2.robjects import NULL
from extract_nucleic_2 import contig_gff3_gffutils_db
from pav_contig_present_for_pan_data_2_blast_gth import multi_alignment
from mafft_gene import extract_gene_gff_mafft
import rpy2.robjects as robjects
from pathlib import Path
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
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
R_code_parse_blast_result='''
parse_blast_result=function(single_copy_gene_vector,multi_copy_gene_Vector,R_blast_vlaue_df,single_err_file_name){
  require(dplyr)
  single_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% single_copy_gene_vector)
  Outliers_stats_df=boxplot.stats(single_blast_vlaue_df$blast_value)
  Outliers_list=Outliers_stats_df$out
  upper=Outliers_stats_df$stats[5]
  lower=Outliers_stats_df$stats[1]
  if (length(Outliers_list)>0){
    write.table(
      dplyr::filter(R_blast_vlaue_df,blast_value %in% c(Outliers_list)),
      single_err_file_name,
      sep =  '\t',
      quote = F,
      row.names = F,
      col.names = T
      )
    return()
  }
  else{
    multi_blast_vlaue_df=R_blast_vlaue_df %>% 
      filter(gene_name %in% multi_copy_gene_Vector) %>% 
      mutate(ok = blast_value >= lower & blast_value <= upper)
    if(!(all(multi_blast_vlaue_df$ok))){
      multi_blast_vlaue_df=multi_blast_vlaue_df[multi_blast_vlaue_df$ok,]
    }
    multi_blast_vlaue_df=mutate(multi_blast_vlaue_df,strain_id=strsplit(gene_name,"_")) 
    dsaf=aggregate(blast_value~strain_id,data = multi_blast_vlaue_df,max)
    best_value_df=filter(multi_blast_vlaue_df,blast_vlaue %in% dsaf)
  }
  return(c(single_blast_vlaue_df$gene_name,best_value_df$gene_name))
}
'''
R_parse_blast_result = SignatureTranslatedAnonymousPackage(R_code_parse_blast_result, "R_parse_blast_result")
def extract_gene_gff(orthogroup_file,strain_db_dir_path,contig_path,all_row_gene_fasta_dir):
    global_names=globals()
    for db_file in strain_db_dir_path.iterdir():
        global_names[db_file.stem.strip()+"_db"]=gffutils.FeatureDB(db_file)
    global_names["MGG_db"]=SeqIO.index("../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta","fasta")
    gene_file_path=all_row_gene_fasta_dir/(orthogroup_file.stem+".fasta")
    if gene_file_path.exists() is True:return
    with gene_file_path.open("w") as out_fl:
        for protein_sequence in SeqIO.parse(orthogroup_file,"fasta"):
            protein_id_list=protein_sequence.id.split("_",1)
            strain_id=protein_id_list[0]
            gff_protein_id=protein_id_list[1]
            if strain_id=="MGG":
                MGG_SeqRecord=globals().get("MGG_db")[protein_sequence.id[0:9]]
                MGG_SeqRecord_amend=SeqRecord(
                    MGG_SeqRecord.seq,
                    id=protein_sequence.id[0:9],
                    description=""
                    )
                SeqIO.write(MGG_SeqRecord_amend,out_fl,"fasta")
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
def blastdb (in_file,db_file):
    make_db_cmd=NcbimakeblastdbCommandline(
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.1+/bin/makeblastdb',
        dbtype='nucl',
        input_file=in_file,
        out=db_file
    )
    make_db_cmd()
def blast(db_file,query_file,blast_out_xml_file,blast_out_asn_file,blast_out_txt_file):
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

def ortholog_blast(
    sequence_file,
    blast_identity_value_dir,
    no_mgg_fl,
    only_itself_fl,
    length_0_fl,
    read_function,
    write_function,
    blast_db_path,
    blast_out_xml_path,
    blast_out_asn_path,
    blast_out_txt_path,
    rbase,
    ):
    #需要修改read
    # sequence_file=Path("OG0000012")
    blast_out_value_tsv_file=(blast_identity_value_dir/(sequence_file.stem+".tsv")).resolve()
    if blast_out_value_tsv_file.exists() is True:
        return(read_function(
            str(blast_out_value_tsv_file),
            header = True,
            sep = "\t",
            **{'stringsAsFactors': False},
            **{'check.names': False}
        ))
    blast_db_file=blast_db_path/(sequence_file.stem+"_blast_head_db")
    blast_out_xml_file=blast_out_xml_path/(sequence_file.stem+"_blast_out.xml")
    blast_out_asn_file=blast_out_asn_path/(sequence_file.stem+"_blast_out.asn")
    blast_out_txt_file=blast_out_txt_path/(sequence_file.stem+"_blast_out.txt")
    # blast_out_xml_file=Path("../Pan_genome_data_2/ortholog_blast_2/blast_general/blast_out/blast_out_xml/OG0000012_blast_out.xml")
    if blast_out_xml_file.exists() is False:
        blastdb(sequence_file,blast_db_file)
        blast(blast_db_file,sequence_file,blast_out_xml_file,blast_out_asn_file,blast_out_txt_file)
    
    with open(blast_out_xml_file) as xml_fl:
        R_blast_vlaue_list=rbase.list()
        i=1
        for record in NCBIXML.parse(xml_fl):
            MGG_flag=False
            gene_name_list=[]
            blast_value_list=[]
            identity_perscent=0
            gene_name=record.query.split()[0]
            if gene_name[0:3]!="MGG":
                continue
            MGG_flag=True
            if record.alignments:
                if len(record.alignments)>1:
                    for alignment in record.alignments:
                        max_flag=-1
                        if alignment.hit_def==gene_name:continue
                        for hsp in alignment.hsps:
                            if max_flag > -1: break
                            identity_perscent=hsp.identities/hsp.align_length
                            max_flag=max_flag+20
                            blast_value_list.append(identity_perscent)
                        gene_name_list.append(alignment.hit_def)
                    if MGG_flag:
                        if len(gene_name_list)==0:
                            only_itself_fl.write("{}".format(sequence_file.stem))
                            continue
                        MGG_head_id=[gene_name]*len(gene_name_list)
                        R_blast_vlaue_list.rx2[i]=robjects.DataFrame({
                            "gene_name":robjects.StrVector(gene_name_list),
                            "blast_value":robjects.FloatVector(blast_value_list),
                            "MGG_head":robjects.StrVector(MGG_head_id)
                        })
                        i=i+1
                    else:
                        no_mgg_fl.write("{}".format(sequence_file.stem))
                        break
    if rbase.length(R_blast_vlaue_list)[0]==0:
        length_0_fl.write("{}".format(sequence_file.stem))
        with blast_out_value_tsv_file.open('w') as blast_out_value_tsv_fl:
            blast_out_value_tsv_fl.write("NA")
        return "NA"
    R_blast_vlaue_df=rbase.do_call("rbind",R_blast_vlaue_list)
    write_function(R_blast_vlaue_df,**{'file': str(blast_out_value_tsv_file)},**{'append': False},**{'quote': False},**{'sep': "\t"},**{'row.names': False},**{'col.names': True})
    return R_blast_vlaue_df


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

def generate_row_list(df_row,out_dir,out_dir_all):
    head_id=df_row.rx(1,1)[0]
    df_row_iter=iter(df_row[1:])
    multicopy_list_file_name=out_dir / (head_id+"_multicopy_list.txt")
    singlecopy_list_file_name=out_dir/(head_id+"_singlecopy_list.txt")
    all_row_gene_list_file_name=out_dir_all/(head_id+".txt")
    with all_row_gene_list_file_name.open("w") as all_row_gene_list_fl:
        with multicopy_list_file_name.open('w+') as multicopy_list_fl:
            with singlecopy_list_file_name.open('w+') as singlecopy_list_fl:
                for R_cell in df_row_iter:
                    if R_cell[0]=='NA' or (str(R_cell[0])=="NA"):
                        continue
                    body_list=R_cell[0].split()
                    for body_id in body_list:
                        body_id=body_id.strip(",")
                        if len(body_list)>1:
                            multicopy_list_fl.write("{}\n".format(body_id))
                            all_row_gene_list_fl.write("{}\n".format(body_id))
                        else:
                            singlecopy_list_fl.write("{}\n".format(body_id))
                            all_row_gene_list_fl.write("{}\n".format(body_id))

def extract_gene(single_copy_list,all_row_gene_fasta_file,single_copy_fasta_file_path):
    all_row_gene_dic=SeqIO.index(all_row_gene_fasta_file,"fasta")
    with single_copy_fasta_file_path.open('w') as single_copy_fasta_fl:
        for single_copy_gene_id in single_copy_list:
            single_copy_gene_id_list=single_copy_gene_id.split("_",1)
            strain_id=single_copy_gene_id_list[0]
            if strain_id=="MGG":
                MGG_SeqRecord=all_row_gene_dic[single_copy_gene_id]
                MGG_SeqRecord_amend=SeqRecord(
                    MGG_SeqRecord.seq,
                    id="70-15",
                    description=""
                    )
                SeqIO.write(MGG_SeqRecord_amend,single_copy_fasta_fl,"fasta")
            else:
                if strain_id=="WD-3-1":
                    strain_id=single_copy_gene_id[0:8]
                OTHER_SeqRecord=all_row_gene_dic[single_copy_gene_id]
                OTHER_SeqRecord_amend=SeqRecord(
                    OTHER_SeqRecord.seq,
                    id=strain_id,
                    description=""
                    )
                SeqIO.write(OTHER_SeqRecord_amend,single_copy_fasta_fl,"fasta")
def extract_ortholog_gene(strain_db_dir_path,contig_path,orthogroup_sequence_path,ortholog_blast_path_name):
    '''
    input 1: strain_db_dir_path
    input 2: contig_path
    input 3: orthogroup_sequence_path
    output 1: ortholog_blast_path_name
    '''
    general_out_path=Path(ortholog_blast_path_name)
    blast_general=general_out_path / "blast_general"
    if blast_general.exists() is False:
        blast_general.mkdir()
    all_row_gene_dir_path=directory_creater(general_out_path/"all_row_gene")
    all_row_gene_fasta_dir=directory_creater(all_row_gene_dir_path/"fasta")
    all_row_gene_list_dir=directory_creater(all_row_gene_dir_path/"list")
    all_row_gene_list_dir_all=directory_creater(all_row_gene_dir_path/"all_list")
    blast_db_path=blast_general/"blast_db"
    if blast_db_path.exists() is False:
        blast_db_path.mkdir()
    blast_out_path=blast_general/"blast_out"
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
    blast_identity_value_dir=blast_general/"blast_identity_value"
    if blast_identity_value_dir.exists() is False:
        blast_identity_value_dir.mkdir()

    base=importr("base")
    utils=importr("utils")
    R_blast_vlaue_list=base.list()
    # Parallel(n_jobs=12)(delayed(generate_row_list)(orthogroup_file,all_row_gene_list_dir_all) for orthogroup_file in orthogroup_sequence_path.iterdir()
    # Parallel(n_jobs=12)(delayed(extract_gene_gff)(orthogroup_file,strain_db_dir_path,contig_path,all_row_gene_fasta_dir) for orthogroup_file in orthogroup_sequence_path.iterdir())
    R_blast_vlaue_list=base.list()
    i=1
    no_mgg_fl=(general_out_path/"no_mgg.txt").open('w')
    only_itself_fl=(general_out_path/"only_itself.txt").open('w')
    length_0_fl=(general_out_path/"length_0.txt").open('w')
    for all_row_gene_fasta_file in blast_identity_value_dir.iterdir():
        R_blast_vlaue_df=ortholog_blast(
                all_row_gene_fasta_file,
                blast_identity_value_dir,
                only_itself_fl,
                length_0_fl,
                no_mgg_fl,
                utils.read_table,
                utils.write_table,
                blast_db_path,
                blast_out_xml_path,
                blast_out_asn_path,
                blast_out_txt_path,
                base
                )
        if R_blast_vlaue_df=="NA" or R_blast_vlaue_df.rx2("V1")=="NA": continue
        R_blast_vlaue_list.rx2[i]=robjects.DataFrame(R_blast_vlaue_df)
        i=i+1
    R_blast_vlaue_df=base.do_call("rbind",R_blast_vlaue_list)
    utils.write_table(
        R_blast_vlaue_df,
        **{'file': str(general_out_path/"all.txt")},
        **{'append': False},
        **{'quote': False},
        **{'sep': "\t"},
        **{'row.names': False},
        **{'col.names': True}
        )

    # for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)):
    #     df_row=ortholog_joined_df.rx(i, True)
    #     df_row_iter=iter(df_row[1:])
    #     head_id=df_row.rx(1,1)[0]
    #     # df_row_iter=next(next(df_row_iter))
    #     all_row_gene_list=[]
    #     single_copy_gene_list=[]
    #     multi_copy_gene_list=[]
    #     all_row_gene_fasta_file=all_row_gene_fasta_dir / (head_id+"_"+"_all_row_gene_fasta.fasta")
    #     generate_row_list(df_row_iter,head_id,all_row_gene_list,single_copy_gene_list,multi_copy_gene_list,all_row_gene_list_dir,na_count,mgg_list)
    #     # if na_count<140:
    #     #     break
    #     extract_gene_gff(all_row_gene_list,strain_db_dir_path,contig_path,all_row_gene_fasta_file)
    #     R_blast_vlaue_df=ortholog_blast(
    #         head_id,
    #         all_row_gene_fasta_file,
    #         blast_identity_value_dir/(head_id+".tsv"),
    #         utils.read_table,
    #         utils.write_table,
    #         blast_db_path,
    #         blast_out_xml_path,
    #         blast_out_asn_path,
    #         blast_out_txt_path,
    #         base
    #         )
    #     R_blast_vlaue_list.rx2[i]=robjects.DataFrame(R_blast_vlaue_df)
    #     i=i+1
    #     #检查结果不是all true？multi_copy_gene_Vector中是否有异常值，blast结果中检查只用一个hsp行不行
    #     # best_blast_value_list=list(R_parse_blast_result.parse_blast_result(
    #     #     single_copy_gene_list,
    #     #     multi_copy_gene_list,
    #     #     R_blast_vlaue_df,
    #     #     str(blast_general/(head_id+"_single_err.txt"))
    #     #     ))
    #     # extract_gene(best_blast_value_list,all_row_gene_fasta_file,blast_general/(head_id+"_single_copy.fasta"))
    # R_blast_vlaue_df=base.do_call("rbind",R_blast_vlaue_list)
    # utils.write_table(R_blast_vlaue_df,**{'file': str(general_out_path/"all.txt")},**{'append': False},**{'quote': False},**{'sep': "\t"},**{'row.names': False},**{'col.names': True})


