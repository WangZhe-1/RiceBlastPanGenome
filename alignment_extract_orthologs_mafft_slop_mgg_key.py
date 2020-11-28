#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from joblib import Memory
import subprocess
from joblib import Parallel,delayed
from distributed import Client
import dask
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio
from joblib import memory
import rpy2.robjects as robjects
from pathlib import Path
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import re
import rpy2
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastformatterCommandline
from Bio.Blast import NCBIXML
from Directory_creater import directory_creater
import json
import re
import gffutils

R_code_read_ortholog_df='''
read_ortholog_df_to_list=function(ortholog_file_name,strain_name,i,ortholog_df_list){
    ortholog_df=read.table(ortholog_file_name,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    ortholog_df_list[[i]]=ortholog_df
    names(ortholog_df_list)[i]=strain_name
    print(colnames(ortholog_df))
  }
'''
R_code_parse_blast_result=r'''
require(dplyr)
extract_max=function(R_blast_vlaue_df,in_vector){
  in_df=filter(R_blast_vlaue_df,gene_name %in% in_vector)
  return(in_df[which.max(in_df$blast_value),]$gene_name)
}


parse_blast_result_have_multi=function(
  single_vector_file_name,
  multi_Vector_file_name,
  in_df_file_name,
  lower_0.95_file_name,
  best_sequence_id_list_file_name
){
  R_blast_vlaue_df=read.table(in_df_file_name,header = T,stringsAsFactors = F)
  head_id=R_blast_vlaue_df$MGG_head[1]
  if(colnames(R_blast_vlaue_df)[1]=="NA."){
    print(c("have_multi_原始输入是NA",head_id))
    return(0)
    }
  single_copy_gene_vector=read.table(single_vector_file_name,stringsAsFactors = F)
  multi_copy_gene_Vector=read.table(multi_Vector_file_name,stringsAsFactors = F)
  R_blast_vlaue_df_lower=R_blast_vlaue_df %>% 
    filter(blast_value<0.95)
  R_blast_vlaue_df=R_blast_vlaue_df %>% 
    filter(blast_value>0.95)
  if(nrow(R_blast_vlaue_df_lower)>0){
    write.table(R_blast_vlaue_df_lower,lower_0.95_file_name,quote = F,sep = '\t',row.names = F)
  }
  if(nrow(R_blast_vlaue_df)==0){
    print(c("have_multi_被过滤为NA",head_id))
    return(0)
    }
  single_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% single_copy_gene_vector$V1)
  multi_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% multi_copy_gene_Vector$V1)
  if(nrow(multi_blast_vlaue_df)==0){
    print(c("have_multi_multi文件为NA",R_blast_vlaue_df$MGG_head[1]))
    print(multi_copy_gene_Vector)
    return(parse_blast_result_no_multi(
      single_vector_file_name,
      in_df_file_name,
      lower_0.95_file_name,
      best_sequence_id_list_file_name
    ))
  }
  multi_blast_vlaue_df=mutate(multi_blast_vlaue_df,strain_id=sub("^(.+?)_.+","\\1",gene_name))
  dsaf=aggregate(gene_name~strain_id,data = multi_blast_vlaue_df,function(x) extract_max(R_blast_vlaue_df,x))
  multi_best_value_df=filter(multi_blast_vlaue_df,gene_name %in% dsaf$gene_name)
  multi_best_value_df=multi_best_value_df[,-4]
  best_value_df=rbind(single_blast_vlaue_df,multi_best_value_df)
  write_gene_id=c(best_value_df$gene_name,R_blast_vlaue_df$MGG_head[1])
  write.table(
    write_gene_id,
    best_sequence_id_list_file_name,
    quote = F,
    row.names = F,
    col.names = F
  )
  return(nrow(best_value_df))
}
parse_blast_result_no_multi=function(
  single_vector_file_name,
  in_df_file_name,
  lower_0.95_file_name,
  best_sequence_id_list_file_name
){
  R_blast_vlaue_df=read.table(in_df_file_name,header = T,stringsAsFactors = F)
  head_id=R_blast_vlaue_df$MGG_head[1]
  if(colnames(R_blast_vlaue_df)[1]=="NA."){
    print(c("single_原始输入是NA",head_id))
    return(0)
    }
  single_copy_gene_vector=read.table(single_vector_file_name,stringsAsFactors = F)
  R_blast_vlaue_df_lower=R_blast_vlaue_df %>% 
    filter(blast_value<0.95)
  R_blast_vlaue_df=R_blast_vlaue_df %>% 
    filter(blast_value>0.95)
  if(nrow(R_blast_vlaue_df_lower)>0){
    write.table(R_blast_vlaue_df_lower,lower_0.95_file_name,quote = F,sep = '\t',row.names = F)
  }
  if(nrow(R_blast_vlaue_df)==0){
    print(c("single_被过滤为NA",head_id))
    return(0)
  }
  write_gene_id=c(R_blast_vlaue_df$gene_name,R_blast_vlaue_df$MGG_head[1])
  write.table(
    write_gene_id,
    best_sequence_id_list_file_name,
    quote = F,
    row.names = F,
    col.names = F
  )
  return(nrow(R_blast_vlaue_df))
}
'''
# R_parse_blast_result = SignatureTranslatedAnonymousPackage(R_code_parse_blast_result, "R_parse_blast_result")
def extract_gene_gff(all_row_gene_list_file,strain_db_dir_path,contig_path,all_row_gene_fasta_dir):
    gene_file_path=all_row_gene_fasta_dir/(all_row_gene_list_file.stem+".fasta")
    if gene_file_path.exists() is True:return
    global_names=globals()
    for db_file in strain_db_dir_path.iterdir():
        global_names[db_file.stem.strip()+"_db"]=gffutils.FeatureDB(db_file)
    global_names["MGG_db"]=SeqIO.index("../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta","fasta")
    with gene_file_path.open("w") as out_fl:
        with all_row_gene_list_file.open() as protein_id_list:
            for protein_id in protein_id_list:
                protein_id=protein_id.strip()
                protein_id_list=protein_id.split("_",1)
                strain_id=protein_id_list[0]
                gff_protein_id=protein_id_list[1]
                if strain_id=="MGG":
                    MGG_SeqRecord=globals().get("MGG_db")[protein_id[0:9]]
                    MGG_SeqRecord_amend=SeqRecord(
                        MGG_SeqRecord.seq,
                        id=protein_id[0:9],
                        description=""
                        )
                    SeqIO.write(MGG_SeqRecord_amend,out_fl,"fasta")
                else:
                    if strain_id=="WD-3-1":
                        strain_id=protein_id[0:8]
                        gff_protein_id=protein_id[9:]
                    gene_plain_sequences=globals().get(strain_id+"_db")[gff_protein_id].sequence(
                        str(contig_path/(strain_id+".fasta")),
                        use_strand=True
                        )
                    record = SeqRecord(
                        Seq(gene_plain_sequences,Bio.Alphabet.IUPAC.unambiguous_dna),
                        id=protein_id,
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
    # all_row_gene_fasta_dir,
    blast_identity_value_dir,
    only_itself_list,
    length_0_list,
    no_mgg_list,
    read_function,
    write_function,
    blast_db_path,
    blast_out_xml_path,
    blast_out_asn_path,
    blast_out_txt_path,
    rlist,
    rlength,
    rdo_call
    ):
    #需要修改read
    # sequence_file=Path("OG0000012")
    blast_out_value_tsv_file=(blast_identity_value_dir/(sequence_file.stem+".tsv")).resolve()
    print(str(blast_out_value_tsv_file))
    # sequence_file=(all_row_gene_fasta_dir/(sequence_file.stem+".fasta")).resolve()
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
        R_blast_vlaue_list=rlist()
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
                            only_itself_list.append(sequence_file.stem)
                            continue
                        MGG_head_id=[gene_name]*len(gene_name_list)
                        R_blast_vlaue_list.rx2[i]=robjects.DataFrame({
                            "gene_name":robjects.StrVector(gene_name_list),
                            "blast_value":robjects.FloatVector(blast_value_list),
                            "MGG_head":robjects.StrVector(MGG_head_id)
                        })
                        i=i+1
                    else:
                        no_mgg_list.append(sequence_file.stem)
                        break
    if rlength(R_blast_vlaue_list)[0]==0:
        length_0_list.append(sequence_file.stem)
        with blast_out_value_tsv_file.open('w') as blast_out_value_tsv_fl:
            blast_out_value_tsv_fl.write("NA")
        return "NA"
    R_blast_vlaue_df=rdo_call("rbind",R_blast_vlaue_list)
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

def generate_row_list(
    df_row,
    out_dir,
    out_dir_all,
    head_id_list,
    strain_num_list
    ):
    df_row_iter=iter(df_row[2:])
    head_id=df_row.rx(1,2)[0]
    strain_count=1
    if re.search(",",head_id) is not None:
        return
    multicopy_list_file_name=out_dir / (head_id+"_multicopy_list.txt")
    singlecopy_list_file_name=out_dir/(head_id+"_singlecopy_list.txt")
    all_row_gene_list_file_name=out_dir_all/(head_id+".txt")
    with all_row_gene_list_file_name.open("w") as all_row_gene_list_fl:
        with multicopy_list_file_name.open('w+') as multicopy_list_fl:
            with singlecopy_list_file_name.open('w+') as singlecopy_list_fl:
                for R_cell in df_row_iter:
                    if R_cell[0]=='NA' or (str(R_cell[0])=="NA"):
                        continue
                    strain_count=strain_count+1
                    body_list=R_cell[0].split()
                    for body_id in body_list:
                        body_id=body_id.strip(",")
                        if len(body_list)>1:
                            multicopy_list_fl.write("{}\n".format(body_id))
                            all_row_gene_list_fl.write("{}\n".format(body_id))
                        else:
                            singlecopy_list_fl.write("{}\n".format(body_id))
                            all_row_gene_list_fl.write("{}\n".format(body_id))
                all_row_gene_list_fl.write("{}\n".format(head_id))
                head_id_list.append(head_id)
                strain_num_list.append(strain_count)


def extract_gene(single_copy_list,all_row_gene_fasta_file,single_copy_fasta_file_path):
    if single_copy_fasta_file_path.exists() is True: return
    all_row_gene_dic=SeqIO.index(str(all_row_gene_fasta_file),"fasta")
    with single_copy_fasta_file_path.open('w') as single_copy_fasta_fl:
        with single_copy_list.open() as single_copy_list_fl:
            for single_copy_gene_id in single_copy_list_fl:
                single_copy_gene_id=single_copy_gene_id.strip()
                single_copy_gene_id_list=single_copy_gene_id.split("_",1)
                strain_id=single_copy_gene_id_list[0]
                if strain_id=="MGG":
                    MGG_SeqRecord=all_row_gene_dic[single_copy_gene_id[0:9]]
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
def merge_df(input_list,rlist,rwrite_table,rdo_call,general_out_path):
    i=1
    R_blast_vlaue_list=rlist()
    for input_df in input_list:
        if input_df=="NA" or input_df.rx2("V1")=="NA": continue
        R_blast_vlaue_list.rx2[i]=robjects.DataFrame(input_df)
        i=i+1
    R_blast_vlaue_df=rdo_call("rbind",R_blast_vlaue_list)
    rwrite_table(
        R_blast_vlaue_df,
        **{'file': str(general_out_path/"all.txt")},
        **{'append': False},
        **{'quote': False},
        **{'sep': "\t"},
        **{'row.names': False},
        **{'col.names': True}
        )
def parse_blast_result(
    single_vector_file_name,
    multi_Vector_file_name,
    in_df_file_name,
    lower_file_name,
    best_sequence_id_list_file_name
):

    if multi_Vector_file_name.stat().st_size == 0: 
        strain_num=R_parse_blast_result.parse_blast_result_no_multi(
            str(single_vector_file_name),
            str(in_df_file_name),
            str(lower_file_name),
            str(best_sequence_id_list_file_name)
            )
        # print (in_df_file_name.stem)
        return in_df_file_name.stem,strain_num[0]
    strain_num=R_parse_blast_result.parse_blast_result_have_multi(
        str(single_vector_file_name),
        str(multi_Vector_file_name),
        str(in_df_file_name),
        str(lower_file_name),
        str(best_sequence_id_list_file_name)
    )
    # print (in_df_file_name.stem)
    return in_df_file_name.stem,strain_num[0]
def extract_ortholog_gene(strain_db_dir_path,contig_path,joined_df_file_name,ortholog_blast_path_name):
    '''
    input 1: strain_db_dir_path
    input 2: contig_path
    input 3: joined_df_file_name
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
    lower_blast_value_dir_path=directory_creater(general_out_path/"lower_blast_value")
    best_sequence_id_list_dir_path=directory_creater(general_out_path/"best_sequence_id_list")
    single_copy_fasta_dir_path=directory_creater(general_out_path.parent/"MAFFT_ortholog_MGG"/"in_put_fasta")
    base=importr("base")
    utils=importr("utils")
    
    # na_count=0
    # R_blast_vlaue_list=base.list()
    # # client = Client(processes=False)
    
    # generate_row_list_calls=[]
    # extract_gene_gff_calls=[]
    # for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)):
    #     df_row=ortholog_joined_df.rx(i, True)
    #     df_row_iter=iter(df_row[2:])
    #     head_id=df_row.rx(1,2)[0]
    #     if re.search(",",head_id) is not None:
    #         continue
    #     all_row_gene_fasta_file=all_row_gene_fasta_dir / (head_id+"_"+"_all_row_gene_fasta.fasta")
        # generate_row_list_calls.append(dask.delayed(generate_row_list)(
        #     df_row_iter,
        #     head_id,
        #     all_row_gene_list_dir,
        #     all_row_gene_list_dir_all
        # ))
        # extract_gene_gff_calls.append(dask.delayed(extract_gene_gff)(
        #     all_row_gene_list_dir_all/(head_id+".txt"),
        #     strain_db_dir_path,
        #     contig_path,
        #     all_row_gene_fasta_dir
        # ))
    # dask.compute(*generate_row_list_calls)
    # dask.compute(*extract_gene_gff_calls)
    '''
    generate_row_list
    '''
    # ortholog_joined_df=utils.read_table(
    #     str(joined_df_file_name),
    #     sep = "\t",
    #     header = True,
    #     **{'stringsAsFactors': False},
    #     **{'check.names': False}
    #     )
    # head_id_list=[]
    # strain_num_list=[]
    # Parallel(n_jobs=1)(delayed(generate_row_list)(
    #     ortholog_joined_df.rx(i, True),
    #     all_row_gene_list_dir,
    #     all_row_gene_list_dir_all,
    #     head_id_list,
    #     strain_num_list
    #     ) for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)))
    # with (general_out_path/"MGG_strain_num.txt").open('w') as MGG_strain_num_fl:
    #     for head_id_1,strain_num_1 in zip(head_id_list,strain_num_list):
    #         MGG_strain_num_fl.write("{}\t{}\n".format(head_id_1,strain_num_1))
    '''
    extract_gene_gff
    '''
    # Parallel(n_jobs=12)(delayed(extract_gene_gff)(all_row_gene_list_file,strain_db_dir_path,contig_path,all_row_gene_fasta_dir) for all_row_gene_list_file in all_row_gene_list_dir_all.iterdir())
    # Parallel(n_jobs=1)(delayed(extract_gene_gff)(all_row_gene_list_file,strain_db_dir_path,contig_path,all_row_gene_fasta_dir) for all_row_gene_list_file in all_row_gene_list_dir_all.glob("MGG_16565T0*"))
    '''
    ortholog_blast
    '''
    # output_list=[]
    # lineList=[]
    # rdocall=base.do_call
    # rread_table=utils.read_table
    # rwrite_table=utils.write_table
    # rlist=base.list
    # rlength=base.length
    # only_itself_list=[]
    # length_0_list=[]
    # no_mgg_list=[]
    # for all_row_gene_fasta_file in all_row_gene_fasta_dir.glob("MGG_01742T0*"):
    # # for all_row_gene_fasta_file in all_row_gene_fasta_dir.iterdir():
    #     R_blast_vlaue_df=dask.delayed(ortholog_blast)(
    #     # R_blast_vlaue_df=ortholog_blast(
    #         all_row_gene_fasta_file,
    #         # all_row_gene_fasta_dir,
    #         blast_identity_value_dir,
    #         only_itself_list,
    #         length_0_list,
    #         no_mgg_list,
    #         rread_table,
    #         rwrite_table,
    #         blast_db_path,
    #         blast_out_xml_path,
    #         blast_out_asn_path,
    #         blast_out_txt_path,
    #         rlist,
    #         rlength,
    #         rdocall
    #         )
    #     output_list.append(R_blast_vlaue_df)
    # total=dask.delayed(merge_df)(output_list,base.list,utils.write_table,base.do_call,general_out_path)
    # # total.visualize()
    # total.compute()
    '''
    parse_blast_result
    '''
    # mem=Memory(general_out_path/"parse_blast_result_cache")
    # parse_blast_result_mem=mem.cache(parse_blast_result,verbose=0)
    # reddd=Parallel(n_jobs=10)(delayed(parse_blast_result_mem)(
    #     all_row_gene_list_dir/(blast_identity_value_tsv.stem+"_singlecopy_list.txt"),
    #     all_row_gene_list_dir/(blast_identity_value_tsv.stem+"_multicopy_list.txt"),
    #     blast_identity_value_tsv,
    #     lower_blast_value_dir_path/(blast_identity_value_tsv.stem+".txt"),
    #     best_sequence_id_list_dir_path/(blast_identity_value_tsv.stem+".txt")
    # ) for blast_identity_value_tsv in blast_identity_value_dir.iterdir())
    # # ) for blast_identity_value_tsv in blast_identity_value_dir.glob("MGG_00010T0*"))
    # with (general_out_path/"best_strain_num.txt").open('w') as best_strain_num_fl:
    #     for head_id_1,strain_num_1 in reddd:
    #         best_strain_num_fl.write("{}\t{}\n".format(head_id_1,strain_num_1))

    '''
    extract_gene
    '''
    Parallel(n_jobs=12)(delayed(extract_gene)(
        best_sequence_id_list,
        all_row_gene_fasta_dir/(best_sequence_id_list.stem+".fasta"),
        single_copy_fasta_dir_path/(best_sequence_id_list.stem+".fasta")
    ) for best_sequence_id_list in best_sequence_id_list_dir_path.iterdir())


    # for blast_identity_value_tsv in blast_identity_value_dir.iterdir():
    #     parse_blast_result(
    #         all_row_gene_list_dir/(blast_identity_value_tsv.stem+"_singlecopy_list.txt"),
    #         all_row_gene_list_dir/(blast_identity_value_tsv.stem+"_multicopy_list.txt"),
    #         blast_identity_value_tsv,
    #         lower_blast_value_dir_path/(blast_identity_value_tsv.stem+".txt"),
    #         best_sequence_id_list_dir_path/(blast_identity_value_tsv.stem+".txt")
    #     )
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
if __name__ == '__main__':
    # client = Client(n_workers=12, threads_per_worker=1)
    extract_ortholog_gene(
        Path("../Pan_genome_data_2/contig_gff3_gffutils_db"),
        Path("../Pan_genome_data_2/156_contig"),
        Path("../Pan_genome_data_2/ortho/joined_df.tsv"),
        directory_creater("../Pan_genome_data_2/ortholog_blast_mgg_key"),
    )
