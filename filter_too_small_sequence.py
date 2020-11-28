from __future__ import with_statement
from qc_contig_read_busco import out_fl
from Bio.Blast.Record import Alignment
from pathlib import Path
from Bio.Blast import NCBIXML
from Bio import SeqIO
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
def filter_pan_id(pan_id_file_name,length_table_file,filtered_pan_id_file_name):
    '''
    input 1: pan_id_file_name
    input 2: length_table_file
    output 1: filtered_pan_id_file_name
    '''
    R_code='''
    filter_pan_id=function(pan_id_file_name,length_table_file,filtered_pan_id_file_name){
    require(dplyr)
    pan_id=read.table(pan_id_file_name)
    length_table=read.table(length_table_file)
    length_table_filter=length_table %>% 
        filter(V2>20)
    pan_id_filter=merge(pan_id,length_table_filter,by.x = 2,by.y = 1,all.y = T)
    pan_id_filter=pan_id_filter[,-3]
    pan_id_filter=pan_id_filter[,c(2,1)]
    write.table(pan_id_filter,filtered_pan_id_file_name,sep = "\t",quote = F,row.names = F,col.names = F)
    }
    '''
    R_filter_pan_id=SignatureTranslatedAnonymousPackage(R_code,"R_filter_pan_id")
    R_filter_pan_id.filter_pan_id(pan_id_file_name,length_table_file,filtered_pan_id_file_name)
def filter_blast(blast_result_file_name,pan_protein_file_path,filtered_pan_protein_file_path,length_table):
    '''
    input 1: blast_result_file_name
    input 2: pan_protein_file_path
    output 1: filtered_pan_protein_file_path
    output 2: length_table
    '''
    pan_protein_dic=SeqIO.index(str(pan_protein_file_path),"fasta")
    with filtered_pan_protein_file_path.open("w+") as out_fl:
        with length_table.open("w") as length_table_fl:
            for record in NCBIXML.parse(open(blast_result_file_name)):
                query_name=record.query
                query_len=record.query_letters
                assert record.query_letters==record.query_length
                length_table_fl.write("{}\t{}\n".format(query_name,query_len))
                if query_len<=20: continue
                SeqIO.write(pan_protein_dic[query_name],out_fl,"fasta")
                # if record.alignments:
                #     max_flag=-1
                #     for Alignment in record.alignments:
                #         if max_flag==-1:
                #             max_flag=max_flag+2