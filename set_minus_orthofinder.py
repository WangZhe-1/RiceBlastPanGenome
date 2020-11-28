'''
@Author: your name
@Date: 2020-08-12 10:24:30
@LastEditTime: 2020-08-12 10:24:30
@LastEditors: your name
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/set_minus_orthofinder.py
'''

import rpy2.robjects as robjects
R_pav_orthofinder_code="""
library(tidyverse)
library(openxlsx)

# 读orthofinder的结果数据
# tsv_in_df=read.csv("..\\protein_host_rice_low_quality_removed\\OrthoFinder\\Results_Jan14\\WorkingDirectory\\OrthoFinder_old\\Results_Feb03_1\\Orthogroups\\Orthogroups.tsv",sep = "\t",nrows = 10,check.names = F,stringsAsFactors = F,na.strings = "")
tsv_assignedGenes_in_df=read.csv(orthofinder_tsv_file_name,sep = "\t",check.names = F,stringsAsFactors = F,na.strings = "")
tsv_UnassignedGenes_in_df=read.csv(orthofinder_unassianed_tsv_file_name,sep = "\t",check.names = F,stringsAsFactors = F,na.strings = "")
# tsv_in_df=read.csv("..\\Results_Feb17_未成功\\Orthogroups\\Orthogroups.tsv",sep = "\t",check.names = F,stringsAsFactors = F,na.strings = "")
pan_id=read.table(pan_id_file,sep="\t",stringsAsFactors = F)
# 把两部分合并起来
tsv_all_in_df=rbind(tsv_assignedGenes_in_df,tsv_UnassignedGenes_in_df)


# 用于把70-15放到前面
# col_7015=grep("70-15.+",colnames(tsv_in_df))
# pav_df=tsv_in_df[,c(col_7015,1:(col_7015-1),(col_7015+1):ncol(tsv_in_df))]

# 把菌株名的_protein去掉，会损失"_"后面的内容
colnames(tsv_all_in_df)=sub("_.+$","",colnames(tsv_all_in_df))
rownames(tsv_all_in_df)=tsv_all_in_df$Orthogroup
pav_df=tsv_all_in_df[,-1]
colnames(pav_df)[colnames(pav_df) == 'WD-3-1'] <- 'WD-3-1_1'

na_tsv=is.na(pav_df)
pav_df[na_tsv]=0
pav_df[!na_tsv]=1


pav_df=data.frame(pav_df,check.names = F,stringsAsFactors = F)
pav_df=pav_df %>% 
  rownames_to_column("Orthogroup")
# 删掉MGG_06302 MGG_06131 MGG_16639，他们三个有两个或以上的转录本，保留在各个菌株分布范围最广的那个

# 把基因id和og号整合好放进来
pav_df=merge(pav_df,pan_id,by.x = 1,by.y = 1,all.y = T)
names(pav_df)[names(pav_df) == 'V2'] <- "protein_id"
WriteXLS::WriteXLS(
    pav_df,
    paste(result_path,"pav_orthofinder.xlsx",sep = ""),
    col.names = T,
    row.names = T
  )

Cut=function(start_point){
  gene_is=pav_df %>% 
    filter((!!sym(start_point))==1) %>%
    # filter(`70-15`)==1
    select("protein_id") 
  
  minus_part=pav_df %>% 
    filter((!!sym(start_point))==1) %>% 
    column_to_rownames("protein_id")
  # pav_df_colsum=colSums(minus_part_num)
  # pav_df_colsum_sort=sort(pav_df_colsum)
  add_part=pav_df %>% 
    filter((!!sym(start_point))==0) %>% 
    column_to_rownames("protein_id")
  # add_part=add_part %>% 
  #   column_to_rownames(pav_df_raw$...2)
  add_part_num=sapply(add_part[2:157], function(x) as.numeric(x))
  pav_df_colsum=colSums(add_part_num)
  pav_df_colsum_sort=sort(pav_df_colsum)
  
  write.table(attributes(pav_df_colsum_sort),paste(result_path,sprintf("set_minus_sort_protein_id_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = F,col.names = F)
  write.table(pav_df_colsum_sort,paste(result_path,sprintf("set_minus_sort_protein_id_num_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = T,col.names = F)
  
  
  WriteXLS::WriteXLS(
    minus_part,
    paste(result_path,sprintf("set_minus_minus_%s.xlsx", start_point),sep = ""),
    col.names = T,
    row.names = T
  )
  WriteXLS::WriteXLS(
    add_part,
    paste(result_path,sprintf("set_minus_add_%s.xlsx", start_point),sep = ""),
    col.names = T,
    row.names = T
  )
  write.table(gene_is,paste(result_path,sprintf("set_minus_gene_id_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = F,col.names = F)
}
Cut('70-15')
"""

def R_set_minus_cut_orthofinder(orthofinder_tsv_file_name,orthofinder_unassianed_tsv_file_name,pan_id_file,set_minus_orthofinder_result):
    '''
    input 1: orthofinder_tsv_file_name
    input 2: orthofinder_unassianed_tsv_file_name
    input 3: pan_id_file
    output 1: set_minus_orthofinder_result
    '''
    robjects.globalenv["result_path"] = set_minus_orthofinder_result
    robjects.globalenv["orthofinder_unassianed_tsv_file_name"] = orthofinder_unassianed_tsv_file_name
    robjects.globalenv["orthofinder_tsv_file_name"] = orthofinder_tsv_file_name
    robjects.globalenv["pan_id_file"] = pan_id_file
    robjects.r(R_pav_orthofinder_code)