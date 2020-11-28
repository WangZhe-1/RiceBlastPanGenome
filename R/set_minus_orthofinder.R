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

na_tsv=is.na(pav_df)
pav_df[na_tsv]=0
pav_df[!na_tsv]=1


pav_df=data.frame(pav_df,check.names = F,stringsAsFactors = F)
pav_df=pav_df %>% 
  rownames_to_column("Orthogroup")
# 删掉MGG_06302 MGG_06131 MGG_16639，他们三个有两个或以上的转录本，保留在各个菌株分布范围最广的那个

# 把基因id和og号整合好放进来
pav_df=merge(pav_df,pan_id,by.x = 1,by.y = 1,all.y = T)
names(pav_df)[names(pav_df) == 'V2'] <- 'protein_id'

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
# 读unassigned中70-15的部分，1860个
pav_unassigned_70_15=openxlsx::read.xlsx("..\\result\\pav_1861.xlsx",rowNames=T,colNames = T)
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) == "70-15_supercontigs"] <- '70-15'
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) == "X2" ] <-  "protein_id" 
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) == "Ina168_contig"] <- "ina168"
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) == "FJ81278_PB"] <- "FJ81278" 
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) == "FR13_PB"  ] <- "FR13" 
names(pav_unassigned_70_15)[names(pav_unassigned_70_15) ==   "Guy11_PB"  ] <- "Guy11" 
protein_id=pav_unassigned_70_15$protein_id
pav_unassigned_70_15=pav_unassigned_70_15[,-1]
pav_unassigned_70_15[pav_unassigned_70_15<=1]=0
pav_unassigned_70_15[pav_unassigned_70_15>1]=1
pav_unassigned_70_15=cbind(pav_unassigned_70_15,protein_id)

# 删掉从orthofinder结果中读取的unassignedGene中的70-15基因未处理部分
TE_need_remove=read.table("../70-15中有1891个unassigned_protein_远多于其他菌株_查看怎么回事/TE_present_need_remove.txt",stringsAsFactors = F)
pav_df=pav_df[!(sub("T.+$","",pav_df$protein_id) %in% pav_unassigned_70_15$protein_id),]
pav_df=pav_df[!(sub("T.+$","",pav_df$protein_id) %in% TE_need_remove$V1),]
# 把基因ID和og号整合放到处理好的unassignedGene的70-15部分
pan_id_unassigned=pan_id
pan_id_unassigned$V2=sub("T.+$","",pan_id$V2)
pav_unassigned_70_15_og=merge(pav_unassigned_70_15,pan_id_unassigned,by.x = 119,by.y = 2,all.x = T)
names(pav_unassigned_70_15_og)[names(pav_unassigned_70_15_og) == "V1"  ] <- "Orthogroup"
# 把处理好的unassignedGene的70-15部分和删除70-15基因未处理部分的pav合并
pav_df=rbind(pav_df,pav_unassigned_70_15_og)
# 写结果到excel
rownames(pav_df)=sprintf("pan_%05d", 1:nrow(pav_df))
write.xlsx(pav_df,"..\\result\\pav_all_final_remove_TE.xlsx",rowNames=T)

# 验证
count_assignedGenes=read.csv("..\\Results_Feb11\\Orthogroups\\Orthogroups.GeneCount.tsv",sep = "\t",check.names = F,stringsAsFactors = F)
count_assignedGenes=count_assignedGenes[,-c(1,120)]
count_UnassignedGenes=tsv_UnassignedGenes_in_df[,-1]
na_UnassignedGenes=is.na(count_UnassignedGenes)
count_UnassignedGenes[na_UnassignedGenes]=0
count_UnassignedGenes[!na_UnassignedGenes]=1
count_UnassignedGenes=sapply(count_UnassignedGenes, function(x) as.numeric(x))
count_all=rbind(count_assignedGenes,count_UnassignedGenes)
count_sum=colSums(count_all)
count_sum[[ "magnaporthe_oryzae_70-15_8_proteins"]]
