library(readxl)
library(tidyverse)
library(data.table)
fl_1=read_xlsx(pav_file_name)
fl_2=fl_1[,-c(1:2)]
fl=transpose(fl_2)
colnames(fl)=fl_1$...1
rownames(fl)=colnames(fl_2)

d=dist(fl,method = "euclidean")
fit=hclust(d,method = "ward.D2")

clus3=cutree(fit,4)

fl_count_gene=fl %>% 
  rownames_to_column("species") %>% 
  mutate(cluster=clus3) 

write.table(fl_count_gene,file = pav_with_cluster_out_file_name,sep = "\t",append=FALSE)
fl_count_gene_write=fl_count_gene %>% 
  select(species,cluster)

write.table(
  fl_count_gene_write,file = cluster_out_file_name,sep="\t", row.names = FALSE, col.names = TRUE,append=FALSE
)

#X>=2意味着基因存在，这个函数自上而下扫描149个物种，如果基因在149个（全部物种）物种中都存在，
#149个true，那么他就是core gene
count_2=function(x)
{
  table_count=table(x>=2)
  if('TRUE' %in% unlist(dimnames(table_count)))
  {
    return(as.numeric(table_count[['TRUE']]))
  }
  else
  {
    return(0)
  }
}
result=sapply(fl_count_gene,count_2)
result_df=as.data.frame(result)

how_many_149=as.data.frame(table(result),stringsAsFactors=FALSE)
num_0.95_ceiling=ceiling(0.95*num_near_rice_specise)
num_0.95_flooring=floor(0.95*num_near_rice_specise)
sum_0.95_s_flooring=how_many_149 %>% 
  filter(as.numeric(result)>=num_0.95_flooring)
sum_0.95_s_ceiling=how_many_149 %>% 
  filter(as.numeric(result)>=num_0.95_ceiling)
sum_0.95_flooring=sum(sum_0.95_s_flooring$Freq)
sum_0.95_ceiling=sum(sum_0.95_s_ceiling$Freq)
num_0.05_ceiling=ceiling(0.05*num_near_rice_specise)
sum_0.05_s_ceiling=how_many_149 %>% 
  filter(as.numeric(result)<=num_0.05_ceiling)
sum_0.05_ceiling=sum(sum_0.05_s_ceiling$Freq)
middle=as.numeric(12388) - as.numeric(sum_0.05_ceiling) - as.numeric(sum_0.95_ceiling)
cat(sprintf('above %d is %d\n',num_0.95_flooring,sum_0.95_flooring),file = ".\\结果\\how_many_149.txt",append = FALSE)
cat(sprintf('above %d is %d\n',num_0.95_ceiling,sum_0.95_ceiling),file = ".\\结果\\how_many_149.txt",append = TRUE)
cat(sprintf('blew %d is %d\n',num_0.05_ceiling,sum_0.05_ceiling),file = ".\\结果\\how_many_149.txt",append = TRUE)
cat(sprintf('middle is %d\n',middle),file = ".\\结果\\how_many_149.txt",append = TRUE)
write.table(
  how_many_149,file = ".\\结果\\how_many_149.txt",sep="\t", row.names = FALSE, col.names = TRUE,append=TRUE
)
core_gene=result_df %>% 
  rownames_to_column("gene") %>% 
  filter(result>=num_0.95_ceiling) %>% 
  select(gene)
write.table(
  core_gene,file = "core_gene.txt", row.names = F, col.names = F,append=FALSE
)


fl_clus3=fl %>% 
  rownames_to_column('newcio') %>%
  mutate(cluster=clus3) %>% 
  column_to_rownames('newcio')
fl_shui=subset(fl_clus3,cluster<4)

d_fei=dist(fl_shui,method = "euclidean")
fit_fei=hclust(d_fei,method = "ward.D2")
clus3_fei=cutree(fit_fei,3)
colors_3 = c("red", "blue", "green")
plot(as.phylo(fit_fei), 
     type = "unrooted",
     label.offset = 0.1,
     cex = 0.6,
     show.tip.label = TRUE,
     no.margin = TRUE,
     tip.color=colors_3[clus3_fei]
     #edge.color=colors[clus3]
)
plot(as.phylo(fit_fei), cex = 1, label.offset = 0.4,tip.color=colors_3[clus3_fei])
