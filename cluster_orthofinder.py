'''
@Author: your name
@Date: 2020-08-12 12:00:00
@LastEditTime: 2020-08-12 16:07:41
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/cluster_orthofinder.py
'''
'''
@Author: your name
@Date: 2020-08-12 12:00:00
@LastEditTime: 2020-08-12 12:00:00
@LastEditors: your name
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/cluster_orthofinder.py
'''
import rpy2.robjects as robjects
from clade_color import create_color_clade
R_cluster_code='''
require(readxl)
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
# require(RColorBrewer)
fl_1=read_xlsx(pav_file_name)
fl=transpose(fl_1[,c(3:158)])
colnames(fl)=fl_1$protein_id
rownames(fl)=colnames(fl_1)[3:158]


d=dist(fl,method = "euclidean")
fit=hclust(d,method = "ward.D2")

clus3=cutree(fit,4)

fl_count_gene=fl %>% 
  rownames_to_column("species") %>% 
  mutate(cluster=clus3) 

fl_count_gene_write=fl_count_gene %>% 
  select(species,cluster)

write_cluster_result=function(){
  write.table(fl_count_gene,file = pav_with_cluster_out_file_name,sep = "\t",append=FALSE,row.names = FALSE, col.names = TRUE)
  write.table(
  fl_count_gene_write,file = cluster_out_file_name,sep="\t", row.names = FALSE, col.names = TRUE,append=FALSE,quote = F
)
}

fit_phy=as.phylo(fit)

ptree <- ggtree(
  fit_phy,
  color="black", 
  size=0.8, 
  linetype=1,  
  right=TRUE,
  layout="equal_angle",
  # branch.length = 'none'
) %<+% 
  fl_count_gene_write + 
  geom_tippoint(aes(color=as.character(cluster)),size=3)+
  scale_color_manual(values=alpha(color_clade$color, .8), 
                     name="Clade",
                     breaks=color_clade$clade
                     # labels=c("Control", "Treatment 1", "Treatment 2")
  )
ggsave(plot = ptree,filename =paste(result_path, "cluster.png",sep = ""))
ggsave(plot = ptree,filename =paste(result_path, "cluster.svg",sep = ""))
'''
def cluster(
  pav_file_name,
  clade_set,
  color_set,
  pav_with_cluster_out_file_name,
  cluster_out_file_name,
  R_result
  ):
    '''
    input 1: pav_file_name
    input 2: clade_set, clade 1;clade 2;clade 3;clade 4,used as break,not lable, must be identitiy with cutree result
    input 3: color_set, https://colorbrewer2.org/
    output 1: pav_with_cluster_out_file_name
    output 2: strain vs cluster: cluster_out_file_name
    output 3: cluster plot, provide R_result
    '''
    robjects.globalenv["pav_file_name"] = pav_file_name
    robjects.globalenv["pav_with_cluster_out_file_name"]=pav_with_cluster_out_file_name
    robjects.globalenv["cluster_out_file_name"]=cluster_out_file_name
    robjects.globalenv["result_path"]=R_result
    robjects.globalenv["color_clade"]=create_color_clade(clade_set,color_set)
    
    robjects.r(R_cluster_code)
def write_cluster_result():
  '''
  write 2 file 
  1:pav_with_cluster_out_file_name
  2:strain id vs clade category ID
  '''
  robjects.r['write_cluster_result']()