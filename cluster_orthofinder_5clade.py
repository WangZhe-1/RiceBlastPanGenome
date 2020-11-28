import rpy2.robjects as robjects
from clade_color_5clade import create_color_clade
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
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


d=dist(fl,method = "binary")
fit=hclust(d,method = "ward.D2")

clus3=cutree(fit,color_num)

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
ggsave(plot = ptree,filename =paste(result_path,color_num, "_","cluster_binary.png",sep = ""))
ggsave(plot = ptree,filename =paste(result_path,color_num, "_","cluster_binary.svg",sep = ""))
'''
def cluster(
  pav_file_name,
  clade_set,
  color_set,
  color_num,
  pav_with_cluster_out_file_name,
  cluster_out_file_name,
  R_result
  ):
    '''
    input 1: pav_file_name
    input 2: clade_set, clade 1;clade 2;clade 3;clade 4,used as break,not lable, must be identitiy with cutree result
    input 3: color_set, https://colorbrewer2.org/
    input 4: color_num how much clade?
    output 1: pav_with_cluster_out_file_name
    output 2: strain vs cluster: cluster_out_file_name
    output 3: cluster plot, provide R_result
    '''
    robjects.globalenv["pav_file_name"] = pav_file_name
    robjects.globalenv["pav_with_cluster_out_file_name"]=pav_with_cluster_out_file_name
    robjects.globalenv["cluster_out_file_name"]=cluster_out_file_name
    robjects.globalenv["result_path"]=R_result
    robjects.globalenv["color_clade"]=create_color_clade(clade_set,color_set,color_num)
    robjects.globalenv["color_num"]=color_num
    
    robjects.r(R_cluster_code)
def write_cluster_result():
  '''
  write 2 file 
  1:pav_with_cluster_out_file_name
  2:strain id vs clade category ID
  '''
  robjects.r['write_cluster_result']()

def pca(pav_file_name,clade_id_file_name,out_png_file_name,out_eig_png_file_name):
  '''
  input 1: pav_matrix
  input 2: clade_id_file_name
  output 1: out_png_name
  output 2: out_eig_png_file_name
  '''
  R_pca_code='''
  pca_fun=function(in_file_name,clade_id_file_name,out_png_name,out_eig_png_file_name){
  require(data.table)
  require(readxl)
  require(factoextra)
  pav_input_raw=read_xlsx(in_file_name)
  pav_input_t=data.table::transpose(pav_input_raw[,3:158])
  colnames(pav_input_t)=pav_input_raw$protein_id
  rownames(pav_input_t)=colnames(pav_input_raw[,3:158])
  clade_num=read.table(clade_id_file_name,header = T)
  pca=prcomp(pav_input_t,scale. = F)
  eig=fviz_eig(pca)
  ggsave(out_eig_png_file_name,eig)

  pca_plot=fviz_pca_ind(
    pca,
    # geom.ind = c("point", "text"),
    geom.ind = c("point"),
    addEllipses = T, 
    legend.title = "Clade",
    habillage = clade_num$cluster
  )
  ggsave(out_png_name,pca_plot)
  }
  '''
  R_pca=SignatureTranslatedAnonymousPackage(R_pca_code,"R_pca")
  R_pca.pca_fun(pav_file_name,clade_id_file_name,out_png_file_name,out_eig_png_file_name)

def cluster_rice(
  pav_file_name,
  clade_set,
  color_set,
  color_num,
  out_path,
  pav_with_cluster_out_file_name,
  cluster_out_file_name,
  pav_rice_out_excel_name):
  '''
  input 1: pav_file_name
  input 2: clade_set, clade 1;clade 2;clade 3;clade 4,used as break,not lable, must be identitiy with cutree result
  input 3: color_set, https://colorbrewer2.org/
  input 4: color_num how much clade?
  output 1: out_path
  output 2: pav_with_cluster_out_file_name
  output 3: cluster_out_file_name
  output 4: pav_rice_out_excel_name
  '''
  R_cluster_rice_code='''
  cluster_rice=function(pav_file_name,color_clade,color_num,result_path,pav_with_cluster_out_file_name,cluster_out_file_name,pav_out_excel_name){
  require(tidyverse)
  require(ape)
  require(ggtree)
  require(data.table)
  require(WriteXLS)
  pav_raw=read.table(pav_file_name,header = T,check.names = F)
  pav_rice=pav_raw %>% 
    filter(cluster<3) %>% 
    column_to_rownames("species")
  pav_rice=pav_rice[,-24486]
  
  other_host_StrainID_df=pav_raw %>% 
    filter(cluster>2) %>% 
    select("species")
  geneID_sub=sub("_.+","",colnames(pav_rice))
  other_host_StrainID_need_remove_list=sapply(other_host_StrainID_df$species, function(x) grep(x,geneID_sub))
  other_host_StrainID_need_remove_unlist=unlist(other_host_StrainID_need_remove_list)
  pav_rice=pav_rice[,-other_host_StrainID_need_remove_unlist]
  pav=data.table::transpose(pav_rice)
  colnames(pav)=rownames(pav_rice)
  rownames(pav)=colnames(pav_rice)
  WriteXLS::WriteXLS(pav,pav_out_excel_name,row.names = T,col.names = T)

  d=dist(pav_rice,method = "binary")
  fit=hclust(d,method = "ward.D2")
  clus3=cutree(fit,color_num)
  
  pav_rice_clade=pav_rice %>% 
    rownames_to_column("species") %>% 
    mutate(cluster=clus3) 

  pav_rice_clade_write=pav_rice_clade %>% 
    select(species,cluster)
  write.table(pav_rice_clade,file = pav_with_cluster_out_file_name,sep = "\t",append=FALSE,row.names = FALSE, col.names = TRUE)
  write.table(
    pav_rice_clade_write,file = cluster_out_file_name,sep="\t", row.names = FALSE, col.names = TRUE,append=FALSE,quote = F
  )
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
    pav_rice_clade_write + 
    geom_tippoint(aes(color=as.character(cluster)),size=3)+
    scale_color_manual(values=alpha(color_clade$color, .8), 
                       name="Clade",
                       breaks=color_clade$clade
                       # labels=c("Control", "Treatment 1", "Treatment 2")
    )
  ggsave(plot = ptree,filename =paste(result_path,"rice_",color_num, "_","cluster_binary.png",sep = ""))
  ggsave(plot = ptree,filename =paste(result_path,"rice_",color_num, "_","cluster_binary.svg",sep = ""))
  }
  '''
  R_cluster_rice=SignatureTranslatedAnonymousPackage(R_cluster_rice_code,"R_cluster_rice")
  R_cluster_rice.cluster_rice(
    pav_file_name,
    create_color_clade(clade_set,color_set,color_num),
    color_num,
    out_path,
    pav_with_cluster_out_file_name,
    cluster_out_file_name,
    pav_rice_out_excel_name
    )