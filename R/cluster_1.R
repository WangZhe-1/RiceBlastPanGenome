require(readxl)
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(RColorBrewer)
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

write.table(fl_count_gene,file = pav_with_cluster_out_file_name,sep = "\t",append=FALSE,row.names = FALSE, col.names = TRUE)

fl_count_gene_write=fl_count_gene %>% 
  select(species,cluster)

write.table(
  fl_count_gene_write,file = cluster_out_file_name,sep="\t", row.names = FALSE, col.names = TRUE,append=FALSE,quote = F
)
fit_phy=as.phylo(fit)

ptree <- ggtree(
  fit_phy,
  color="black", 
  size=0.8, 
  linetype=1,  
  right=TRUE,
  layout="daylight",
  # branch.length = 'none'
) %<+% 
  fl_count_gene_write + 
  geom_tippoint(aes(color=as.character(cluster)),size=3)+
  scale_color_manual(values=alpha(color_clade$color, .8), 
                     name="Clade",
                     breaks=color_clade$clade
                     # labels=c("Control", "Treatment 1", "Treatment 2")
  )
  
print(ptree)
