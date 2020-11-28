in_file=read.table('../../Pan_genome_data/cluster.txt',header = T,sep = "\t")
rice_host=in_file %>% 
  filter(!(cluster==3))

MGG_unpresent_Augustus_locate_in_pav_orthofinder=function(MGG_unpresent_Augustus_list_file_name,pav_orthofinder_file_name){
  MGG_unpresent_Augustus_list=read.table(MGG_unpresent_Augustus_list_file_name,stringsAsFactors = F)
  pav_orthofinde=read_xlsx(pav_orthofinder_file_name)
  pav_orthofinde=pav_orthofinde[,-1]
  
  gene_protein_mapping_table=read.table(gene_protein_mapping_table_file_name)
  
  pav_orthofinde_protein=merge(MGG_unpresent_Augustus_list,gene_protein_mapping_table,by.x = 1,by.y = 1)
  pav_MGG_unpresent_Augustus =pav %>% 
    filter(protein_id %in% pav_orthofinde_protein$V2)
  pav_MGG_unpresent_Augustus=pav_MGG_unpresent_Augustus[,c(158,1:157)]
  WriteXLS::WriteXLS(pav_MGG_unpresent_Augustus,"../../Pan_genome_data/70-15_MGG_Augustus/MGG_unpresent_Augustus_pav_orthofinder.xlsx")
}
library(readxl)
library(dplyr)

orthofinder_unassianed=read.table(orthofinder_unassianed_tsv_file_name,sep = "\t",header = T,check.names = F)
MGG_unpresent_Augustus_unassianed_list=intersect(pav_MGG_unpresent_Augustus$protein_id,orthofinder_unassianed$`70-15_protein`)
write.table(MGG_unpresent_Augustus_unassianed_list,MGG_unpresent_Augustus_unassianed_list_file_name,quote = F,row.names = F,col.names = F)

length(intersect(pav_MGG_unpresent_Augustus$protein_id,orthofinder_unassianed$`70-15_protein`))

  MGG_unpresent_Augustus_unassianed_list=read.table(MGG_unpresent_Augustus_unassianed_list_file_name)
  pav_orthofinder=read_xlsx(pav_orthofinder_file_name)
  pav_orthofinder_1574=pav_orthofinder %>% 
    filter(!(protein_id %in% MGG_unpresent_Augustus_unassianed_list$V1))
  pav_orthofinder_1574=pav_orthofinder_1574[,-1]
  pav_orthofinder_1574=pav_orthofinder_1574[,c(1,158,2:157)]
  WriteXLS::WriteXLS(pav_orthofinder_1574,pav_orthofinder_1574)
pav_df=read_xlsx(pav_orthofinder_file_name)
pav_orthofinder_1574_file_name=

library(tidyverse)
df <- data.frame(x = c("a","b","c","D"), y = c(1,2,3,4))
df_2 <- data.frame(x = c("a","b","c"), y = c(1,2,3))
[[match(df_2$x,df$x)]]
ggplot(data = df) +
  geom_rect(aes(x=x,
                y=y,
    xmin = match(df_2$x,x) - 0.3,
            xmax = match(df_2$x,x) + 0.3),
            ymin = 0, ymax = 2
    )

ggplot(data = df) +
  geom_rect(data = df, aes(x = x, y=y), xmin = as.numeric(df$x[[2]]) - 0.3,
            xmax = as.numeric(df$x[[3]]) + 0.3,
            ymin = 0, ymax = 2)

mydat_1 <- tibble(
  mymsmt = rep(c("bio", "bio", "den", "den"), 2),
  mylvl = c("NT", "till", "NT", "till", "no", "yes", "no", "yes"),
  mytrt = c(rep("tillage", 4), rep("herbicides", 4)),
  est = c(-60, -13, -65, -40, -16, -24, -49, -50),
  cilow = c(-85, -48, -78, -56, -61, -60, -68, -64),
  ciup = c(8, 45, -44, -18, 79, 42, -20, -31)) 
mydat=mydat_1 %>% 
  mutate(mylvln = rep(c(1, 2), 4))
ggplot(mydat, aes(est, mylvl)) + 
  geom_tile(aes(width = ciup-cilow, height=0.1),  fill="red", color="black") +
  geom_point() + 
  facet_grid(mytrt ~ mymsmt, scales = "free")
ggplot()+
geom_tile(
  data = blast_same_strain,
  aes(
    x=V1,
    y=1
  ),
  height=0.1,
  fill="red",
  color="black"
)
as.numeric(blast_same_strain$V1[[1]])
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
ggplot()+
geom_rect(
  data = blast_same_strain,
  aes(
    xmin=V1,
    xmax=V1,
  ),
  ymax=Inf,
  ymin=0,
  fill="red",
  color="black"
)
