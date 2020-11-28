require(data.table)
require(tidyverse)
require(WriteXLS)
# View(pav_t[1:10,1:10])
cluster_clade=read.table(cluster_clade_file_name,header = T)
pav=readxl::read_xlsx(pav_orthofinder_file_name)
pav_t=data.table::transpose(pav[,3:158])
colnames(pav_t)=pav$protein_id
rownames(pav_t)=colnames(pav[,3:158])
pav_t=pav_t %>% 
  rownames_to_column("strain_id")
pav_t_clade=merge(pav_t,cluster_clade,by.x = 1,by.y = 1)
setdiff(cluster_clade$species,pav_t$strain_id)
avr_list=c("MGG_18041T0", "MGG_13283T0", "MGG_03685T0", "MGG_04301T0", "MGG_12655T0")
clade_avr=function(srecre){
  clade_avr=numeric()
  for (i in c(1,2,4)){
    clade=pav_t_clade %>% 
      filter(cluster==i) 
    strain_sum=clade%>% 
      select(!!sym(srecre)) %>% 
      sum
    clade_avr=c(clade_avr,strain_sum/nrow(clade))
  }
  return(clade_avr)
}
ds=sapply(avr_list,clade_avr)
ds_df=as.data.frame(ds)
WriteXLS::WriteXLS(ds_df,"srepro.xlsx")
clade=pav_t_clade %>% 
  filter(cluster==4) %>% 
  select(MGG_18041T0) %>% 
  sum
