strain_id=read.table("../../Pan_genome_data/R_result/strain_clade_category_ID.txt",header = T)
ncbi_blast=read.csv("../../Pan_genome_data/avr_2/HZVJZGPH01R-Alignment-Descriptions.csv",stringsAsFactors = F)

split_df=unlist(strsplit(ncbi_blast$Description,"\\s+"))

pib=(intersect(strain_id$species,split_df))
setdiff(strain_id$species,split_df)
pita=read.csv("../../Pan_genome_data/avr_2/HZVGZ62J01R-Alignment-Descriptions.csv",stringsAsFactors = F)
pita_df=unlist(strsplit(pita$Description,"\\s+"))

pita_inter=unlist(intersect(strain_id$species,pita_df))
pita_setd=setdiff(strain_id$species,pita_df)
rept=rep(3,times=10)
duplicated(rept)
pin_df=as.data.frame(table(pita_df))
