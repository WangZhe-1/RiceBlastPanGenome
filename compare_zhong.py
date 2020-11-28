R="""
zhong_clade=read_xlsx("../../Pan_genome_data/R_result/zhong_cluster_clade.xlsx",skip = 2,n_max = 135)
zhong_clade=zhong_clade[,c(1,7)]
zhong_clade$Isolates=sub("\\(.+","",zhong_clade$Isolates)

zhe_clade=read.table("../../Pan_genome_data/cluster.txt",header = T)
all_clade=merge(zhong_clade,zhe_clade,by.x = 1,by.y = 1,all.x = T)
busco=read.table("../../Pan_genome_data/read_busco.txt")
all_clade_NA=all_clade[which(is.na(all_clade$cluster)),]
all_busco=merge(all_clade_NA,busco,by.x = 1,by.y = 1,all.x = T)
all_busco=all_busco[-c(5,14,37),]
all_busco_na=all_busco[which(is.na(all_busco$V2)),]
"""