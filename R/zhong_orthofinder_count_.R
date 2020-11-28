zhong_clade=read_xlsx("../../Pan_genome_data/R_result/zhong_cluster_clade.xlsx",skip = 2,n_max = 135)
zhong_clade=zhong_clade[,c(1,7)]
zhong_clade$Isolates=sub("\\(.+","",zhong_clade$Isolates)

zhe_old_clade=read.table("../../Pan_genome_data/cluster.txt",header = T)
zhe_old_r_clade=read.table("../../Pan_genome_data/R_result/strain_clade_category_ID.txt",header = T)
zhe_old_oldr=merge(zhe_old_clade,zhe_old_r_clade,by.x = 1,by.y = 1,all = T)
all_clade=merge(zhong_clade,zhe_clade,by.x = 1,by.y = 1,all = T)
busco=read.table("../../Pan_genome_data/read_busco.txt")

all_clade_NA=all_clade[which(is.na(all_clade$cluster)),]
all_busco=merge(all_clade_NA,busco,by.x = 1,by.y = 1,all.x = T)
all_busco=all_busco[-c(5,14,37),]
all_busco_na=all_busco[which(is.na(all_busco$V2)),]
zhe_new_clade=read.table("../../Pan_genome_data_2/set_minus_orthofinder_result/strain_clade_category_ID.txt",header = T)
zhe=merge(zhe_new_clade,zhe_old_clade,by.x = 1,by.y = 1)
all_clade_zhe_old_new=merge(zhong_clade,zhe,by.x = 1,by.y = 1,all=T)
know_class=c("EI9604","EI9411","SV9623","SV9610","MG04","WHTQ","Br130","BdBar16-1","BdJes16-1","WBSS","4091-5-8")
all_clade_zhe_old_new[match(know_class,all_clade_zhe_old_new$Isolates),]
busco[match(know_class,busco$V1),]
zhe_5_clade=read.table("../../Pan_genome_data_2/set_minus_orthofinder_result/strain_clade_category_ID_5.txt",header = T)
zhe_5_old=merge(zhe_5_clade,zhe_old_clade,by.x = 1,by.y = 1,all = T)
colnames(zhe_5_old)=c("fgenesh_5","augustus")
zhong_zhe_5_old=merge(zhong_clade,zhe_5_old,by.x = 1,by.y = 1,all = T)
zhe_rice=read.table("../../Pan_genome_data_2/set_minus_orthofinder_result_2/rice_strain_clade_category_ID_3.txt",header = T)
zhong_zhe_5_old_rice=merge(zhong_zhe_5_old,zhe_rice,by.x = 1,by.y = 1,all = T)
zhe_rice_2=read.table("../../Pan_genome_data_2/set_minus_orthofinder_result_2/rice_strain_clade_category_ID_3_binary.txt",header = T)
zhong_zhe_5_old_rice_rice_2=merge(zhong_zhe_5_old_rice,zhe_rice_2,by.x = 1,by.y = 1,all = T)
write.table(zhong_zhe_5_old_rice_rice_2,"../../Pan_genome_data_2/set_minus_orthofinder_result_2/zhong_zhe_5_old_rice_rice_2.txt",sep = "\t",row.names = F)
fl=read_xlsx("../../Pan_genome_data_2/set_minus_orthofinder_result/pav_orthofinder.xlsx")
Sums(fl$FR13)
sum(fl$FR13)
orthofinder_tsv_file_name="../../Pan_genome_data_2/Results_Aug04/Orthogroups/Orthogroups.tsv"
orthofinder_unassianed_tsv_file_name="../../Pan_genome_data_2/Results_Aug04/Orthogroups/Orthogroups_UnassignedGenes.tsv"
tsv_assignedGenes_in_df=read.csv(orthofinder_tsv_file_name,sep = "\t",check.names = F,stringsAsFactors = F,na.strings = "")
tsv_UnassignedGenes_in_df=read.csv(orthofinder_unassianed_tsv_file_name,sep = "\t",check.names = F,stringsAsFactors = F,na.strings = "")
tsv_all_in_df=rbind(tsv_assignedGenes_in_df,tsv_UnassignedGenes_in_df)
str_detect(tsv_all_in_df,"MGG")

count_assignedGenes=read.csv("../../Pan_genome_data_2/Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv",sep = "\t",check.names = F,stringsAsFactors = F)
count_assignedGenes=count_assignedGenes[,-c(1,158)]
count_UnassignedGenes=tsv_UnassignedGenes_in_df[,-1]
na_UnassignedGenes=is.na(count_UnassignedGenes)
count_UnassignedGenes[na_UnassignedGenes]=0
count_UnassignedGenes[!na_UnassignedGenes]=1
count_UnassignedGenes=sapply(count_UnassignedGenes, function(x) as.numeric(x))
count_all=rbind(count_assignedGenes,count_UnassignedGenes)
count_sum=colSums(count_all)
count_sum[[ "FR13"]]


extract_strain_protein=function(strain)
  {
  assign(strain,fl$protein_id[as.logical(fl[[(sym(strain))]])],envir = .GlobalEnv)
}
clade_4=c("G17","G22","Arcadia","FR13")
clade_2=c("13FM-16-1","13FM-24-1","13FM-3-2","13FM-5-1","13FM-9-1","70-15","87-120")
FR13=extract_strain_protein("FR13")
lapply(clade_4,extract_strain_protein)
lapply(clade_2,extract_strain_protein)
intersect_fun=function(one,two){
  assign(paste(one,two,sep = "-"),intersect(get(one),get(two)),envir = .GlobalEnv)
}
by_sub=c("G17","G22","Arcadia","13FM-16-1","13FM-24-1","13FM-3-2","13FM-5-1","13FM-9-1","70-15","87-120")
FR13_subtract=function(by_subtract){
  intersect_fun("FR13",by_subtract)
}
lapply(by_sub, FR13_subtract)
intersect_fun("FR13","G17")

clade_blast=read.table("../../Pan_genome_data_2/set_minus_orthofinder_result_2/strain_clade_category_ID_4_blastn.txt",header = T)
zhong_blast=merge(zhong_clade,clade_blast,by.x = 1,by.y = 1)
write.table(zhong_blast,"../../Pan_genome_data_2/set_minus_orthofinder_result_2/zhong_blast.txt",row.names = F)
