# clade_cut_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/strain_clade_category_ID_4_blastn.txt"
clade_cut_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/strain_clade_category_ID_4_blastn.txt"
# mafft_result_dir_name="../../Pan_genome_data_2/mafft_result_test/"
# mafft_result_dir_name="/mnt/d/zhes_learning_space/software_in_ubuntu/PopGenome_data/fasta/"
mafft_result_dir_name="/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/mafft_result/"
rice_clade=c("1","2","3","4")
library(PopGenome)
clade_cut=function(in_df,clade_category){
  require(dplyr)
  return((dplyr::filter(in_df,cluster==clade_category))$species)
}
clade_cut_df=read.table(clade_cut_file_name,header = T,stringsAsFactors = F)
clade_cut_list=lapply(rice_clade,function(x) clade_cut(clade_cut_df,x))
mafft_result=readData(mafft_result_dir_name,populations = clade_cut_list)
mafft_result=neutrality.stats(mafft_result)
tajima=mafft_result@Tajima.D
mafft_result=F_ST.stats(mafft_result)
fst=mafft_result@nucleotide.F_ST
mafft_result=diversity.stats(mafft_result,pi = T)
pi_value=mafft_result@Pi
table(pi_value[,1]==0)
table(is.nan(tajima[,1])|is.na(tajima[,1]))
table(pi_value[,2]==0)
table(is.nan(tajima[,2])|is.na(tajima[,2]))
