pav_file_name_2="../../Pan_genome_data_2/set_minus_orthofinder_result/pav_orthofinder.xlsx"
clade_id_file_name_2_5="../../Pan_genome_data_2/set_minus_orthofinder_result/strain_clade_category_ID_5.txt"
clade_id_file_name_2_4="../../Pan_genome_data_2/set_minus_orthofinder_result/strain_clade_category_ID.txt"
pav_file_name_1="../../Pan_genome_data/set_minus_orthofinder_result/pav_orthofinder.xlsx"
clade_id_file_name_1_4="../../Pan_genome_data_2/Augustus_based_orthofinder/strain_clade_category_ID.txt"
pca_fun=function(in_file_name,out_png_name,clade_id_file_name_2){
  require(data.table)
  require(readxl)
  require(factoextra)
  pav_input_raw=read_xlsx(in_file_name)
  pav_input_t=data.table::transpose(pav_input_raw[,3:158])
  colnames(pav_input_t)=pav_input_raw$protein_id
  rownames(pav_input_t)=colnames(pav_input_raw[,3:158])
  clade_num=read.table(clade_id_file_name_2,header = T)
  pca=prcomp(pav_input_t,scale. = F)
  fviz_eig(pca)
  
  pca_plot=fviz_pca_ind(
    pca,
    geom.ind = c("point", "text"),
    addEllipses = T, 
    legend.title = "Clade",
    habillage = clade_num$cluster
  )
  ggsave(out_png_name,pca_plot)
}
# pca_fun(pav_file_name_2,paste("../../Pan_genome_data_2/","pav_2_5.png"),clade_id_file_name_2_5)
# pca_fun(pav_file_name_2,paste("../../Pan_genome_data_2/","pav_2_4.png"),clade_id_file_name_2_4)
pca_fun(pav_file_name_1,paste("../../Pan_genome_data_2/Augustus_based_orthofinder/","pav_1_4.png"),clade_id_file_name_1_4)
        