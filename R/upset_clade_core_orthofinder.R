pav_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/pav_orthofinder.xlsx"
category_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/strain_clade_category_ID_4_blastn.txt"
require(readxl)
require(ComplexHeatmap)
pav_raw=readxl::read_xlsx(pav_file_name)
pav_df=pav_raw[,-c(1:2,159)]

pav_logical=as.data.frame(apply(pav_df,2,as.logical))
strain_clade=read.table(category_file_name,header = T,stringsAsFactors = F)

union_clade_element=function(clade_set){
  result=pav_logical[[clade_set[1]]]
  for (clade in clade_set) {
    result=result|pav_logical[[clade]]
  }
  return(result)
}
clade_binary_df_raw=aggregate(species ~ cluster, data = strain_clade, union_clade_element)

clade_binary_df=data.frame(
  "1"=clade_binary_df_raw$species[1,],
  "2"=clade_binary_df_raw$species[2,],
  "3"=clade_binary_df_raw$species[3,],
  "4"=clade_binary_df_raw$species[4,]
)
colnames(clade_binary_df)=seq(1,4)
comb_matrix=make_comb_mat(clade_binary_df,mode = "distinct")
UpSet(comb_matrix)
comb_size(comb_matrix)
