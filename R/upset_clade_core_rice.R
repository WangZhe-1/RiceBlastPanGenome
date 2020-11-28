pav_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/pav_blastn_gene.xlsx"
category_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/strain_clade_category_ID_4_blastn.txt"
out_dir="../../Pan_genome_data_2/set_minus_orthofinder_result_2/"
zhong_zhe_merge=data.frame(
  zhong=c("clade_1","clade_2","clade_3","other"),
  zhe=c("1","4","2","3")
)
zhong_zhe_merge=data.frame(
  zhong=c("clade_1","clade_2","clade_3"),
  zhe=c("1","4","2")
)
which_mode="distinct"

# read_all_pav=function(pav_file_name){
  require(readxl)
  pav_raw=readxl::read_xlsx(pav_file_name)
  pav_df=pav_raw[,-c(1:2)]
  na_tsv=pav_df>=2
  pav_df[na_tsv]=1
  pav_df[!na_tsv]=0
  pav_logical=as.data.frame(apply(pav_df,2,as.logical))
  return(pav_logical)
}

  require(readxl)
  pav_raw=readxl::read_xlsx(pav_file_name)
  strain_clade<<-read.table(category_file_name,header = T,stringsAsFactors = F)
  rice_strainID=dplyr::filter(strain_clade,cluster %in% c("1","4","2"))
  pav_rice=pav_raw[rice_strainID$species]
  
  geneID_sub=sub("_.+","",pav_raw$...2)
  other_host_StrainID_df=dplyr::filter(strain_clade,cluster ==3)
  other_host_StrainID_need_remove_list=sapply(other_host_StrainID_df$species, function(x) grep(x,geneID_sub))
  other_host_StrainID_need_remove_unlist=unlist(other_host_StrainID_need_remove_list)
  pav_rice=pav_rice[-other_host_StrainID_need_remove_unlist,]
  general_protein=pav_raw$...2
  general_protein<<-general_protein[-other_host_StrainID_need_remove_unlist]
  write.table(general_protein,'../../Pan_genome_data_2/pan/pan_gene_no_shorter_20.fasta/rice_pan_gene.txt'
,row.names = F,col.names = F,quote = F)
  write.table(pav_rice,"",row.names = F,col.names = T,quote = F,sep = "/t")
  
  na_tsv=pav_rice>=2
  pav_rice[na_tsv]=1
  pav_rice[!na_tsv]=0
  pav_logical=as.data.frame(apply(pav_rice,2,as.logical))
  return(pav_logical)


extract_comb_gene=function(which_comb,out_dir){
  result=general_protein[extract_comb(comb_matrix,comb_label[which_comb])]
  write.table(result,paste(out_dir,which_comb,"_",which_mode,".txt"),row.names = F,col.names = F,quote = F)
}

  require(ComplexHeatmap)
  strain_clade=read.table(category_file_name,header = T,stringsAsFactors = F)
  strain_clade=merge(strain_clade,zhong_zhe_merge,by.x = 2,by.y = 2)
  
  union_clade_element=function(clade_set){
    result=pav_logical[[clade_set[1]]]
    for (clade in clade_set) {
      result=result|pav_logical[[clade]]
    }
    return(result)
  }
  clade_binary_df_raw=aggregate(species ~ zhong, data = strain_clade, union_clade_element)
  
  datalist = list()
  for (i in 1:nrow(zhong_zhe_merge)) {
    datalist[[i]]=clade_binary_df_raw$species[i,]
    names(datalist)[i]=as.character(clade_binary_df_raw$zhong[i])
  }
  clade_binary_df = do.call(cbind, datalist)
  comb_matrix=make_comb_mat(clade_binary_df,mode = which_mode)
  comb_label=structure(comb_name(comb_matrix), names = comb_name(comb_matrix,readable = T))
  general_protein=pav_raw$...2
  clade_combination=c("clade_1&clade_2","clade_2&clade_3","clade_1&clade_3")
  sapply(clade_combination,function(x) extract_comb_gene(x,out_dir))
}
