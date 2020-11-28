# in_df_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/blast_general/blast_identity_value/MGG_00010T0.tsv"
in_df_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/blast_general/blast_identity_value/MGG_09379T0.tsv"
# in_df_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/blast_general/blast_identity_value/MGG_00016T0.tsv"

# single_vector_file_name="/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/ortholog_blast_mgg_key/all_row_gene/list/MGG_00010T0_singlecopy_list.txt"
single_vector_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/all_row_gene/list/MGG_09379T0_singlecopy_list.txt"
# multi_Vector_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/all_row_gene/list/MGG_00010T0_multicopy_list.txt"
multi_Vector_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/all_row_gene/list/MGG_09379T0_multicopy_list.txt"
lower_0.95_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/lower_0.95_file.txt"
best_sequence_id_list_file_name="../../Pan_genome_data_2/ortholog_blast_mgg_key/best_sequence_id_list.txt"
require(dplyr)
extract_max=function(R_blast_vlaue_df,in_vector){
  in_df=filter(R_blast_vlaue_df,gene_name %in% in_vector)
  return(in_df[which.max(in_df$blast_value),]$gene_name)
}


parse_blast_result_have_multi=function(
  single_vector_file_name,
  multi_Vector_file_name,
  in_df_file_name,
  lower_0.95_file_name,
  best_sequence_id_list_file_name
){
  R_blast_vlaue_df=read.table(in_df_file_name,header = T,stringsAsFactors = F)
  head_id=R_blast_vlaue_df$MGG_head[1]
  if(colnames(R_blast_vlaue_df)[1]=="NA."){
    print(c("fsad",head_id))
    return(0)
    }
  single_copy_gene_vector=read.table(single_vector_file_name,stringsAsFactors = F)
  multi_copy_gene_Vector=read.table(multi_Vector_file_name,stringsAsFactors = F)
  R_blast_vlaue_df_lower=R_blast_vlaue_df %>% 
    filter(blast_value<0.95)
  R_blast_vlaue_df=R_blast_vlaue_df %>% 
    filter(blast_value>0.95)
  if(nrow(R_blast_vlaue_df_lower)>0){
    write.table(R_blast_vlaue_df_lower,lower_0.95_file_name,quote = F,sep = '\t',row.names = F)
  }
  if(nrow(R_blast_vlaue_df)==0){
    print(c("sdaf",head_id))
    return(0)
    }
  single_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% single_copy_gene_vector$V1)
  multi_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% multi_copy_gene_Vector$V1)
  if(nrow(multi_blast_vlaue_df)==0){
    print(R_blast_vlaue_df$MGG_head[1])
    print(multi_copy_gene_Vector)
    return(parse_blast_result_no_multi(
      single_vector_file_name,
      in_df_file_name,
      lower_0.95_file_name,
      best_sequence_id_list_file_name
    ))
  }
  multi_blast_vlaue_df=mutate(multi_blast_vlaue_df,strain_id=sub("^(.+?)_.+","\\1",gene_name))
  dsaf=aggregate(gene_name~strain_id,data = multi_blast_vlaue_df,function(x) extract_max(R_blast_vlaue_df,x))
  multi_best_value_df=filter(multi_blast_vlaue_df,gene_name %in% dsaf$gene_name)
  multi_best_value_df=multi_best_value_df[,-4]
  best_value_df=rbind(single_blast_vlaue_df,multi_best_value_df)
  write_gene_id=c(best_value_df$gene_name,R_blast_vlaue_df$MGG_head[1])
  write.table(
    write_gene_id,
    best_sequence_id_list_file_name,
    quote = F,
    row.names = F,
    col.names = F
  )
  return(nrow(best_value_df))
}
parse_blast_result_no_multi=function(
  single_vector_file_name,
  in_df_file_name,
  lower_0.95_file_name,
  best_sequence_id_list_file_name
){
  R_blast_vlaue_df=read.table(in_df_file_name,header = T,stringsAsFactors = F)
  head_id=R_blast_vlaue_df$MGG_head[1]
  if(colnames(R_blast_vlaue_df)[1]=="NA."){
    print(c("fsad",head_id))
    return(0)
    }
  single_copy_gene_vector=read.table(single_vector_file_name,stringsAsFactors = F)
  R_blast_vlaue_df_lower=R_blast_vlaue_df %>% 
    filter(blast_value<0.95)
  R_blast_vlaue_df=R_blast_vlaue_df %>% 
    filter(blast_value>0.95)
  if(nrow(R_blast_vlaue_df_lower)>0){
    write.table(R_blast_vlaue_df_lower,lower_0.95_file_name,quote = F,sep = '\t',row.names = F)
  }
  if(nrow(R_blast_vlaue_df)==0){
    print(c("sdaf",head_id))
    return(0)
  }
  write_gene_id=c(R_blast_vlaue_df$gene_name,R_blast_vlaue_df$MGG_head[1])
  write.table(
    write_gene_id,
    best_sequence_id_list_file_name,
    quote = F,
    row.names = F,
    col.names = F
  )
  return(nrow(R_blast_vlaue_df))
}
parse_blast_result_have_multi(single_vector_file_name,multi_Vector_file_name,in_df_file_name,lower_0.95_file_name,best_sequence_id_list_file_name)
parse_blast_result_no_multi(single_vector_file_name,in_df_file_name,lower_0.95_file_name,best_sequence_id_list_file_name)
