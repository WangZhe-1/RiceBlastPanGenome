R_blast_vlaue_df=read.table("/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/ortholog_blast/blast_general/blast_identity_value/OG0000006.tsv",header = T)
single_copy_gene_vector=read.table("/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/ortholog_blast/all_row_gene/list/OG0000006_singlecopy_list.txt")
multi_copy_gene_Vector=read.table("/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/ortholog_blast/all_row_gene/list/OG0000006_multicopy_list.txt")
lower_0.95_file_name="/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data_2/ortholog_blast/all_row_gene/list/lowr.txt"
parse_blast_result=function(single_copy_gene_vector,multi_copy_gene_Vector,R_blast_vlaue_df,single_err_file_name){
  single_blast_vlaue_df=dplyr::filter(R_blast_vlaue_df,gene_name %in% single_copy_gene_vector$V1)
  Outliers_stats_df=boxplot.stats(single_blast_vlaue_df$blast_value)
  Outliers_list=Outliers_stats_df$out
  upper=Outliers_stats_df$stats[5]
  lower=Outliers_stats_df$stats[1]
  if (length(Outliers_list)>0 &any(Outliers_list<0.95)){
    write.table(
      filter(R_blast_vlaue_df,blast_value %in% c(Outliers_list)),
      lower_0.95_file_name,
      quote = F,
      sep = "\t",
      row.names =F
      )
  }else{
    multi_blast_vlaue_df=R_blast_vlaue_df %>% 
      filter(gene_name %in% multi_copy_gene_Vector) %>% 
      mutate(ok = blast_value >= lower & blast_value <= upper)
    if(!(all(multi_blast_vlaue_df$ok))){
      multi_blast_vlaue_df=multi_blast_vlaue_df[multi_blast_vlaue_df$ok,]
    }
    multi_blast_vlaue_df=mutate(multi_blast_vlaue_df,strain_id=strsplit(gene_name,"_")) 
    dsaf=aggregate(blast_value~strain_id,data = multi_blast_vlaue_df,max)
    best_value_df=filter(multi_blast_vlaue_df,blast_vlaue %in% dsaf)
  }
  return(c(single_blast_vlaue_df$gene_name,best_value_df$gene_name))
}
main=function(single_copy_gene_vector,multi_copy_gene_Vector,R_blast_vlaue_df,single_err_file_name){
  require(stringr)
  require(dplyr)
  if (length(unique(R_blast_vlaue_df$MGG_head))==1){
    return(parse_blast_result(single_copy_gene_vector,multi_copy_gene_Vector,R_blast_vlaue_df,single_err_file_name))
  }else if(all(R_blast_vlaue_df[str_detect(R_blast_vlaue_df$gene_name,"MGG"),]$blast_value==1)){
    return(parse_blast_result(
      single_copy_gene_vector,
      multi_copy_gene_Vector,
      filter(R_blast_vlaue_df,MGG_head==MGG_head[1]),
      single_err_file_name
      ))
  }else{
    return("NA")
  }
}
main(single_copy_gene_vector,multi_copy_gene_Vector,R_blast_vlaue_df,single_err_file_name)