`clade_1&clade_3_distinct`=read.table("../../Pan_genome_data_2/combination_matrix_distinct_core/rice/clade_1&clade_3_distinct.txt")
clade_combination_set=c("clade_1&clade_3")
what_kind="distinct"
in_df=read.table("../../Pan_genome_data_2/clade_different_enrichment/kegg_annotation.txt",header = T)
classic_secreted_protein_file_name="../../Pan_genome_data_2/annotation/classic_secreted_protein.txt"
new_secreted_protein_file_name="../../Pan_genome_data_2/annotation/new_secreted_protein.txt"
universal_gene=read.table("../../Pan_genome_data_2/pan/rice_pan_gene.txt",stringsAsFactors = F)

library(clusterProfiler)
library(tidyverse)
library(data.table)
library(ggrepel)
enrich_result_list=list()
i=1
sp_list_classic <- read.csv(classic_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
sp_list_new_pre <- read.csv(new_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
sp_list_new <- sp_list_new_pre %>%
  mutate(V2 = V2 - 1)
sp_list_all <- rbind(sp_list_classic, sp_list_new)
sp_list_all$V1=sub("T0$","",sp_list_all$V1)
term2gene_all <- as.data.frame(universal_gene, stringsAsFactors = F)
term2gene_all[, 2] <- 3

term2gene_all[match(sp_list_all$V1,term2gene_all$V1,nomatch = 0), 2] <- 1

term2gene_all <- term2gene_all[, c(2, 1)]
term2name_all=data.frame(
  term=c(1,2),
  name=c("secretory_protein","secretory_protein")
)

fun_enrich=function(clade_comb){
  require(clusterProfiler)
  # print(head(gene_taxonomy))
  # print(head(annotation_case))
  # print(annotation_case)
  # print(str(gene_taxonomy))
  # print(str(annotation_case))
  gene_taxonomy=get(paste(clade_comb,"_",which_mode,sep=""))
  enrich_result=enricher(gene_taxonomy$V1,
                         TERM2GENE = term2gene_all,
                         TERM2NAME = term2name_all,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1, 
                         universe = universal_gene$V1,
                         minGSSize = 1,
                         maxGSSize=5000
  )
  # return(counter(enrich_result@result))
  enrich_result_list[[i]] <<- enrich_result@result
  names(enrich_result_list)[i]<<-clade_comb
  i <<- i + 1
  # write.table(
  # enrich_result@result,file = paste(result_path,sprintf('%s_%s_enrich_result.txt',as.character(substitute(gene_taxonomy)),which_mode),sep = ""), row.names = F, col.names = T,append=FALSE,quote = F,sep="\t"
  # )
  assign(sprintf('%s_sum_enrich',clade_comb),sum(enrich_result@result$Count),envir = .GlobalEnv)
  # print("name")
  # print(deparse(substitute(gene_taxonomy)))
  # print(str(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy)))))
  # print(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy))))
  # list.append(enrich_result_list,substitute(gene_taxonomy)=enrich_result@result)
  return(enrich_result)
}

run=function(){
  sum_enrich=0
  for (clade_comb in clade_combination_set){
    # assign(sprintf('no_present_%s_%s_enrich',kind,level_name),fun_enrich(no_present_id,in_df[,c(2,1)],in_df[,c(2,3)],which_mode),envir = .GlobalEnv)
    # assign(sprintf('no_present_%s_%s_count',kind,level_name),counter(get(sprintf('no_present_%s_%s_enrich',kind,level_name))@result,"no_present"),envir = .GlobalEnv)
    assign(sprintf('%s_enrich',deparse(substitute(clade_comb))),fun_enrich(clade_comb),envir = .GlobalEnv)
    sum_enrich=sum_enrich+get(sprintf('%s_sum_enrich',clade_comb))
  }
  # str(enrich_result_list)
  print(names(enrich_result_list))
  WriteXLS::WriteXLS(
    enrich_result_list,
    paste(result_path,sprintf('clade_different_%s_secreted_protein_enrich_result.xlsx',which_mode),sep = ""),
    # SheetNames=names(enrich_result_list),
    col.names = T,
    row.names = F
  )
  print(sum_enrich)
}
run()