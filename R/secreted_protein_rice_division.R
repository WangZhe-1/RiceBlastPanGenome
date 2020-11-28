pav_0_1_tsv_file_name="../../Pan_genome_data_2/set_minus_orthofinder_result_2/rice_pav_with_cluster_3_binary.tsv"
classic_secreted_protein_file_name="../../Pan_genome_data_2/annotation/classic_secreted_protein.txt"
new_secreted_protein_file_name="../../Pan_genome_data_2/annotation/new_secreted_protein.txt"
result_path="../../Pan_genome_data_2/annotation/"
n=10
s=1
library(clusterProfiler)
library(tidyverse)
library(data.table)
library(ggrepel)
library(readxl)

fl=read.table(pav_0_1_tsv_file_name,header = T,check.names = F)
fl=fl[,-18380]
fl=fl %>% 
  column_to_rownames("species")
unverse_id_list <- colnames(fl)
# 用于读取哪些ID对应的是分泌蛋白

sp_list_classic <- read.csv(classic_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
sp_list_new_pre <- read.csv(new_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
sp_list_new <- sp_list_new_pre %>%
  mutate(V2 = V2 - 1)
sp_list_all <- rbind(sp_list_classic, sp_list_new)
term2gene_all <- as.data.frame(unverse_id_list, stringsAsFactors = F)
term2gene_all[, 2] <- 3

term2gene_all[match(sp_list_all$V1,term2gene_all$unverse_id_list,nomatch = 0), 2] <- 1

term2gene_all <- term2gene_all[, c(2, 1)]
term2name_all=data.frame(
  term=c(1,2),
  name=c("secretory_protein","secretory_protein")
)
result <- colSums(fl)
result_df <- as.data.frame(result)
n_species <- nrow(fl)


location <- result_df[match(known_secreted_proteins, rownames(result_df)), 1]
known_se <- data.frame(known_secreted_proteins, location)
write.table(known_se,paste(result_path, "known_secreted_proteins_location.tsv",sep = ""))

general_BgRatio <- NULL
enrich <- function(n,s) {
  result_list_BgRatio <- NULL
  result_list_p <- NULL
  result_name <- NULL
  result_detail=NULL
  Step <- ceiling((1 / n) * n_species)
  for (i in 1:(s*Step)) {
    Step_list <- result_df %>%
      rownames_to_column("gene") %>%
      filter(result == i) %>%
      select("gene")

    x <- enricher(
      Step_list$gene,
      TERM2GENE = term2gene_all,
      TERM2NAME = term2name_all,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      universe = unverse_id_list,
      minGSSize = 1,
      maxGSSize = 5000
    )

    if (is.null(x)) {
      result_list_p <- c(result_list_p, as.numeric(0))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(0))
    }
    else {
      result_list_p <- c(result_list_p, as.numeric(x@result$p.adjust))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(eval(parse(text = x@result$GeneRatio))))
    }
    result_name <- c(result_name, paste(i, i, sep = "-"))
    x@result$ID=c(paste(i, i, sep = "-"))
    result_detail=rbind(result_detail,x@result)
  }
  for (i in seq((s*Step)+1,nrow(fl)-s*Step,by=Step)) {
    from <- i
    End <- i + Step
    # if (End == n_species + 1) {
    #   break
    # }
    Step_list <- result_df %>%
      rownames_to_column("gene") %>%
      filter(result >= from & result < End) %>%
      select("gene")
    x <- enricher(
      Step_list$gene,
      TERM2GENE = term2gene_all,
      TERM2NAME = term2name_all,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      universe = unverse_id_list,
      minGSSize = 1,
      maxGSSize = 5000
    )
    
    if (is.null(x)) {
      result_list_p <- c(result_list_p, as.numeric(0))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(0))
    }
    else {
      result_list_p <- c(result_list_p, as.numeric(x@result$p.adjust))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(eval(parse(text = x@result$GeneRatio))))
    }
    result_name <- c(result_name, paste(from, End-1, sep = "-"))
    x@result$ID=c(paste(from, End-1, sep = "-"))
    result_detail=rbind(result_detail,x@result)
    
  }
  for (i in (nrow(fl)-s*Step+1):nrow(fl)){
    Step_list <- result_df %>%
      rownames_to_column("gene") %>%
      filter(result == i) %>%
      select("gene")

    x <- enricher(
      Step_list$gene,
      TERM2GENE = term2gene_all,
      TERM2NAME = term2name_all,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      universe = unverse_id_list,
      minGSSize = 1,
      maxGSSize = 5000
    )

    if (is.null(x)) {
      result_list_p <- c(result_list_p, as.numeric(0))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(0))
    }
    else {
      result_list_p <- c(result_list_p, as.numeric(x@result$p.adjust))
      result_list_BgRatio <- c(result_list_BgRatio, as.numeric(eval(parse(text = x@result$GeneRatio))))
    }
    result_name <- c(result_name, paste(i, i, sep = "-"))
    x@result$ID=c(paste(i, i, sep = "-"))
    result_detail=rbind(result_detail,x@result)
  }
  general_BgRatio <<- eval(parse(text = x@result$BgRatio))
  result_enrich <- data.frame(result_name, result_list_p, result_list_BgRatio)
  write.table(x = result_detail,file = paste(result_path, "secreted_proteins_rice_detail_",n,"_",s,".txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = F)
  write.table(result_enrich,paste(result_path, "secreted_proteins_rice_calculated_",n,"_",s,".txt",sep = ""),sep="\t",quote = F,row.names = F)
  return(result_enrich)
}

enrich_result <- enrich(n,s)

enrich_result$result_name <- factor(enrich_result$result_name, levels = enrich_result$result_name)

# general_BgRatio <- NULL
# enrich <- function(n) {
#   result_list_BgRatio <- NULL
#   result_list_p <- NULL
#   result_name <- NULL
#   result_detail=NULL
#   
#   general_BgRatio <<- eval(parse(text = x@result$BgRatio))
#   result_enrich <- data.frame(result_name, result_list_p, result_list_BgRatio)
#   write.table(x = result_detail,file = paste(result_path, "secreted_proteins_rice_detail_",n,"_qu",".txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = F)
#   return(result_enrich)
# }
# enrich(40)
# find_known_se=function()
# {
#   known_se
#   for (i in 1:nrow(known_se)){
#     if
#     enrich_known_se[i,]=enrich_result[grepl(paste("^",known_se$location[i],sep = ""),enrich_result$result_name),]
#   }
# }
# find_known_se()

p_value <- function() {
  p_result_sub <- subset(enrich_result, result_list_p < 0.05 & result_list_p > 0)
  p_value_line <- ggplot(
    data = enrich_result,
    aes(
      x = result_name,
      y = result_list_p
    ) 
  ) +
    geom_line(group = 1) +
    geom_point(
      data = p_result_sub,
      aes(
        x = result_name,
        y = result_list_p
      ),
      colour = "red"
    ) +
    geom_hline(yintercept = 0.05, lty = 2, lwd = 1, alpha = 0.8, colour = "red")+
    labs(
      x="Number of isolates",
      y="p.adjust"
    )+
    theme(
      panel.background=element_blank(),
      panel.border = element_blank(),
      panel.grid= element_blank(),
      axis.line = element_line(size=1, colour = "black")
    )
  # geom_label_repel(data = p_result_sub,
  #                  aes(
  #                    label=result_list_p,
  #                    x=result_name,
  #                    y=result_list_p
  #                  )
  #                  )
  ggsave(paste(result_path, "secreted_proteins_P_value_rice_",n,"_",s,".png",sep = ""), p_value_line, height = 5, width = 15)
  ggsave(paste(result_path, "secreted_proteins_P_value_rice_",n,"_",s,".svg",sep = ""), p_value_line, height = 5, width = 15)
}
p_value()

gene_radio <- function() {
  BgRatio_result_sub <- subset(enrich_result, result_list_BgRatio > general_BgRatio)
  BgRatio_line <- ggplot(
    data = enrich_result,
    aes(
      x = result_name,
      y = result_list_BgRatio
    )
  ) +
    geom_line(group = 1) +
    geom_point(
      data = BgRatio_result_sub,
      aes(
        x = result_name,
        y = result_list_BgRatio
      ),
      colour = "blue"
    ) +
    geom_hline(yintercept = general_BgRatio, lty = 2, lwd = 1, alpha = 0.8, colour = "blue")+
    theme(
      panel.background=element_blank(),
      panel.border = element_blank(),
      panel.grid= element_blank(),
      axis.line = element_line(size=1, colour = "black")
    )+
    labs(
      x="Number of isolates",
      y="GeneRatio"
    )
  # geom_label_repel(data = p_result_sub,
  #                  aes(
  #                    label=result_list_p,
  #                    x=result_name,
  #                    y=result_list_p
  #                  )
  #                  )
  # print(BgRatio_line)
  ggsave(paste(result_path, "secreted_proteins_BgRatio_line_rice_",n,"_",s,".png",sep = ""), BgRatio_line, height = 5, width = 15)
  ggsave(paste(result_path, "secreted_proteins_BgRatio_line_rice_",n,"_",s,".svg",sep = ""), BgRatio_line, height = 5, width = 15)
}
gene_radio()
