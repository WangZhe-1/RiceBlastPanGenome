#!/usr/bin/env python3
# -*- coding: utf-8 -*-
     
from pathlib import Path
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
import re
rgx_TMH_len=re.compile('Length:\s+(\d+)')
rgx_TMH_num=re.compile("^#\s(\w+)\sNumber of predicted TMHs:\s+(\d+)")
rgx_TMH_site=re.compile("TMhelix\s+(\d+)\s+(\d+)")
rgx_ProtComp_Secreted=re.compile("Integral.+Secreted")
rgx_ProtComp_id_pre=re.compile("Seq name:")
rgx_ProtComp_id_mgg=re.compile("MGG_.+T0")
rgx_ProtComp_id_other=re.compile(":\s+(.+),")
rgx_category_id=re.compile("out_(.+)_")
def get_Protcomp_ID_name(line):
    result_mgg=rgx_ProtComp_id_mgg.search(line)
    if result_mgg is not None:
        return result_mgg.group(0)
    else:
        return rgx_ProtComp_id_other.search(line).group(1)
def secreted_protein(
    Protcomp_file,
    tmhmm_file,
    signalp_file,
    new_secreted_protein_file_name,
    classic_secreted_protein_file_name
    ):
    '''
    input 1: Protcomp_file
    input 2: tmhmm_file
    input 3: signalp_file
    output 1: new_secreted_protein_file_name
    output 2: classic_secreted_protein_file_name
    '''
    

    dic_TMH_num={}
    dic_TMH_stat={}
    dic_TMH_end={}
    dic_ProtComp={}
    dic_signalP_end={}


    with open(tmhmm_file) as tmhmm_fl:
        for line in tmhmm_fl:
            TMH_len_result=rgx_TMH_len.search(line)
            if TMH_len_result is not None:
                TMH_len_result=re.split("\s",line)
                assert len(TMH_len_result)==5
                protein_len_T=int(TMH_len_result[3])
                if protein_len_T > 400:
                    continue
                line=next(tmhmm_fl)
                TMH_num_result=re.split("\s",line)
                assert len(TMH_num_result)==9
                protein_id_T=TMH_num_result[1]
                TMH_num=TMH_num_result[7]
                if int(TMH_num)==0:
                    dic_TMH_num.setdefault(str(protein_id_T),int(TMH_num))   
                elif int(TMH_num) ==1:
                    dic_TMH_num.setdefault(str(protein_id_T),int(TMH_num))
                    result_TMH_site=rgx_TMH_site.search(line)
                    while result_TMH_site is None:
                        line=next(tmhmm_fl)
                        result_TMH_site=rgx_TMH_site.search(line)
                    if result_TMH_site is not None:
                        TMH_start=result_TMH_site.group(1)
                        TMH_end=result_TMH_site.group(2)
                        dic_TMH_stat.setdefault(str(protein_id_T),int(TMH_start))
                        dic_TMH_end.setdefault(str(protein_id_T),int(TMH_end))

    with open(signalp_file) as signalp_fl:
        next(signalp_fl)
        for line_S in signalp_fl:
            protein_id_S=line_S.split('\t')[0]
            signalp_end=line_S.split('\t')[4]
            dic_signalP_end.setdefault(str(protein_id_S),int(signalp_end))

    with open(new_secreted_protein_file_name,"w+") as fl_out_n:
        with open(classic_secreted_protein_file_name,"w+") as fl_out_c:
            with open(Protcomp_file) as Protcomp_fl:
                for line_ProtComp in Protcomp_fl:
                    result_ProtComp_Secreted=rgx_ProtComp_Secreted.search(line_ProtComp)
                    result_ProtComp_id_pre=rgx_ProtComp_id_pre.search(line_ProtComp)
                    if result_ProtComp_id_pre is not None:
                        ProtComp_id_pre=line_ProtComp
                    if result_ProtComp_Secreted is not None:
                        ProtComp_id=get_Protcomp_ID_name(ProtComp_id_pre)
                        if ProtComp_id not in dic_TMH_num:
                            continue
                        if dic_TMH_num[ProtComp_id]==0:
                            if ProtComp_id in dic_signalP_end:
                                fl_out_c.write("{}\t{}\n".format(ProtComp_id,1))
                            else:
                                fl_out_n.write("{}\t{}\n".format(ProtComp_id,2))
                        elif dic_TMH_num[ProtComp_id]==1:
                            if ProtComp_id in dic_signalP_end:
                                if(int(dic_signalP_end[ProtComp_id])>int(dic_TMH_stat[ProtComp_id])):
                                    fl_out_c.write("{}\t{}\n".format(ProtComp_id,1))

def drawer_secreted_protein(
    pav_0_1_tsv_file_name,
    new_secreted_protein_file_name,
    classic_secreted_protein_file_name,
    known_secreted_proteins,
    R_result
    ):
    '''
    input 1: pav_0_1_tsv_file_name
    input 2: new_secreted_protein_file_name
    input 3: classic_secreted_protein_file_name
    input 5: known_secreted_proteins
    output are all in R_result Directory
    output 1: secreted_proteins_detail_file_name
    output 2: secreted_proteins_P_value_picture_name
    output 3: secreted_proteins_BgRatio_line_picture_name
    output 4: known_secreted_proteins_location
    '''
    R_code_draw_secreted_protein='''
    library(clusterProfiler)
    library(tidyverse)
    library(data.table)
    library(ggrepel)
    library(readxl)

    fl=readxl::read_xlsx(pav_0_1_tsv_file_name)
    fl_tr <- transpose(fl[,3:158])
    colnames(fl_tr) <- fl$protein_id
    rownames(fl_tr) <- colnames(fl[,3:158])
    unverse_id_list <- fl$protein_id
    # 用于读取哪些ID对应的是分泌蛋白

    sp_list_classic <- read.csv(classic_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
    sp_list_new_pre <- read.csv(new_secreted_protein_file_name, header = F, stringsAsFactors = F, sep = "\t")[, c(2, 1)]
    sp_list_new <- sp_list_new_pre %>%
    mutate(V2 = V2 - 1)
    sp_list_all <- rbind(sp_list_classic, sp_list_new)
    term2gene_all <- as.data.frame(unverse_id_list, stringsAsFactors = F)
    term2gene_all[, 2] <- 3

    term2gene_all[match(sp_list_all$V1,term2gene_all$unverse_id_list), 2] <- 1

    term2gene_all <- term2gene_all[, c(2, 1)]
    term2name_all=data.frame(
    term=c(1,2),
    name=c("secretory_protein","secretory_protein")
    )
    result <- colSums(fl_tr)
    result_df <- as.data.frame(result)
    n_species <- nrow(fl_tr)


    location <- result_df[match(known_secreted_proteins, rownames(result_df)), 1]
    known_se <- data.frame(known_secreted_proteins, location)
    write.table(known_se,paste(result_path, "known_secreted_proteins_location.tsv",sep = ""))

    general_BgRatio <- NULL
    enrich <- function(n) {
    result_list_BgRatio <- NULL
    result_list_p <- NULL
    result_name <- NULL
    result_detail=NULL
    Step <- ceiling((1 / n) * n_species)
    for (i in 0:n_species) {
        from <- i
        End <- i + Step
        if (End == n_species + 1) {
        break
        }
        if (End == n) {
        Step_list <- result_df %>%
            rownames_to_column("gene") %>%
            filter(result >= from & result <= End) %>%
            select("gene")
        }
        else {
        Step_list <- result_df %>%
            rownames_to_column("gene") %>%
            filter(result >= from & result < End) %>%
            select("gene")
        }
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
        result_name <- c(result_name, paste(from, End, sep = "-"))
        x@result$ID=c(paste(from, End, sep = "-"))
        if(from==0){
        result_detail=x@result
        }
        else{
        result_detail=rbind(result_detail,x@result)
        }
    }
    general_BgRatio <<- eval(parse(text = x@result$BgRatio))
    result_enrich <- data.frame(result_name, result_list_p, result_list_BgRatio)
    write.table(x = result_detail,file = paste(result_path, "secreted_proteins_detail.txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = F)
    return(result_enrich)
    }
    enrich_result <- enrich(10)
    write.table(enrich_result,paste(result_path, "secreted_proteins_calculated.txt",sep = ""),sep="\t",quote = F,row.names = F)
    enrich_result$result_name <- factor(enrich_result$result_name, levels = enrich_result$result_name)


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
        x = as.numeric(rownames(enrich_result)),
        y = result_list_p
        ) 
    ) +
        geom_line(group = 1) +
        geom_point(
        data = p_result_sub,
        aes(
            x = as.numeric(rownames(p_result_sub)),
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
    ggsave(paste(result_path, "secreted_proteins_P_value.png",sep = ""), p_value_line, height = 5, width = 15)
    ggsave(paste(result_path, "secreted_proteins_P_value.svg",sep = ""), p_value_line, height = 5, width = 15)
    }
    p_value()

    gene_radio <- function() {
    BgRatio_result_sub <- subset(enrich_result, result_list_BgRatio > general_BgRatio)
    BgRatio_line <- ggplot(
        data = enrich_result,
        aes(
        x = as.numeric(rownames(enrich_result)),
        y = result_list_BgRatio
        )
    ) +
        geom_line(group = 1) +
        geom_point(
        data = BgRatio_result_sub,
        aes(
            x = as.numeric(rownames(BgRatio_result_sub)),
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
    ggsave(paste(result_path, "secreted_proteins_BgRatio_line.png",sep = ""), BgRatio_line, height = 5, width = 15)
    ggsave(paste(result_path, "secreted_proteins_BgRatio_line.svg",sep = ""), BgRatio_line, height = 5, width = 15)
    }
    gene_radio()
    '''
    robjects.globalenv["pav_0_1_tsv_file_name"] = pav_0_1_tsv_file_name
    robjects.globalenv["new_secreted_protein_file_name"] = new_secreted_protein_file_name
    robjects.globalenv["classic_secreted_protein_file_name"] = classic_secreted_protein_file_name
    robjects.globalenv["result_path"] = R_result
    robjects.globalenv["known_secreted_proteins"] = robjects.StrVector(known_secreted_proteins)
    robjects.r(R_code_draw_secreted_protein)
