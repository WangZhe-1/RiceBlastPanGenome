from pathlib import Path
import re
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

def annotation_eggnog_cog(egg_cog_input_file_name,pan_categories_path_name,R_result):
    '''
    input 1: egg_cog_input_file_name
    input 2: pan categories path name
    output 1: enrichment analysis result and picture
    '''

    robjects.globalenv["in_fl_name"]=egg_cog_input_file_name
    robjects.globalenv["result_path"]=R_result
    pan_categories_path=Path(pan_categories_path_name)
    for category_file in pan_categories_path.iterdir():
        robjects.globalenv[category_file.stem+"_id"] = robjects.r['read.table'](str(category_file),sep="\t",**{'stringsAsFactors': False})

    run_enrich_code="""
    library(openxlsx)
    library(tidyverse)
    library(clusterProfiler)
    library(RColorBrewer)
    library(reshape)

    in_fl <- read.xlsx(
    in_fl_name,
    startRow = 5,
    colNames = T
    )
    enrich_result_list<<-vector("list", 6)
    i<<-1
    des <- data.frame(unique(in_fl$description))
    des <- na.omit(des)
    des <- des %>%
    rownames_to_column("id")
    colnames(des) <- c("id", "description")
    gene2desID <- na.omit(in_fl[, c(1, 22)])
    gene2desID <- merge(gene2desID, des, by.x = 2, by.y = 2, all.x = T)
    gene2desID <- gene2desID[, -1]
    gene2desID <- gene2desID[, c(2, 1)]

    counter <- function(input, nameS) {
    sum_result <- aggregate(Count ~ Description, data = input[, c(2, 9)], sum)
    colnames(sum_result) <- c("Description", nameS)
    return(sum_result)
    }
    fun_enrich <- function(gene_taxonomy, annotation_case, levelS) {
    gene_set <- intersect(gene_taxonomy$V1, annotation_case$Query)
    enrich_result <- enricher(gene_set,
        TERM2GENE = annotation_case,
        TERM2NAME = levelS,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        universe = gene2desID$Query,
        minGSSize = 1,
        maxGSSize = 5000
    )
    # return(counter(enrich_result@result))
    # write.table(enrich_result@result, paste(result_path,sprintf("%s_%s_enrich_result.txt", as.character(substitute(gene_taxonomy)), "cog"),sep = ""), append = F, col.names = T, row.names = F, sep = "\t", quote = F)
    enrich_result_list[[i]] <<- enrich_result@result
    names(enrich_result_list)[i]<<-as.character(substitute(gene_taxonomy))
    i <<- i + 1
    assign(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy))),sum(enrich_result@result$Count),envir = .GlobalEnv)
    return(enrich_result)
    }
    run <- function(kind, level_name) {
    assign(sprintf("core_%s_%s_enrich", kind, level_name), fun_enrich(core_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("core_%s_%s_count", kind, level_name), counter(get(sprintf("core_%s_%s_enrich", kind, level_name))@result, "core"), envir = .GlobalEnv)
    assign(sprintf("soft_core_%s_%s_enrich", kind, level_name), fun_enrich(soft_core_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("soft_core_%s_%s_count", kind, level_name), counter(get(sprintf("soft_core_%s_%s_enrich", kind, level_name))@result, "soft_core"), envir = .GlobalEnv)
    assign(sprintf("middle_%s_%s_enrich", kind, level_name), fun_enrich(middle_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("middle_%s_%s_count", kind, level_name), counter(get(sprintf("middle_%s_%s_enrich", kind, level_name))@result, "middle"), envir = .GlobalEnv)
    assign(sprintf("soft_specific_%s_%s_enrich", kind, level_name), fun_enrich(soft_specific_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("soft_specific_%s_%s_count", kind, level_name), counter(get(sprintf("soft_specific_%s_%s_enrich", kind, level_name))@result, "soft_specific"), envir = .GlobalEnv)
    assign(sprintf("specific_%s_%s_enrich", kind, level_name), fun_enrich(specific_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("specific_%s_%s_count", kind, level_name), counter(get(sprintf("specific_%s_%s_enrich", kind, level_name))@result, "specific"), envir = .GlobalEnv)
    # assign(sprintf("no_present_%s_%s_enrich", kind, level_name), fun_enrich(no_present_id, gene2desID, des), envir = .GlobalEnv)
    # assign(sprintf("no_present_%s_%s_count", kind, level_name), counter(get(sprintf("no_present_%s_%s_enrich", kind, level_name))@result, "no_present"), envir = .GlobalEnv)
    assign(sprintf("pan_new_%s_%s_enrich", kind, level_name), fun_enrich(pan_new_id, gene2desID, des), envir = .GlobalEnv)
    assign(sprintf("pan_new_%s_%s_count", kind, level_name), counter(get(sprintf("pan_new_%s_%s_enrich", kind, level_name))@result, "pan_new"), envir = .GlobalEnv)
    WriteXLS::WriteXLS(
    enrich_result_list,
    paste(result_path,'annotation_eggnog_cog_enrich_result.xlsx',sep = ""),
    # SheetNames=names(enrich_result_list),
    col.names = T,
    row.names = F
    )
    print(sum(core_id_sum_enrich,soft_core_id_sum_enrich,middle_id_sum_enrich,soft_specific_id_sum_enrich,specific_id_sum_enrich))
    }

    drawer <- function(main_order, other, main_id) {
    main_order_index <- sort(main_order[[2]], index.return = T, decreasing = TRUE)
    main_order <- main_order[c(0, main_order_index$ix + 0), ]
    main_order_10 <- main_order[1:10, ]

    for (i in 1:length(other)) {
        main_order_10 <- merge(main_order_10, other[[i]], by.x = 1, by.y = 1, all.x = T)
    }
    main_order_10[is.na(main_order_10)] <- 0
    Per <- (as.matrix(sapply(main_order_10[1:nrow(main_order_10), main_id], as.numeric))) / as.matrix(rowSums(sapply(main_order_10[1:nrow(main_order_10), 2:ncol(main_order_10)], as.numeric)))
    main_order_10_in <- sort(as.numeric(Per), index.return = T, decreasing = TRUE)
    main_order_10 <- main_order_10[c(0, main_order_10_in$ix + 0), ]
    main_order_10$Description <- factor(main_order_10$Description, levels = main_order_10$Description)

    main_order_10_melt <- melt(main_order_10, id.vars = "Description")
    main_order_10_melt$variable <- factor(main_order_10_melt$variable, levels = c("core", "soft_core", "middle", "soft_specific", "specific")[5:1])
    p <- ggplot(data = main_order_10_melt, aes(Description, as.numeric(value), fill = variable)) +
        geom_bar(stat = "identity", position = "fill", color = "black", width = 0.5, size = 0.25) +
        scale_fill_manual(values = brewer.pal(5, "RdYlGn")) +
        coord_flip() +
        # scale_y_continuous(labels = percent_format())+
        scale_y_continuous(
        # labels = percent
        expand = c(0, 0, 0, 0)
        ) +
        labs(
        x = "",
        y = "percent",
        title = ""
        ) +
        theme(
        axis.title = element_text(size = 15, face = "plain", color = "black"),
        axis.text = element_text(size = 12, face = "plain", color = "black"),
        legend.title = element_text(size = 14, face = "plain", color = "black"),
        legend.position = "right",
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.x = element_line(size = 0.5, color = "black"),
        axis.ticks = element_blank()
        # panel.grid =element_blank()
        )
    print(p)
    return(p)
    }

    run("pathway", "C")
    plot_core_pathway_C <- drawer(core_pathway_C_count, list(soft_core_pathway_C_count, middle_pathway_C_count, soft_specific_pathway_C_count, specific_pathway_C_count), 2)
    ggsave(paste(result_path,"annotation_cog_core.png",sep = ""), plot_core_pathway_C, width = 10, height = 5)
    plot_specific_pathway_C <- drawer(specific_pathway_C_count, list(soft_core_pathway_C_count, middle_pathway_C_count, soft_specific_pathway_C_count, core_pathway_C_count), 2)
    ggsave(paste(result_path,"annotation_cog_specific.png",sep = ""), plot_specific_pathway_C, width = 10, height = 5)
    """

    run_enrich = robjects.r(run_enrich_code)
