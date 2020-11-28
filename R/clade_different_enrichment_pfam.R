input_pfam_df=read.table("../../Pan_genome_data_2/clade_different_enrichment/pfam_annotation.txt",header = T)
input_pfam_df$geneID=sub("T0$","",input_pfam_df$geneID)
`clade_1&clade_2_distinct`=read.table("../../Pan_genome_data_2/combination_matrix_distinct_core/rice/clade_1&clade_2_distinct.txt")
clade_combination_set=c("clade_1&clade_2")
result_path="../../Pan_genome_data_2/clade_different_enrichment/"
counter=function(input,nameS){
  sum_result=aggregate(Count ~ Description, data = input[,c(2,9)], sum)
  colnames(sum_result)=c("Description",nameS)
  return(sum_result)
}
fun_enrich=function(clade_comb,annotation_case,levelS,annotation){
  require(clusterProfiler)
  # print(head(gene_taxonomy))
  # print(head(annotation_case))
  # print(annotation_case)
  # print(str(gene_taxonomy))
  # print(str(annotation_case))
  gene_taxonomy=get(paste(clade_comb,"_distinct",sep=""))
  gene_set=intersect(gene_taxonomy$V1,annotation_case$geneID)
  enrich_result=enricher(gene_set,
                         TERM2GENE = annotation_case,
                         TERM2NAME = levelS,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1, 
                         universe = annotation_case$geneID,
                         minGSSize = 1,
                         maxGSSize=5000
  )
  # return(counter(enrich_result@result))
  enrich_result_list[[i]] <<- enrich_result@result
  names(enrich_result_list)[i]<<-clade_comb
  i <<- i + 1
  # write.table(
  # enrich_result@result,file = paste(result_path,sprintf('%s_%s_enrich_result.txt',as.character(substitute(gene_taxonomy)),annotation),sep = ""), row.names = F, col.names = T,append=FALSE,quote = F,sep="\t"
  # )
  assign(sprintf('%s_sum_enrich',clade_comb),sum(enrich_result@result$Count),envir = .GlobalEnv)
  # print("name")
  # print(deparse(substitute(gene_taxonomy)))
  # print(str(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy)))))
  # print(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy))))
  # list.append(enrich_result_list,substitute(gene_taxonomy)=enrich_result@result)
  return(enrich_result)
}
run=function(kind,level_name){
  sum_enrich=0
  for (clade_comb in clade_combination_set){
    # assign(sprintf('no_present_%s_%s_enrich',kind,level_name),fun_enrich(no_present_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    # assign(sprintf('no_present_%s_%s_count',kind,level_name),counter(get(sprintf('no_present_%s_%s_enrich',kind,level_name))@result,"no_present"),envir = .GlobalEnv)
    assign(sprintf('%s_%s_%s_enrich',deparse(substitute(clade_comb)),kind,level_name),fun_enrich(clade_comb,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('%s_%s_%s_count',deparse(substitute(clade_comb)),kind,level_name),counter(get(sprintf('%s_%s_%s_enrich',deparse(substitute(clade_comb)),kind,level_name))@result,deparse(substitute(clade_comb))),envir = .GlobalEnv)
    sum_enrich=sum_enrich+get(sprintf('%s_sum_enrich',clade_comb))
  }
  # str(enrich_result_list)
  print(names(enrich_result_list))
  WriteXLS::WriteXLS(
    enrich_result_list,
    paste(result_path,sprintf('clade_different_%s_enrich_result.xlsx',annotation),sep = ""),
    # SheetNames=names(enrich_result_list),
    col.names = T,
    row.names = F
  )
  print(sum_enrich)
}

drawer=function(main_order,other,main_id){
  require(tidyverse)
  require(RColorBrewer)
  require(reshape)
  main_order_index=sort(main_order[[2]],index.return=T,decreasing = TRUE)
  main_order=main_order[c(0,main_order_index$ix+0),]
  main_order_10=main_order[1:10,]
  
  for(i in 1:length(other)){
    main_order_10=merge(main_order_10,other[[i]],by.x = 1,by.y=1,all.x=T)
  }
  main_order_10[is.na(main_order_10)]=0
  Per<-(as.matrix(sapply(main_order_10[1:nrow(main_order_10),main_id], as.numeric))) / as.matrix(rowSums(sapply(main_order_10[1:nrow(main_order_10),2:ncol(main_order_10)],as.numeric)))
  main_order_10_in<-sort(as.numeric(Per),index.return=T,decreasing = TRUE)                                                               
  main_order_10<-main_order_10[c(0,main_order_10_in$ix+0),] 
  main_order_10$Description=factor(main_order_10$Description,levels =main_order_10$Description)
  
  main_order_10_melt=melt(main_order_10,id.vars = "Description")
  main_order_10_melt$variable=factor(main_order_10_melt$variable,levels =c( "core","soft_core","middle","soft_specific","specific" )[5:1])
  p=ggplot(data=main_order_10_melt,aes(Description,as.numeric(value) ,fill=variable))+    
    geom_bar(stat="identity", position="fill",color="black", width=0.5,size=0.25)+
    scale_fill_manual(values=brewer.pal(5,"RdYlGn"))+
    coord_flip()+
    #scale_y_continuous(labels = percent_format())+
    scale_y_continuous(
      #labels = percent
      expand = c(0,0,0,0)
    )+
    labs(
      x="",
      y="percent",
      title = ""
    )+
    theme(
      axis.title=element_text(size=15,face="plain",color="black"),
      axis.text = element_text(size=12,face="plain",color="black"),
      legend.title=element_text(size=14,face="plain",color="black"),
      legend.position = "right",
      panel.background=element_rect(fill = "transparent",color=NA),
      panel.grid.major.x = element_line(size = 0.5,color = 'black'),
      axis.ticks = element_blank()
      #panel.grid =element_blank()
    )
  print(p)
  return(p)
}
run_enrich=function(in_df,annotation){
  enrich_result_list<<-list()
  i<<-1
  print(i)
  in_df<<-in_df
  annotation<<-annotation
  print(head(`clade_1&clade_2_distinct`))
  run("pathway","D")
  plot_core_pathway_D=drawer(core_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,specific_pathway_D_count),2)
  ggsave(paste(result_path,sprintf('clade_different_%s_core.png',annotation),sep = ""),plot_core_pathway_D,width = 10,height = 5)
  plot_specific_pathway_D=drawer(specific_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,core_pathway_D_count),2)
  ggsave(paste(result_path,sprintf('clade_different_%s_specific.png',annotation),sep = ""),plot_specific_pathway_D,width = 10,height = 5)
}
run_enrich(input_pfam_df,"pfam")

