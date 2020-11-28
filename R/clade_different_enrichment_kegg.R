`clade_1&clade_3_distinct`=read.table("../../Pan_genome_data_2/combination_matrix_distinct_core/rice/clade_1&clade_3_distinct.txt")
in_df=read.table("../../Pan_genome_data_2/clade_different_enrichment/kegg_annotation.txt",header = T)
clade_combination_set=c("clade_1&clade_3")


pathway_level=read.table("../../kegg_brite/hierarchy/pathway_merge.txt",sep = "\t",header = T)
brite_level=read.table("../../kegg_brite/hierarchy/brite_merge.txt",sep = "\t",header = T)
counter=function(input,nameS){
  sum_result=aggregate(Count ~ Description, data = input[,c(2,9)], sum)
  colnames(sum_result)=c("Description",nameS)
  return(sum_result)
}

fun_enrich=function(clade_comb,annotation_case,levelS){
  require(clusterProfiler)
  gene_taxonomy=get(paste(clade_comb,"_distinct",sep=""))
  gene_set=intersect(gene_taxonomy$V1,annotation_case$geneID)
  # str(annotation_case$geneID)
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
  enrich_level=merge(enrich_result@result,pathway_level,by.x=1,by.y = 1,all.x=T)
  enrich_level=merge(enrich_level,brite_level,by.x=1,by.y = 1,all.x=T)
  assign(sprintf('%s_sum_enrich',clade_comb),sum(enrich_result@result$Count),envir = .GlobalEnv)
  # write.table(enrich_level,paste(result_path,sprintf('kegg_%s_enrich_result.txt',as.character(substitute(gene_taxonomy))),sep = ""),append = F,col.names = T,row.names = F,sep = '\t',quote = F)
  enrich_result_list[[i]] <<- enrich_level
  names(enrich_result_list)[i]<<-clade_comb
  i <<- i + 1
  # print(str(enrich_result_list))
  # print(i)
  # print(names(enrich_result_list))
  return(enrich_result)
}

run=function(kind,level_name){
  # print(1)
  sum_enrich=0
  for (clade_comb in clade_combination_set){
    # assign(sprintf('no_present_%s_%s_enrich',kind,level_name),fun_enrich(no_present_id,in_df[,c(2,1)],in_df[,c(2,3)]),envir = .GlobalEnv)
    # assign(sprintf('no_present_%s_%s_count',kind,level_name),counter(get(sprintf('no_present_%s_%s_enrich',kind,level_name))@result,"no_present"),envir = .GlobalEnv)
    assign(sprintf('%s_%s_%s_enrich',clade_comb,kind,level_name),fun_enrich(clade_comb,in_df[,c(2,1)],in_df[,c(2,3)]),envir = .GlobalEnv)
    assign(
      sprintf('%s_%s_%s_count',clade_comb,kind,level_name),
      counter(get(sprintf('%s_%s_%s_enrich',clade_comb,kind,level_name))@result,
              clade_comb),
      envir = .GlobalEnv)
    
    sum_enrich=sum_enrich+get(sprintf('%s_sum_enrich',clade_comb))
  }
  # assign(sprintf('pan_new_%s_%s_enrich',kind,level_name),fun_enrich(pan_new_id,in_df[,c(2,1)],in_df[,c(2,3)]),envir = .GlobalEnv)
  # assign(sprintf('pan_new_%s_%s_count',kind,level_name),counter(get(sprintf('pan_new_%s_%s_enrich',kind,level_name))@result,"pan_new"),envir = .GlobalEnv)
  WriteXLS::WriteXLS(
    enrich_result_list,
    paste(result_path,'annotation_kegg_enrich_result.xlsx',sep = ""),
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
kegg=function(in_df){
  in_df$geneID=sub("T0$","",in_df$geneID)
  in_df<<-in_df
  enrich_result_list<<-list()
  i<<-1
  run("pathway","D")
  plot_core_pathway_D=drawer(core_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,specific_pathway_D_count),2)
  ggsave(paste(result_path,"annotation_pathway_core_D.png",sep = ""),plot_core_pathway_D,width = 10,height = 5)
  plot_specific_pathway_D=drawer(specific_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,core_pathway_D_count),2)
  ggsave(paste(result_path,"annotation_pathway_specific_D.png",sep = ""),plot_specific_pathway_D,width = 10,height = 5)
}
kegg(in_df)