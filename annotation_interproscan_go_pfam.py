#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from pathlib import Path
import pickle
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
def annotation_go_pfam(annotation_go_pfam,pan_categories_path_name,R_result,go_out_file,pfam_out_file):
  '''
  input 1: interproscan result
  input 2: pan categories path name
  output 1: enrichment analysis result and picture
  output 2: go_out_file, go functional annotation output
  output 3: pfam_out_file，pfam functional annotation output
  瑕疵：
  判断基因名字的正则表达式，不能适应所有的基因名字
  应该使用位置判断，
  比如，第一列是基因名字

  检查：
  interpro gff->go、pfam list
  在gff文件中用正则表达式计数，和go、pfam list中元素总数对比

  list->enrich result
  划分
    各个部分（core，specific）的GeneRatio之和应等于BgRatio

  各个部分（core，specific）的count之和应等于go、pfam list中元素总数

  总：
  在gff文件中用正则表达式计数=各个部分（core，specific）的count之和
  '''
  ###create go dic
  dic_go_name={}
  dic_go_namespace={}
  rgx_go_id=re.compile("^id: (GO:\d+)")
  rgx_go_name=re.compile("^name: (.+)")
  rgx_go_namespace=re.compile("^namespace: (.+)")
  rgx_go_alt_id=re.compile("^alt_id: (GO:\d+)")

  with open("../annotation/go.obo") as fh_go_ob:
      for go_line in fh_go_ob:
          result_go_id=rgx_go_id.search(go_line)
          if result_go_id is not None:
              go_id_obo=result_go_id.group(1)
          result_go_name=rgx_go_name.search(go_line)
          if result_go_name is not None:
              go_name_obo=result_go_name.group(1)
              dic_go_name.setdefault(str(go_id_obo),str(go_name_obo))
          result_go_namespace=rgx_go_namespace.search(go_line)
          if result_go_namespace is not None:
              go_namespace_obo=result_go_namespace.group(1)
              dic_go_namespace.setdefault(str(go_id_obo),str(go_namespace_obo))
          result_go_alt_id=rgx_go_alt_id.search(go_line)
          if result_go_alt_id is not None:
              go_alt_id_obo=result_go_alt_id.group(1)
              dic_go_name.setdefault(str(go_alt_id_obo),str(go_name_obo))
              dic_go_namespace.setdefault(str(go_alt_id_obo),str(go_namespace_obo))
  ##################################################################################

  ###read gff to get go pfam
  rgx_gff_go=re.compile("GO:\d{7}")
  rgx_gff_pfam_id=re.compile("Name=(PF\d{5})")
  rgx_gene_1=re.compile("^(?P<TBX>\w.+?_protein_.+?aa)\s")
  rgx_gene_2=re.compile("^((?:mg)*g.+?)\s",re.IGNORECASE)
  rgx_gff_pfam_annotation=re.compile("signature_desc=(.+?);")

  dic_gff_go={}
  dic_gff_pfam={}
  r_pfam=2
  pfam_annotation_list=[]
  pfam_id_list=[]
  pfam_protein_id_list=[]
  go_name_list=[]
  go_namespace_list=[]
  go_protein_id_list=[]
  go_id_list=[]
  def get_gene_id(line):
      result_line_1=rgx_gene_1.search(line)
      result_line_2=rgx_gene_2.search(line)
      if result_line_1 is not None:
          return result_line_1.group('TBX')
      elif result_line_2 is not None:
          return result_line_2.group(1)
      else:
          return None
  with open(annotation_go_pfam,'r') as fl_gff:
      for gff_line in fl_gff:
          result_protein_id_gff=get_gene_id(gff_line)
          result_go_gff=rgx_gff_go.finditer(gff_line)
          result_pfam_id_gff=rgx_gff_pfam_id.search(gff_line)
          result_pfam_annotation_gff=rgx_gff_pfam_annotation.search(gff_line)
          if result_protein_id_gff is None:
              continue
          protein_id_gff=result_protein_id_gff
          if result_go_gff is not None:
              for go_id_gff_match in result_go_gff:
                  go_id_gff=go_id_gff_match.group(0)
                  if protein_id_gff in dic_gff_go:
                      dic_gff_go[protein_id_gff].add(go_id_gff)
                  else:
                      dic_gff_go.setdefault(protein_id_gff,set())
                      dic_gff_go[protein_id_gff].add(go_id_gff)
          if result_pfam_id_gff is not None:
              pfam_id_gff=result_pfam_id_gff.group(1)
              pfam_annotation_gff=result_pfam_annotation_gff.group(1)
              pfam_protein_id_list.append(protein_id_gff)
              pfam_annotation_list.append(pfam_annotation_gff)
              pfam_id_list.append(pfam_id_gff)

  for protein_id_go in dic_gff_go:
      for go_id_write in dic_gff_go[protein_id_go]:
          go_protein_id_list.append(protein_id_go)
          go_id_list.append(go_id_write)
          go_name_list.append(dic_go_name[go_id_write])
          go_namespace_list.append(dic_go_namespace[go_id_write])

  # with open("../eggnog_annotation/pfam_pfam_id_list.txt", 'w') as p_id_fl:
  #     p_id_fl.writelines("%s\n" % place for place in pfam_id_list)
  # with open("../eggnog_annotation/pfam_protein_id_list.txt", 'w') as p_protein_fl:
  #     p_protein_fl.writelines("%s\n" % place for place in pfam_protein_id_list)

  input_pfam_df=robjects.DataFrame({
      "geneID":robjects.StrVector(pfam_protein_id_list),
      "KOID":robjects.StrVector(pfam_id_list),
      "description":robjects.StrVector(pfam_annotation_list)
  })
  input_go_df=robjects.DataFrame({
      "geneID":robjects.StrVector(go_protein_id_list),
      "KOID":robjects.StrVector(go_id_list),
      "description":robjects.StrVector(go_name_list),
      "namespace":robjects.StrVector(go_namespace_list)
  })

  robjects.r['write.table'](input_pfam_df,pfam_out_file,**{'append': False},**{'row.names': False},**{'sep': "\t"})
  robjects.r['write.table'](input_go_df,go_out_file,**{'append': False},**{'row.names': False},**{'sep': "\t"})
  
  pan_categories_path=Path(pan_categories_path_name)
  for category_file in pan_categories_path.iterdir():
    robjects.globalenv[category_file.stem+"_id"] = robjects.r['read.table'](str(category_file),sep="\t",**{'stringsAsFactors': False})
  
  robjects.globalenv["result_path"]=R_result

  run_enrich_code="""
  counter=function(input,nameS){
    sum_result=aggregate(Count ~ Description, data = input[,c(2,9)], sum)
    colnames(sum_result)=c("Description",nameS)
    return(sum_result)
  }
  fun_enrich=function(gene_taxonomy,annotation_case,levelS,annotation){
    require(clusterProfiler)
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
    names(enrich_result_list)[i]<<-as.character(substitute(gene_taxonomy))
    i <<- i + 1
    # write.table(
    # enrich_result@result,file = paste(result_path,sprintf('%s_%s_enrich_result.txt',as.character(substitute(gene_taxonomy)),annotation),sep = ""), row.names = F, col.names = T,append=FALSE,quote = F,sep="\t"
    # )
    assign(sprintf('%s_sum_enrich',as.character(substitute(gene_taxonomy))),sum(enrich_result@result$Count),envir = .GlobalEnv)
    # list.append(enrich_result_list,substitute(gene_taxonomy)=enrich_result@result)
    return(enrich_result)
  }
  run=function(kind,level_name){
    assign(sprintf('core_%s_%s_enrich',kind,level_name),fun_enrich(core_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('core_%s_%s_count',kind,level_name),counter(get(sprintf('core_%s_%s_enrich',kind,level_name))@result,"core"),envir = .GlobalEnv)
    assign(sprintf('soft_core_%s_%s_enrich',kind,level_name),fun_enrich(soft_core_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('soft_core_%s_%s_count',kind,level_name),counter(get(sprintf('soft_core_%s_%s_enrich',kind,level_name))@result,"soft_core"),envir = .GlobalEnv)
    assign(sprintf('middle_%s_%s_enrich',kind,level_name),fun_enrich(middle_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('middle_%s_%s_count',kind,level_name),counter(get(sprintf('middle_%s_%s_enrich',kind,level_name))@result,"middle"),envir = .GlobalEnv)
    assign(sprintf('soft_specific_%s_%s_enrich',kind,level_name),fun_enrich(soft_specific_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('soft_specific_%s_%s_count',kind,level_name),counter(get(sprintf('soft_specific_%s_%s_enrich',kind,level_name))@result,"soft_specific"),envir = .GlobalEnv)
    assign(sprintf('specific_%s_%s_enrich',kind,level_name),fun_enrich(specific_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('specific_%s_%s_count',kind,level_name),counter(get(sprintf('specific_%s_%s_enrich',kind,level_name))@result,"specific"),envir = .GlobalEnv)
    assign(sprintf('no_present_%s_%s_enrich',kind,level_name),fun_enrich(no_present_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('no_present_%s_%s_count',kind,level_name),counter(get(sprintf('no_present_%s_%s_enrich',kind,level_name))@result,"no_present"),envir = .GlobalEnv)
    assign(sprintf('pan_new_%s_%s_enrich',kind,level_name),fun_enrich(pan_new_id,in_df[,c(2,1)],in_df[,c(2,3)],annotation),envir = .GlobalEnv)
    assign(sprintf('pan_new_%s_%s_count',kind,level_name),counter(get(sprintf('pan_new_%s_%s_enrich',kind,level_name))@result,"pan_new"),envir = .GlobalEnv)
    
    # str(enrich_result_list)
    print(names(enrich_result_list))
    WriteXLS::WriteXLS(
    enrich_result_list,
    paste(result_path,sprintf('annotation_%s_enrich_result.xlsx',annotation),sep = ""),
    # SheetNames=names(enrich_result_list),
    col.names = T,
    row.names = F
  )
    print(sum(core_id_sum_enrich,soft_core_id_sum_enrich,middle_id_sum_enrich,soft_specific_id_sum_enrich,specific_id_sum_enrich,no_present_id_sum_enrich))
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
  function(in_df,annotation){
    enrich_result_list<<-vector("list", 7)
    i<<-1
    head(core_id,10)
    in_df<<-in_df
    annotation<<-annotation
    run("pathway","D")
    plot_core_pathway_D=drawer(core_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,specific_pathway_D_count),2)
    ggsave(paste(result_path,sprintf('annotation_%s_core.png',annotation),sep = ""),plot_core_pathway_D,width = 10,height = 5)
    plot_specific_pathway_D=drawer(specific_pathway_D_count,list(soft_core_pathway_D_count,middle_pathway_D_count,soft_specific_pathway_D_count,core_pathway_D_count),2)
    ggsave(paste(result_path,sprintf('annotation_%s_specific.png',annotation),sep = ""),plot_specific_pathway_D,width = 10,height = 5)
  }
  """
  run_enrich = robjects.r(run_enrich_code)
  run_enrich(input_pfam_df,"pfam")
  run_enrich(input_go_df,"go")