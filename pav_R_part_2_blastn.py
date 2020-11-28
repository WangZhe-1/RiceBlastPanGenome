import rpy2.robjects as robjects

R_pav_code="""
library(tidyverse)
library(readxl)
library(data.table)
library(RColorBrewer)
library(WriteXLS)
library(stringr)
pav_df_raw=read_xlsx(pav_name)

pav_df=pav_df_raw[,-c(1:2)]

na_tsv=pav_df>=2
pav_df[na_tsv]=1
pav_df[!na_tsv]=0

pav_df_gene_name=cbind(pav_df,pav_df_raw$...2)
pav_tr=data.table::transpose(pav_df)
colnames(pav_tr)=pav_df_raw$...2
rownames(pav_tr)=colnames(pav_df)
count_present=colSums(pav_tr)
count_present=as.data.frame(count_present)
Classify=function(Min,Max,Name){
  core_df=count_present %>% 
    rownames_to_column("gene") %>% 
    filter(count_present>=Min & count_present<=Max) %>% 
    select(gene)
  write.table(
    core_df,file = paste(result_path,"/pan_categories_gene_blastn/",sprintf("%s.txt", Name),sep = ""), row.names = F, col.names = F,append=FALSE,quote = F
  )
  print(c(Min,Max))
  return(core_df)
}
spcies_num=nrow(pav_tr)
core=Classify(spcies_num,spcies_num,"core")
soft_core=Classify( ceiling(spcies_num*0.95),spcies_num-1,"soft_core")
middle=Classify( ceiling(spcies_num*0.05)+1, ceiling(spcies_num*0.95)-1,"middle")
soft_specific=Classify(2, ceiling(spcies_num*0.05),"soft_specific")
specific=Classify(1,1,"specific")
no_present=Classify(0,0,"no_present")

write_class=function(){
  class_num=rbind(
    core=nrow(core),
    soft_core=nrow(soft_core),
    middle=nrow(middle),
    soft_specific=nrow(soft_specific),
    specific=nrow(specific),
    no_present=nrow(no_present)
  )
  write.table(class_num,file = paste(result_path,"/pan_categories_gene_blastn/","class_num.txt",sep = ""),row.names = T,col.names = F,append = F,quote = F)
  pan_new_bool=!str_detect(colnames(pav_tr),"MGG")
  write.table(
    colnames(pav_tr)[pan_new_bool],file = paste(result_path,"/pan_categories_gene_blastn/","pan_new.txt",sep = ""), row.names = F, col.names = F,append=FALSE,quote = F
  )
}
heatmap_function=function(){
  pav_sort_index=sort(rowSums(pav_df),index.return=T,decreasing=T)
  pav_df_sorted=pav_df[pav_sort_index$ix,]
  heatmap=pheatmap(
    pav_df_sorted,
    scale="none",
    color = c("blue","red"),
    #color =brewer.pal(5,"RdYlBu")[c(5,1)],
    # breaks = c(-1,0,4),
    cluster_rows = F,
    cluster_cols=T,
    show_rownames=F,
    show_colnames=F,
    border=FALSE,
    #main = "PAV Matrix",
    legend=F,
    silent = T
    #legend_breaks = c(1,4),
    #legend_labels = c("Partial","Indel")
  )
  # ggsave(plot = heatmap,filename = paste(result_path,"heatmap.svg",sep = ""))
  ggsave(plot = heatmap,filename =paste(result_path, "heatmap_blastn.png",sep = ""))
}

draw_stack=function(){
  distributed=c(middle$gene,soft_specific$gene)
  
  core=cbind(core$gene,cat=rep("core",times=nrow(core)))
  soft_core=cbind(soft_core$gene,cat=rep("soft_core",times=nrow(soft_core)))
  distributed=cbind(distributed,cat=rep("distributed",times=length(distributed)))
  specific=cbind(specific$gene,cat=rep("strain_specific",times=nrow(specific)))
  category=as.data.frame(rbind(core,soft_core,distributed,specific),stringsAsFactors = F)
  pav_df_cat=merge(pav_df_gene_name,category,by.x = "pav_df_raw$...2",by.y = "V1")
  # pav_df_cat=pav_df_cat[,-c(1:2)]
  
  stat_result=data.frame(specise_name=character(0),category=character(0),num=numeric(0))
  core_sum=data.frame(specise_name=character(0),num=numeric(0))
  all_sum=data.frame(specise_name=character(0),num=numeric(0))
  
  for (i in 2:(ncol(pav_df_cat)-1)) {
    sha=aggregate(pav_df_cat[,i]~cat,pav_df_cat,sum)
    stat_result=rbind(stat_result, data.frame(specise_name=colnames(pav_df_cat)[i],category=as.character(sha[[1,1]]),num=sha[1,2]))
    stat_result=rbind(stat_result, data.frame(specise_name=colnames(pav_df_cat)[i],category=as.character(sha[[2,1]]),num=sha[2,2]))
    stat_result=rbind(stat_result, data.frame(specise_name=colnames(pav_df_cat)[i],category=as.character(sha[[3,1]]),num=sha[3,2]))
    stat_result=rbind(stat_result, data.frame(specise_name=colnames(pav_df_cat)[i],category=as.character(sha[[4,1]]),num=sha[4,2]))
    with(
      sha[sha$cat== "soft_core"|sha$cat=="core",],
      {
        core_sum <<- rbind(core_sum,data.frame(specise_name=colnames(pav_df_cat)[i],num=sum(`pav_df_cat[, i]`)))
      }
    )
    all_sum=rbind(all_sum,data.frame(specise_name=colnames(pav_df_cat)[i],num=sum(sha$`pav_df[, i]`)))
  }
  sum_sort=sort(core_sum$num,index.return=T,decreasing=T)
  core_sum=core_sum[sum_sort$ix,]
  stat_result$specise_name=factor(stat_result$specise_name,levels =core_sum$specise_name )
  stat_result$category=factor(stat_result$category,levels = c("strain_specific","distributed","soft_core","core"))
  
  bar_plot=ggplot(
    # data = stat_result_sub,
    data = stat_result,
    aes(
      x=specise_name,
      y=num,
      fill=category
    )
  )+
    geom_bar(stat="identity",position="stack")+
    scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    # scale_fill_manual(values=c("black",brewer.pal(5,"RdYlGn")[c(5,2:1)]))+
    scale_fill_manual(values=c("black",brewer.pal(3,"RdYlGn")))+
    labs(x="Strain",y="Number")+
    theme(
      panel.grid =element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(paste(result_path,"stack_bar_blastn.png",sep = ""),bar_plot)
  ggsave(paste(result_path,"stack_bar_blastn.svg",sep = ""),bar_plot)
}
"""
def R_run_pav_instance_maker_blastn(pav_excel,R_result):
    '''
    input 1: pav_name
    output 1: result_path_name
    creat a R instance and clasify pan genome into core,specific ... write them into file 
    '''
    robjects.globalenv["pav_name"] = pav_excel
    robjects.globalenv["result_path"]=R_result
    robjects.r(R_pav_code)
def write_class_blastn():
    '''
    write_class
    '''
    robjects.r['write_class']()
def heatmap_blastn ():
    '''
    heatmap function
    '''
    robjects.r['heatmap_function']()
def draw_stack_blastn():
    '''
    '''
    robjects.r["draw_stack"]()

