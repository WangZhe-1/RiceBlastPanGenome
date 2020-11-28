import rpy2.robjects as robjects
from clade_color_5clade import create_color_clade
from rpy2.robjects.vectors import StrVector
R_code_set_minus_with_cluster='''
require(readxl)
require(tidyverse)
in_fl_raw=read_xlsx(
  paste(result_R,"set_minus_set_minus.xlsx",sep = "")
)
in_fl_raw$species=factor(in_fl_raw$species,levels = in_fl_raw$species)
cluster_category=read.table(cluster_file_name,header = T)
in_fl=merge(in_fl_raw,cluster_category,by.x = 1,by.y = 1,sort = F)

in_fl$cluster <- factor(in_fl$cluster , levels = order_levels)
in_fl_sorted=in_fl[order(in_fl$cluster), ]
in_fl_sorted$species=factor(in_fl_sorted$species,levels = in_fl_sorted$species)
row.names(in_fl_sorted)=NULL

draw_curve=function(input_data,plot_name){
  plot_line=ggplot(
    data = input_data,
    mapping=aes(
      x=as.numeric(row.names(input_data)),
      # y=minus
      group=1
    )
  )+
    geom_line(
      aes(
        y=minus
      ),
      size=1
    )+
    geom_line(
      aes(
        y=add,
      ),
      size=1
    )+
    # scale_y_continuous(
    #   breaks = c(8511,9000,10000,11000,12000,12110),
    #   labels = c(8511,9000,10000,11000,12000,12110)
    # )+
    # geom_area(
    #   aes(
    #     y=add,
    #     fill=factor(cluster)
    #     ),
    #   alpha=0.4
    # )+
  # geom_area(
  #   aes(colour = factor(cluster), fill = factor(cluster),y=add)
  #           )+
  # geom_tile(aes(y=add,fill = factor(cluster)))+
  geom_rect(
    aes(
      xmin=as.numeric(row.names(input_data))-0.5,
      xmax=as.numeric(row.names(input_data))+0.5,
      ymax=add,
      fill=cluster
    ),
    ymin=0
  )+
    # scale_fill_manual(values = alpha(brewer.pal(4,"Set2"), .8))+
    scale_fill_manual(values=alpha(color_clade$color, .8), 
                      name="Clade",
                      breaks=color_clade$clade,
                      # labels=c("Control", "Treatment 1", "Treatment 2")
                      )+
    theme(
      panel.background = element_blank(),
      axis.line = element_line(size=1, colour = "black"),
      axis.title = element_text(size=20),
      axis.text = element_text(size = 20)
    )+
    xlab(label = "Number of strain")+
    ylab(label = "Number of gene")+
    geom_hline(
      yintercept = c(range(in_fl$add)[2],range(input_data$minus)[1]),lty=3,lwd=1,alpha=0.8
    )+
    annotate(
      "text",
      x = 150,
      y = range(input_data$minus)[1]+1200,
      label="Core",
      size=10
    )+
    annotate(
      "text",
      x = 150,
      y = range(input_data$add)[2]-1000,
      label="Pan",
      size=10
    )
  ggsave(
    paste(result_R,plot_name,sep = ""),
    plot = plot_line,
    width = 8,
    height = 5
  )
}
'''
def set_minus_with_cluster(result_R,cluster_file_name,clade_set,color_set,color_num,order_levels,file_name_suffix):
  '''
  input 1: set_minus_set_minus.xlsx in result_R
  input 2: cluster_file_name
  input 3: clade_set
  input 4: color_set
  input 5: color_num 
  input 6: order_levels 决定几号对应的位置
  output in result_R directory: png, svg plot
  '''
  robjects.globalenv["result_R"]=result_R
  robjects.globalenv["color_clade"]=create_color_clade(clade_set,color_set,color_num)
  robjects.globalenv["cluster_file_name"]=cluster_file_name
  robjects.globalenv["order_levels"]=StrVector(order_levels)
  robjects.r(R_code_set_minus_with_cluster)
  robjects.r['draw_curve'](robjects.r['in_fl'],"set_minus_clade_mix_"+str(color_num)+"_binary.svg")
  robjects.r['draw_curve'](robjects.r['in_fl'],"set_minus_clade_mix_"+str(color_num)+"_binary.png")
  robjects.r['draw_curve'](robjects.r['in_fl_sorted'],"set_minus_clade_break_clearly_"+str(color_num)+"_"+file_name_suffix+".svg")
  robjects.r['draw_curve'](robjects.r['in_fl_sorted'],"set_minus_clade_break_clearly_"+str(color_num)+"_"+file_name_suffix+".png")
