require(RColorBrewer)
in_fl_raw=read_xlsx(
  paste(result_R,"set_minus_set_minus.xlsx",sep = "")
)
in_fl_raw$species=factor(in_fl_raw$species,levels = in_fl_raw$species)
cluster_category=read.table(cluster_file_name,header = T)
in_fl=merge(in_fl_raw,cluster_category,by.x = 1,by.y = 1,sort = F)

in_fl$cluster <- factor(in_fl$cluster , levels = c("4","1","2","3"))
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
  print(plot_line)
  # ggsave(
  #   paste(result_R,plot_name,sep = ""),
  #   plot = plot_line,
  #   width = 8,
  #   height = 5
  # )
}
draw_curve(in_fl,"safdklj.png")
draw_curve(in_fl_sorted,"dfsa")
