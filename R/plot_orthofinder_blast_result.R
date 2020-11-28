require(tidyverse)
blast_value=read.table(blast_per_MGG_file_name,header = F,sep = '\t')
blast_same_strain=read.table(blast_same_strain_file_name,header = F,sep = '\t')
blast_line=ggplot(
  data = blast_value,
  aes(x=V2)
  )+
  geom_line(
    aes(
      y=V3
    ),
    group=1
  )+
  geom_line(
    aes(
      y=V4
    ),
    group=1
  )+
  # geom_rect(
  #   data = blast_same_strain,
  #   aes(
  #     xmin=as.numeric(V1)-10,
  #     xmax=as.numeric(V1)+10,
  #   ),
  #   ymax=Inf,
  #   ymin=0,
  #   fill="red",
  #   color="black"
  # )+
  # geom_tile(
  #   data = blast_same_strain,
  #   aes(
  #     x=V1,
  #     y=0.5
  #   ),
  #   height=1,
  #   fill="red",
  #   color="black",
  #   alpha=0.5
  # )+
  theme(
    axis.text.x = element_text(angle = 90)
  )
print(blast_line)
