create_color_clade=function(clade_set,color_set){
  require(RColorBrewer)
  return(data.frame(clade=I(clade_set),color=brewer.pal(4,color_set)))
}
dsaf=create_color_clade(parse(text="as.character(1:4)"),"Set2")
ccceval(parse(text="5+5"))
