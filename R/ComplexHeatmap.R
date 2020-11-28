require(ComplexHeatmap)
random_mat = function(nr) {
  m = matrix(rnorm(10*nr), nc = 10)
  colnames(m) = letters[1:10]
  return(m)
}
y = NULL
for(nr in c(1, 20)) {
  ht = draw(Heatmap(random_mat(nr), height = unit(5, "mm")*nr, 
                    column_title = "foo", # one line text
                    top_annotation = HeatmapAnnotation(bar = 1:10)))
  ht_height = sum(component_height(ht)) + unit(4, "mm")
  ht_height = convertHeight(ht_height, "inch", valueOnly = TRUE)
  y = NA
}
