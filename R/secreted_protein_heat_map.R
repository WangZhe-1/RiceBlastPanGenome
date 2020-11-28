require(pheatmap)
require(readxl)
require(ggplot2)
new_tpye_secreted_protein=read.table(new_tpye_secreted_protein_file_name)
classic_secreted_protein=read.table(classic_secreted_protein_file_name)
all_secreted_protein=rbind(new_tpye_secreted_protein,classic_secreted_protein)
pav=readxl::read_xlsx(pav_orthofinder_file_name)
secreted_protein_pav=pav[pav$protein_id %in% all_secreted_protein$V1,]
secreted_protein_pav=secreted_protein_pav[,-c(1:2)]
# heatmap=pheatmap(
#   secreted_protein_pav[,1:156],
#   scale="none",
#   color = c("blue","red"),
#   #color =brewer.pal(5,"RdYlBu")[c(5,1)],
#   # breaks = c(-1,0,4),
#   cluster_rows = F,
#   cluster_cols=T,
#   show_rownames=F,
#   show_colnames=F,
#   border=FALSE,
#   treeheight_col=0,
#   #main = "PAV Matrix",
#   legend=F,
#   silent = T,
#   gaps_row = c(10,11),
#   # cellheight = 1
#   #legend_breaks = c(1,4),
#   #legend_labels = c("Partial","Indel")
# )
# print(heatmap)
# ggsave(plot = heatmap,filename =paste(result_path, "heatmap.png",sep = ""))
# ,height = 50,limitsize = F
# he=rep("yellow",nrow(secreted_protein_pav))

row_labels = structure(c("avr"), names = c("MGG_18041"))

row_labels = structure(ifelse(), names = secreted_protein_pav$protein_id)
whar=function(r){
  return (ifelse(r=="MGG_18041","avr",""))
}
wh=sapply(secreted_protein_pav$protein_id, whar)
wh=rep("df",nrow(secreted_protein_pav))
ha = rowAnnotation(foo = anno_mark(at = match("MGG_18041T0",secreted_protein_pav$protein_id), labels = c("arv")))
require(ComplexHeatmap)
ht <- Heatmap(
  secreted_protein_pav[,1:156],
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  right_annotation = ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
              if(secreted_protein_pav[i,157]=="MGG_18041T0")
                grid.rect(x = x, y = y, width = width, height = unit(3, "mm"),
                        gp = gpar(col = "grey", fill = fill))
            }
  )
draw(ht)
match("MGG_18041T0",secreted_protein_pav$protein_id)
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))
set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
require(circlize)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

small_mat = mat[1:9, 1:9]
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
Heatmap(small_mat, name = "mat", col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          print(height)
          grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
# Heatmap(small_mat, name = "mat",  col = col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           if(small_mat[i, j] > 0)
#             grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# 
# Heatmap(small_mat, name = "mat",  col = col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           height=ifelse(rownames(small_mat)[i]=="row4",unit(66, "mm"),unit(2, "mm"))
#         })
# Heatmap(small_mat, name = "mat",  col = col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           if(rownames(small_mat)[i]=="row4")
#             grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# Heatmap(small_mat, name = "mat",  col = col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           if(rownames(small_mat)[i]=="row4")
#             grid.rect(x = x, y = y, width = width, height = unit(20, "mm"), 
#                     gp = gpar(col = "grey", fill = fill))
#         })
# "row9" %in% rownames(small_mat)
# rownames(8)
# rownames[8]
# rownames(small_mat)[4]

m = matrix(rnorm(1000), nrow = 100)
ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = month.name[1:10]))
Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha)
