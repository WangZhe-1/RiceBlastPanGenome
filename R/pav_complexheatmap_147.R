pav_name="pav_orthofinder.xlsx"

require(readxl)
require(ComplexHeatmap)
pav_df=read_xlsx(pav_name)
pdf("heatmap_raster_by_png.pdf", width = 8, height = 8)
# ht=ComplexHeatmap::Heatmap(pav_df[14000:nrow(pav_df),3:158],  raster_device = "png",column_km = 3)
ht=ComplexHeatmap::Heatmap(pav_df[10000:nrow(pav_df),3:158],column_km = 3,raster_device = "png")
draw(ht)
dev.off()
