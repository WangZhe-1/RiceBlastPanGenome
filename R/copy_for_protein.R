
require(tidyverse)
in_fl <- read.table("../Pan_genome_data/a_qc_contig_busco_read_busco.txt")
in_fl_95 <- in_fl %>%

    filter(V2 >= 95)
pb_protein_id=in_fl_95[grep("_PB",in_fl_95$V1,ignore.case = F),1]
need_remove=gsub("_PB","",pb_protein_id)
in_fl_95_removed_pb=in_fl_95[!(in_fl_95$V1%in%need_remove),]