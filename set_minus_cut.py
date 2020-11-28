from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
def set_minus_cut(start_point,pav_df_file_name,result_path):
    R_code_set_minus_cut='''
    Cut=function(start_point,pav_df_file_name,result_path){
    require(readxl)
    require(WriteXLS)
    require(tidyverse)
    pav_df=read_xlsx(pav_df_file_name)
    gene_is=pav_df %>% 
        filter((!!sym(start_point))==1) %>%
        # filter(`70-15`)==1
        select("protein_id") 
    
    minus_part=pav_df %>% 
        filter((!!sym(start_point))==1) %>% 
        column_to_rownames("protein_id")
    # pav_df_colsum=colSums(minus_part_num)
    # pav_df_colsum_sort=sort(pav_df_colsum)
    add_part=pav_df %>% 
        filter((!!sym(start_point))==0) %>% 
        column_to_rownames("protein_id")
    # add_part=add_part %>% 
    #   column_to_rownames(pav_df_raw$...2)
    add_part_num=sapply(add_part[2:157], function(x) as.numeric(x))
    pav_df_colsum=colSums(add_part_num)
    pav_df_colsum_sort=sort(pav_df_colsum)
    
    write.table(attributes(pav_df_colsum_sort),paste(result_path,sprintf("set_minus_sort_protein_id_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = F,col.names = F)
    write.table(pav_df_colsum_sort,paste(result_path,sprintf("set_minus_sort_protein_id_num_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = T,col.names = F)
    
    
    WriteXLS::WriteXLS(
        minus_part,
        paste(result_path,sprintf("set_minus_minus_%s.xlsx", start_point),sep = ""),
        col.names = T,
        row.names = T
    )
    WriteXLS::WriteXLS(
        add_part,
        paste(result_path,sprintf("set_minus_add_%s.xlsx", start_point),sep = ""),
        col.names = T,
        row.names = T
    )
    write.table(gene_is,paste(result_path,sprintf("set_minus_gene_id_%s.txt", start_point),sep = ""),append = F,quote = F,row.names = F,col.names = F)
    }
    '''
    R_set_minus_cut = SignatureTranslatedAnonymousPackage(R_code_set_minus_cut, "R_set_minus_cut")
    R_set_minus_cut.Cut(start_point,str(pav_df_file_name),result_path)
