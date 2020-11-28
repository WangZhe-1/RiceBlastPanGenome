'''
Author: your name
Date: 2020-10-02 11:59:10
LastEditTime: 2020-10-07 12:50:42
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/upset.py
'''
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

R_code='''
read_all_pav=function(pav_file_name,category_file_name){
    require(readxl)
    pav_raw=readxl::read_xlsx(pav_file_name)
    pav_df=pav_raw[,-c(1:2)]
    na_tsv=pav_df>=2
    pav_df[na_tsv]=1
    pav_df[!na_tsv]=0
    pav_logical=as.data.frame(apply(pav_df,2,as.logical))
    general_protein<<-pav_raw$...2
    strain_clade<<-read.table(category_file_name,header = T,stringsAsFactors = F)
    return(pav_logical)
}
read_rice_pav=function(
    pav_file_name,
    category_file_name,
    rice_pan_protein_file_name,
    rice_pav_file_name
    ){
    require(readxl)
    pav_raw=readxl::read_xlsx(pav_file_name)
    strain_clade<<-read.table(category_file_name,header = T,stringsAsFactors = F)
    rice_strainID=dplyr::filter(strain_clade,cluster %in% c("1","4","2"))
    pav_rice=pav_raw[rice_strainID$species]
    
    geneID_sub=sub("_.+","",pav_raw$...2)
    other_host_StrainID_df=dplyr::filter(strain_clade,cluster ==3)
    other_host_StrainID_need_remove_list=sapply(other_host_StrainID_df$species, function(x) grep(x,geneID_sub))
    other_host_StrainID_need_remove_unlist=unlist(other_host_StrainID_need_remove_list)
    pav_rice=pav_rice[-other_host_StrainID_need_remove_unlist,]
    general_protein=pav_raw$...2
    general_protein<<-general_protein[-other_host_StrainID_need_remove_unlist]
    general_protein=general_protein[-other_host_StrainID_need_remove_unlist]
    write.table(general_protein,rice_pan_protein_file_name,row.names = F,col.names = F,quote = F)
    write.table(pav_rice,rice_pav_file_name,row.names = F,col.names = T,quote = F,sep = "\t")

    na_tsv=pav_rice>=2
    pav_rice[na_tsv]=1
    pav_rice[!na_tsv]=0
    pav_logical=as.data.frame(apply(pav_rice,2,as.logical))
    return(pav_logical)
}

extract_comb_gene=function(which_comb,which_mode,comb_label,comb_matrix,out_dir){
  code=comb_label[which_comb]
  if (which_comb==""){
    code=unname(comb_label[1])
  }
  result=general_protein[extract_comb(comb_matrix,code)]
  # result=extract_comb(comb_matrix,comb_label[which_comb])
  write.table(result,paste(out_dir,which_comb,"_",which_mode,".txt",sep = ""),row.names = F,col.names = F,quote = F)
}
make_combination_matrix=function(pav_logical,zhong_zhe_merge,which_mode,clade_combination,out_dir){
    require(ComplexHeatmap)
    strain_clade=merge(strain_clade,zhong_zhe_merge,by.x = 2,by.y = 2)
    
    union_clade_element=function(clade_set){
        result=pav_logical[[clade_set[1]]]
        for (clade in clade_set) {
        result=result&pav_logical[[clade]]
        }
        return(result)
    }
    clade_binary_df_raw=aggregate(species ~ zhong, data = strain_clade, union_clade_element)
    
    datalist = list()
    for (i in 1:(nrow(zhong_zhe_merge))) {
        datalist[[i]]=clade_binary_df_raw$species[i,]
        names(datalist)[i]=as.character(clade_binary_df_raw$zhong[i])
    }
    clade_binary_df = do.call(cbind, datalist)
    comb_matrix=make_comb_mat(clade_binary_df,mode = which_mode)
    comb_label=structure(comb_name(comb_matrix), names = comb_name(comb_matrix,readable = T))
    sapply(clade_combination,function(x) extract_comb_gene(x,which_mode,comb_label,comb_matrix,out_dir))
    comb_size_readable=cbind(as.data.frame(structure(comb_name(comb_matrix,readable = T),names=comb_name(comb_matrix))),as.data.frame(comb_size(comb_matrix)))
    colnames(comb_size_readable)=c("name","count")
    require(tibble)
    comb_size_readable=tibble::rownames_to_column(comb_size_readable,"code")
    write.table(comb_size_readable,paste(out_dir,which_mode,"_combination_count.txt",sep = ""),quote = F,sep = "\t",col.names = T,row.names = F)
}
'''

R_upset_combination_matrix = SignatureTranslatedAnonymousPackage(R_code, "R_upset_combination_matrix")
def make_combination_matrix_all(pav_all_file_name,category_file_name,zhong_zhe_merge,which_mode,clade_combination,out_dir):
    '''
    input 1: all pav file
    input 2: strain_clade_file
    input 3: 振辉师兄的clade和自己的clade的对应关系
    input 4: which mode,"distinct", "intersect", "union"
    input 5: 你关心的clade combination
    output 1: output dir
    '''
    pav_logical=R_upset_combination_matrix.read_all_pav(pav_all_file_name,category_file_name)
    R_upset_combination_matrix.make_combination_matrix(
        pav_logical,
        zhong_zhe_merge,
        which_mode,
        clade_combination,
        out_dir
        )
def make_combination_matrix_rice(pav_all_file_name,category_file_name,zhong_zhe_merge,which_mode,clade_combination,out_dir,rice_pan_protein_file_name,rice_pav_file_name):
    '''
    input 1: all pav file
    input 2: strain_clade_file
    input 3: 振辉师兄的clade和自己的clade的对应关系
    input 4: which mode,"distinct", "intersect", "union"
    input 5: 你关心的clade combination
    output 1: output dir
    output 2: rice protein 
    output 3: rice pav
    '''
    pav_logical=R_upset_combination_matrix.read_rice_pav(
        pav_all_file_name,
        category_file_name,
        rice_pan_protein_file_name,
        rice_pav_file_name
        )
    R_upset_combination_matrix.make_combination_matrix(
        pav_logical,
        zhong_zhe_merge,
        which_mode,
        clade_combination,
        out_dir
        )
