import rpy2.robjects as robjects
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

def extract_strain_id(busco_result):
    '''
    provide busco_result as input
    input:busco_result
    return:>95% busco result list
    '''
    R_extract_95_code="""
    extract_95 <- function(busco_result) {
        require(tidyverse)
        in_fl <- read.table(busco_result)
        in_fl_95 <- in_fl %>%
        filter(V2 >= 95)
        # print(in_fl_95$V1)
        pb_protein_id=in_fl_95[grep("_PB",in_fl_95$V1,ignore.case = F),1]
        need_remove=gsub("_PB","",pb_protein_id)
        in_fl_95_removed_pb=in_fl_95[!(in_fl_95$V1%in%need_remove),]
        # print(in_fl_95_removed_pb$V1)
        return(as.character(in_fl_95_removed_pb$V1))
    }
    """
    R_extract_95 = SignatureTranslatedAnonymousPackage(R_extract_95_code, "R_extract_95")
    R_95_strain=R_extract_95.extract_95(busco_result)
    strain_95_list=list(R_95_strain)
    strain_95_list.remove("70-15")
    strain_95_list.remove("HO_busco")
    strain_95_list.remove("PH42_busco")
    strain_95_list.append("magnaporthe_oryzae_70-15_8_proteins_T0")
    strain_95_list.append("HO")
    strain_95_list.append("PH42")
    return strain_95_list
