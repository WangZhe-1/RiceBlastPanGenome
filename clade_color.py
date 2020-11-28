from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.vectors import StrVector

R_code_clade_color='''
create_color_clade=function(clade_set,color_set){
  require(RColorBrewer)
  return(data.frame(clade=I(clade_set),color=brewer.pal(4,color_set)))
}
'''

R_clade_color = SignatureTranslatedAnonymousPackage(R_code_clade_color, "R_clade_color")

def create_color_clade(clade_set,color_set):
    '''
    input 1: clade_set, clade 1;clade 2;clade 3;clade 4
        used as break,not lable, must be identitiy with cutree result

    input 2: color_set, https://colorbrewer2.org/

    output: R data.frame clade vs color
    '''
    return R_clade_color.create_color_clade(StrVector(clade_set),color_set)