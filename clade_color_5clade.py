'''
@Author: your name
@Date: 2020-08-13 16:50:23
@LastEditTime: 2020-08-13 16:55:31
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/clade_color_5clade.py
'''
'''
@Author: your name
@Date: 2020-08-13 16:50:23
@LastEditTime: 2020-08-13 16:50:23
@LastEditors: your name
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/clade_color_5clade.py
'''
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.vectors import StrVector

R_code_clade_color='''
create_color_clade=function(clade_set,color_set,color_num){
  require(RColorBrewer)
  return(data.frame(clade=I(clade_set),color=brewer.pal(color_num,color_set)))
}
'''

R_clade_color = SignatureTranslatedAnonymousPackage(R_code_clade_color, "R_clade_color")

def create_color_clade(clade_set,color_set,color_num):
    '''
    input 1: clade_set, clade 1;clade 2;clade 3;clade 4
        used as break,not lable, must be identitiy with cutree result

    input 2: color_set, https://colorbrewer2.org/

    output: R data.frame clade vs color
    '''
    return R_clade_color.create_color_clade(StrVector(clade_set),color_set,color_num)