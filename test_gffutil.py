'''
@Author: your name
@Date: 2020-08-01 22:22:00
LastEditTime: 2020-09-23 15:42:21
LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/test_gffutil.py
'''
import gffutils
# gffutils.create_db(
#     "../Pan_genome_data_2/blast_intermediate_out/gth_out/1106.2.gff", 
#     dbfn='../Pan_genome_data_2/blast_intermediate_out/1106.2.db', 
#     force=True,
#     keep_order=True,
#     id_spec=None
#     )
db = gffutils.FeatureDB('../Pan_genome_data_2/blast_intermediate_out/1106.2.db', keep_order=True)
print(db["gene1"])