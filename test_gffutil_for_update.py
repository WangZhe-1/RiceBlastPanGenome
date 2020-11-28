'''
@Author: your name
@Date: 2020-08-01 22:22:00
LastEditTime: 2020-09-23 18:51:46
LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/test_gffutil.py
'''
import gffutils
def iterator(db):
    for gene in db.features_of_type('gene'):

        # modify attributes

        # add a new attribute for exon id
        for target_id in gene.attributes["Target"][0].split(","):
            gene.attributes['ID'] = target_id.split(" ")[0]+"_"+gene.chrom
            yield gene


gffutils.create_db(
    "../Pan_genome_data_2/blast_intermediate_out/gth_out/1106.2.gff", 
    dbfn='../Pan_genome_data_2/blast_intermediate_out/1106.2.db', 
    force=True,
    keep_order=True,
    id_spec=None
    )
db = gffutils.FeatureDB('../Pan_genome_data_2/blast_intermediate_out/1106.2.db', keep_order=True)
# print(db["gene1"])
db.update(iterator(db),merge_strategy="create_unique")
print(db["FR13_UCOG01000026.1_gene_14_JMKE01000001.1"])