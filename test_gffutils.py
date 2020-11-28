import gffutils
import pybedtools
from gffutils.helpers import asinterval
from pybedtools.helpers import chromsizes
import pysam
import json
import re
from Bio import SeqIO
# gffutils.create_db("../../GFF/13FM-16-1.gff","../Pan_genome_data/13FM-16-1.gff.db",id_spec=None,force=True)
# dfs,afas=re.search("(.+)_protein_(.+)_","13FM-24-1_protein_6136_156aa").group(1,2)
# gffutils.create_db("../../70-15_refference_genome/70-15_Gff/magnaporthe_oryzae_70-15_8_genome_summary_per_gene_amend.txt","../Pan_genome_data/fdsa.db",id_spec=':source:',force=True)
# a = pybedtools.BedTool("../../70-15_refference_genome/70-15_Gff/magnaporthe_oryzae_70-15_8_genome_summary_per_gene.txt")
# a_re=a.remove_invalid()
# print(a_re[1])
# print(a_re["MGG_16078"])
db=gffutils.FeatureDB("../Pan_genome_data/13FM-16-1.gff.db")
# db=gffutils.FeatureDB("../Pan_genome_data/fdsa.db")
# print(db["MGG_16078"])
# print(db["MGG_16078"].strand)
# print(db["MGG_16078"].start)
# db=gffutils.FeatureDB("../Pan_genome_data/13FM-16-1.gff.db")
# print(gffutils.constants._gffkeys)
# for dfs in db.all_features():
#     print(dfs)
print(db['gene_2'])
# print (featuretype)
# id_spec={'gene': 'gene_id'
# id_spec=[':seqid:']
# dm3_chromsizes = chromsizes("../../contig/13FM-16-1.fasta")
# for i in sorted(dm3_chromsizes.items()):
#     print(i)



dsfa=pysam.FastaFile("../../contig/13FM-16-1.fasta")
# for sdf in dsfa:
#     print (sdf)
# print(dsfa.lengths)
# print(dsfa.fetch(reference="MQNR01000010.1").len)
with open('../Pan_genome_data/contig_length_json/13FM-16-1.json', 'r') as f:
    json_file=json.load(f)
# cov_chrsizes = a.genome_coverage(bg=True, genome=chromsizes)
# x=pybedtools.BedTool(asinterval(db['gene_2']))
x_feature=db['gene_2']
def fsdu(which_id):
    yield asinterval(db[which_id])
x_fun = pybedtools.BedTool(fsdu('gene_2'))
print("{}\t{}\t{}\n".format(x_feature.chrom,x_feature.start,x_feature.stop))
x_interval_string="{}\t{}\t{}\n".format(x_feature.chrom,x_feature.start,x_feature.stop)
x_interval = pybedtools.BedTool(x_interval_string, from_string=True)
x_interval_0 = pybedtools.BedTool(x_interval_string, from_string=True)[0]
x_interval_list=pybedtools.create_interval_from_list(str(db['gene_2']).split("\t"))

# x_iter=iter(x)
# extended_genes_list = x_interval_list.flank(b=1000,genome=json_file)
# extended_genes = x_fun.flank(b=1000,genome=json_file)
# extended_genes.sequence(fi="../../contig/13FM-16-1.fasta", fo="../Pan_genome_data/test_2.fasta", s=True)
extended_genes_slop = x_fun.slop(b=1000,genome=json_file,s=True)
sequence=extended_genes_slop.sequence(fi="../../contig/13FM-16-1.fasta",s=True)
print(open(sequence.seqfn).read())
SeqIO.read(sequence.seqfn)
print(sequence)
# for dsd,dsf in zip(dsfa.lengths,dsfa.references):
#     print(dsd)
#     print(dsf)
a = pybedtools.example_bedtool('a.bed')
print(a[0])
print(a.flank(genome='hg19', b=100)) #doctest: +NORMALIZE_WHITESPACE