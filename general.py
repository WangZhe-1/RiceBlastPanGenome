#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from rpy2.robjects.packages import importr
from extract_sequence_from_Orthogroups import extract_sequence_from_Orthogroups
from trim_N import sequence_cleaner_detect
from extract_nucleic import *
from pav_contig_present import blast
import rpy2.robjects as robjects
from pav_R_part import *
from set_minus import set_minus
from pathlib import Path
from set_minus import drawer_curve
from annotation_interproscan_go_pfam import annotation_go_pfam
from annotation_kegg_id import annotation_kegg
from annotation_eggnog_cog import annotation_eggnog_cog
from cluster import cluster
from cluster import write_cluster_result
from set_minus_with_cluster import set_minus_with_cluster
from annotation_secreted_proteins import secreted_protein
from annotation_secreted_proteins import drawer_secreted_protein
from alignment_extract_orthologs import extract_ortholog_gene
from Directory_creater import directory_creater
from set_minus_orthofinder import R_set_minus_cut_orthofinder
from blast_70_15_mgg_Augustus import blast
from blast_70_15_mgg_Augustus import read_blast_result
from blast_70_15_mgg_Augustus import blastdb
from blast_70_15_mgg_Augustus import MGG_unpresent_Augustus_locate_in_pav_orthofinder
from blast_70_15_mgg_Augustus import remove_MGG_unpresent_Augustus
from set_minus_cut import set_minus_cut
from alignment_extract_fst_prepare_contig_length import contig_length2json
from alignment_call_snp_HaplotypeCaller import call_snp_HaplotypeCaller
from alignment_call_snp_HaplotypeCaller import GenomicsDBImport
from alignment_call_snp_HaplotypeCaller import CombineGVCFs
from alignment_call_snp_HaplotypeCaller import GenotypeGVCFs
from alignment_call_snp_mummer import call_snp_mummer
from nucleotide_diversity import calculate_nucleotide_diversity
from calculate_dnds import calculate_dnds, paraAT, parse_ParaAT_result
from calculate_dnds import prepare_for_ParaAT

orthogroup_sequence_path="../Pan_genome_data/Results_May22_1/Orthogroup_Sequences/"
pan_protein_file="../Pan_genome_data/pan_protein.fasta"
pan_id_file="../Pan_genome_data/pan_id.txt"

# extract_sequence_from_Orthogroups(orthogroup_sequence_path,pan_protein_file,pan_id_file)

X_count_table_file="../Pan_genome_data/X_count_table.txt"
X_protein_id_file="../Pan_genome_data/X_protein_id.txt"
# sequence_cleaner_detect(pan_protein_file,X_count_table_file,X_protein_id_file,"X")

GFF_path="../../GFF/"
strain_id="../Pan_genome_data/read_busco.txt"
list_file="../Pan_genome_data/list_file.txt"
gene_base="../Pan_genome_data/gene_base.fasta"
cat_err_file="../Pan_genome_data/cat_err.txt"
# merge_to_one(GFF_path,strain_id,list_file,gene_base,cat_err_file)

pan_gene="../Pan_genome_data/pan_gene.fasta"
gene_protein_mapping_table_file_name="../Pan_genome_data/gene_protein_mapping_table.txt"
# extract_gene(pan_protein_file,gene_base,pan_gene,gene_protein_mapping_table_file_name)
# don't froget remove

contig_path="../../contig/"
blast_specise_out_path="../Pan_genome_data/c_blast_present_contig/"
pav_excel="../Pan_genome_data/pav.xlsx"
# blast(contig_path,strain_id,pan_gene,blast_specise_out_path,pav_excel)
# blast('../contig_base/test_2/',strain_id,pan_gene,blast_specise_out_path,"../Pan_genome_data/pav_70-15.xlsx")

R_result=Path("../Pan_genome_data/R_result/")
# R_result.mkdir()
pan_categories_gene_path_name="../Pan_genome_data/pan_categories_gene"
pan_categories_protein_path_name="../Pan_genome_data/pan_categories_protein/"
# pan_categories_gene_path_name.mkdir()
# pan_categories_protein_path_name.mkdir()
pav_0_1_tsv="../Pan_genome_data/pav_0_1.tsv"
# R_run_pav_instance_maker(
#     pav_excel,
#     gene_protein_mapping_table_file_name,
#     str(R_result)+'/',
#     pav_0_1_tsv
#     )
# # heatmap()
# set_minus_prepare("70-15")
# draw_stack()
# set_minus(str(R_result)+"/")
# drawer_curve(str(R_result)+"/")


pav_with_cluster_out_file_name="../Pan_genome_data/pav_with_cluster.tsv"
clade_set=["1","2","3","4"]
color_set="Set2"
cluster_clade_file_name="../Pan_genome_data/R_result/strain_clade_category_ID.txt"   

# cluster(pav_excel,clade_set,color_set,pav_with_cluster_out_file_name,cluster_clade_file_name,str(R_result)+"/")
# write_cluster_result()
# set_minus_with_cluster(str(R_result)+"/",cluster_clade_file_name,clade_set,color_set)

orthofinder_assianed_tsv_file_name="../Pan_genome_data/Results_May22_1/Orthogroups/Orthogroups.tsv"
orthofinder_unassianed_tsv_file_name="../Pan_genome_data/Results_May22_1/Orthogroups/Orthogroups_UnassignedGenes.tsv"
set_minus_orthofinder_result=directory_creater("../Pan_genome_data/set_minus_orthofinder_result")
pav_orthofinder_file_name="../Pan_genome_data/set_minus_orthofinder_result/pav_orthofinder.xlsx"
# R_set_minus_cut_orthofinder(orthofinder_assianed_tsv_file_name,orthofinder_unassianed_tsv_file_name,pan_id_file,str(set_minus_orthofinder_result)+"/")
# set_minus(str(set_minus_orthofinder_result)+"/")
# drawer_curve(str(set_minus_orthofinder_result)+"/")

MGG_70_15="../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta"
Augustus_70_15="../../GFF/70-15_gene.fasta"
MGG_Augustus_70_15_dir=directory_creater("../Pan_genome_data/70-15_MGG_Augustus")
MGG_Augustus_70_15_blast_db_dir=directory_creater(MGG_Augustus_70_15_dir/"blast_db")
blast_out_asn_file=MGG_Augustus_70_15_dir/"blast_out.asn"
blast_out_xml_file=MGG_Augustus_70_15_dir/"blast_out.xml"
blast_out_txt_file=MGG_Augustus_70_15_dir/"blast_out.txt"
MGG_unpresent_Augustus=MGG_Augustus_70_15_dir/"MGG_unpresent_Augustus_list.txt"
MGG_Augustus_70_15_blast_db=MGG_Augustus_70_15_blast_db_dir/"Augustus_70_15"
pav_MGG_unpresent_Augustus_file_name=MGG_Augustus_70_15_dir/"pav_MGG_unpresent_Augustus.xlsx"
MGG_unpresent_Augustus_unassianed_list_file_name=MGG_Augustus_70_15_dir/"MGG_unpresent_Augustus_unassianed_list.txt"
# blastdb(Augustus_70_15,MGG_Augustus_70_15_blast_db)
# blast(MGG_Augustus_70_15_blast_db,MGG_70_15,blast_out_asn_file,blast_out_xml_file,blast_out_txt_file)
# read_blast_result(blast_out_xml_file,MGG_70_15,str(MGG_unpresent_Augustus))
# MGG_unpresent_Augustus_locate_in_pav_orthofinder(
#     str(MGG_unpresent_Augustus),
#     pav_orthofinder_file_name,
#     gene_protein_mapping_table_file_name,
#     str(pav_MGG_unpresent_Augustus_file_name),
#     orthofinder_unassianed_tsv_file_name,
#     str(MGG_unpresent_Augustus_unassianed_list_file_name)
#     )
pav_result_dir=directory_creater("../Pan_genome_data/pav_result")
pav_orthofinder_removed_1574=pav_result_dir/"pav_orthofinder_removed_1574.xlsx"
# remove_MGG_unpresent_Augustus(pav_orthofinder_file_name,MGG_unpresent_Augustus_unassianed_list_file_name,pav_orthofinder_removed_1574)

MGG_unpresent_Augustus_set_minus_dir=directory_creater(MGG_Augustus_70_15_dir/"set_minus")
# set_minus_cut("70-15",pav_orthofinder_removed_1574,str(MGG_unpresent_Augustus_set_minus_dir)+'/')
# set_minus(str(MGG_unpresent_Augustus_set_minus_dir)+'/')
# drawer_curve(str(MGG_unpresent_Augustus_set_minus_dir)+'/')

annotation_interproscan_go_pfam_file_name="../Pan_genome_data/annotation/raw_input/pan_protein.fasta.gff3"
# annotation_go_pfam(
#     annotation_interproscan_go_pfam_file_name,
#     pan_categories_protein_path_name,
#     str(R_result)+"/",
#     str(R_result)+"/go_annotation.txt",
#     str(R_result)+"/pfam_annotation.txt"
#     )
kegg_input_path_name="../Pan_genome_data/annotation/raw_input/"
# annotation_kegg(
#     kegg_input_path_name,
#     pan_categories_protein_path_name,
#     str(R_result)+"/",
#     str(R_result)+"/kegg_annotation.txt"
#     )
egg_cog_input_file_name="../Pan_genome_data/annotation/raw_input/query_seqs.fa.emapper.annotations.xlsx"
# annotation_eggnog_cog(
#     egg_cog_input_file_name,
#     pan_categories_protein_path_name,
#     str(R_result)+"/"
#     )
#分泌蛋白
Protcomp_file_name="../Pan_genome_data/annotation/secreted_proteins/input/Protcomp_result.txt"
signalp_file_name="../Pan_genome_data/annotation/secreted_proteins/input/signalp_pan_protein.gff3"
tmhmm_file_name="../Pan_genome_data/annotation/secreted_proteins/input/tmhmm_out.txt"
new_secreted_protein_file_name="../Pan_genome_data/annotation/secreted_proteins/new_tpye_secreted_protein.txt"
classic_secreted_protein_file_name="../Pan_genome_data/annotation/secreted_proteins/classic_secreted_protein.txt"
known_secreted_proteins=["MGG_18041", "MGG_13283", "MGG_03685", "MGG_04301", "MGG_12655"]
# secreted_protein(Protcomp_file_name,tmhmm_file_name,signalp_file_name,new_secreted_protein_file_name,classic_secreted_protein_file_name)
# drawer_secreted_protein(
#     pav_0_1_tsv,
#     new_secreted_protein_file_name,
#     classic_secreted_protein_file_name,
#     gene_protein_mapping_table_file_name,
#     known_secreted_proteins,
#     str(R_result)+"/"
#     )
# ortholog_path_name="../Pan_genome_data/Results_May22_1/Orthologues/Orthologues_70-15_protein/"
ortholog_path_name="../Pan_genome_data/ortholog_test/"
joined_df_file_name="../Pan_genome_data/ortholog/joined_df.tsv"
ortholog_blast_path_name="../Pan_genome_data/ortholog/ortholog_blast/"
# read_file(ortholog_path_name,joined_df_file_name)
contig_length_json_dir=directory_creater("../Pan_genome_data/contig_length_json")
# contig_length2json(contig_path,strain_id,contig_length_json_dir)
# extract_ortholog_gene(gene_base,"../Pan_genome_data/ortholog/test.txt",strain_id,ortholog_blast_path_name)
# extract_ortholog_gene(gene_base,joined_df_file_name,strain_id,ortholog_blast_path_name)
snp_gatk_result_path=directory_creater("../Pan_genome_data/call_snp_long")
# call_snp(strain_id,contig_path,snp_result_path)
# GenomicsDBImport(snp_result_path)
# CombineGVCFs(snp_result_path)
# GenotypeGVCFs(snp_result_path)

snp_mummer_result_path=directory_creater("../Pan_genome_data/call_snp_mummer")
# call_snp_mummer(strain_id,contig_path,snp_mummer_result_path)

nucleotide_diversity_dir_path=directory_creater("../Pan_genome_data/nucleotide_diversity")
# calculate_nucleotide_diversity(pan_categories_gene_path_name,nucleotide_diversity_dir_path)

# calculate_dnds(
#     pan_gene,
#     pan_protein_file,
#     directory_creater("../Pan_genome_data/calculate_dnds")
#     )
prepare_for_ParaAT_out_dir=directory_creater("../Pan_genome_data/prepare_for_ParaAT")
# prepare_for_ParaAT(joined_df_file_name,strain_id,prepare_for_ParaAT_out_dir)
# paraAT(prepare_for_ParaAT_out_dir)
parse_ParaAT_result(prepare_for_ParaAT_out_dir)