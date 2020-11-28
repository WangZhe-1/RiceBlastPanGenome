'''
@Author: your name
@Date: 2020-08-02 16:59:47
LastEditTime: 2020-11-11 09:56:19
LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/general_2.py
'''
from alignment_extract_orthologs_mafft_slop_mgg_key import extract_ortholog_gene
# from mafft_gene import mafft_run
# from calculate_dnds_2 import prepare_for_ParaAT, run_paraAT
# from alignment_extract_orthologs_2 import read_file
# from clade_different_secreted_protein import clade_different_secreted_protein
# from enrichment_clade_different_eggnog_cog import clade_different_eggnog_cog
# from enrichment_clade_different_kegg import clade_different_kegg
# from enrichment_clade_different_interproscan_go_pfam import clade_different_go_pfam
import rpy2.robjects as robjects
# from fgenesh_gff_to_targer_dir import fgenesh_gff_to_targer_dir
# from cluster_orthofinder_5clade_blastn import cluster_blastn, write_cluster_result_blastn
# from pav_R_part_2_blastn import R_run_pav_instance_maker_blastn, draw_stack_blastn, write_class_blastn
# from pav_contig_present_for_pan_data_2_blast_only import pav_contig_present_blast_gth_main
# from annotation_secreted_proteins_2_rice_division import drawer_secreted_protein_rice_division #用于纯水稻宿主完全划分，n是划分份数，s是前后多少份详细表示
# from alignment_extract_orthologs import extract_gene
# from extract_nucleic_2 import contig_gff3_gffutils_db, extract_gene_gff, extract_mRNA, merge_to_one  # 使用于新版的
# from annotation_interproscan_go_pfam_2 import annotation_go_pfam #annotation_interproscan_go_pfam_2适用于没有no_present的情况
# from annotation_kegg_id_2 import annotation_kegg #annotation_kegg_id_2适用于没有no_present的情况
# from annotation_eggnog_cog_2 import annotation_eggnog_cog #annotation_eggnog_cog_2适用于没有no_present的情况
# from cluster_orthofinder_5clade import cluster_rice
# from cluster_orthofinder_5clade import pca
# from filter_too_small_sequence import filter_blast
# from filter_too_small_sequence import filter_pan_id
# from pav_R_part_2 import R_run_pav_instance_maker, draw_stack, write_class
# from upset import make_combination_matrix_all, make_combination_matrix_rice
# from set_minus_with_cluster import set_minus_with_cluster
# from cluster_orthofinder_5clade import cluster, write_cluster_result
from pathlib import Path
# from set_minus import drawer_curve, set_minus
# from set_minus_orthofinder import R_set_minus_cut_orthofinder
# from copy_contig import Copy_contig
from Directory_creater import directory_creater
# from annotation_secreted_proteins_2 import drawer_secreted_protein, secreted_protein # annotation_secreted_proteins_2是没有mRNA_protein_mapping_table_file_name的
# from fgenesh_result_to_fasta_pro import fgenesh_result_to_fasta
# from CD_HIT import cd_hit
# from fgenesh_GFF import writen_name
# from extract_sequence_from_Orthogroups import extract_sequence_from_Orthogroups
'''
copy
'''
contig_path=Path("../../contig/")
MGG_70_15_contig="../../70-15_refference_genome/70-15_supercontigs.fasta"
ina168_contig="../../contig/ina168.fasta"
general_out=directory_creater("../Pan_genome_data_2/")
copy_std_out_err=directory_creater(general_out/"copy_out_err")
contig_156_path=directory_creater(general_out/"156_contig")
# Copy(contig_path,MGG_70_15_contig,ina168_contig,contig_156_path,copy_std_out_err)

'''
phase2：分析预测基因结果
70-15需要跳过，从molquest结果里面删掉，已删掉
ina168需要单独处理，已处理
手动把70-15放入protein_base和mRNA_base
FR13,gene,protein, 手动全部替换，把genome换成FR13，文件名也换成FR13

唯三需要注意的：
70-15
ina168
FR13
'''
mRNA_base_path=directory_creater(general_out/"mRNA_base")
protein_base_path=directory_creater(general_out/"protein_base")
# fgenesh_result_to_fasta(
#     Path("/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000001/"),
#     mRNA_base_path,
#     protein_base_path
# )
'''
提取gff文件
'''
gff_cp_stdout_err_dir_path=directory_creater(general_out/"gff_cp_stdout_err")
gff_path=directory_creater(general_out/"contig_gff")
# fgenesh_gff_to_targer_dir(
#     Path("/mnt/d/windows_program/Molquest2_work_path/MolQuest2/projects/default/tasks/000002/"),
#     gff_path,
#     gff_cp_stdout_err_dir_path
# )

'''
去重
'''
mRNA_base_no_duplicate=directory_creater(general_out/"mRNA_base_no_duplicate")
protein_base_no_duplicate=directory_creater(general_out/"protein_base_no_duplicate")
# cd_hit(
#     mRNA_base_path,
#     protein_base_path,
#     mRNA_base_no_duplicate,
#     protein_base_no_duplicate,
#     directory_creater(general_out/"cd_hit_std_out"),
#     directory_creater(general_out/"cd_hit_std_err"),
#     )

'''
验证contig和fasta文件名是否一致，这是菌株名，后面要相互调用，需保证一致
把contig里面的Guy11换成GUY11
'''
fasta_stem_name_list=general_out/"fasta_stem_name_list.txt"
# writen_name(protein_base_path,fasta_stem_name_list)
contig_stem_name_list=general_out/"contig_stem_name_list.txt"
# writen_name(contig_156_path,contig_stem_name_list)

'''
从orthofinder中提取出pan protein
提取后，泛基因组总数由35405降到35399，是由于在unassigned 基因中有6个含有X的基因，被删掉了。
'''
orthogroup_result_path=general_out/"Results_Aug04"
orthogroup_sequence_path=orthogroup_result_path/"Orthogroup_Sequences"
pan_dir_path=directory_creater(general_out/"pan")
pan_protein_file=pan_dir_path/"pan_protein.fasta"
pan_id_file=pan_dir_path/"pan_id.txt"
# extract_sequence_from_Orthogroups(orthogroup_sequence_path,pan_protein_file,pan_id_file)

'''
过滤掉长度小于20的
'''
self_blast_result=pan_dir_path/"blast.xml"
filtered_pan_protein_file_path=pan_dir_path/"pan_protein_no_shorter_20.fasta"
filtered_length_table=pan_dir_path/"length_table.txt"
filtered_pan_id_file_name=pan_dir_path/"filtered_pan_id.txt"
# filter_blast(self_blast_result,pan_protein_file,filtered_pan_protein_file_path,filtered_length_table)
# filter_pan_id(str(pan_id_file),str(filtered_length_table),str(filtered_pan_id_file_name))
'''
提取mRNA，把蛋白转为mRNA
'''
mRNA_base_file_name=general_out/"mRNA_base.fasta"
cat_err_file=general_out/"cat_err.txt"
pan_mRNA_file_name=pan_dir_path/"pan_mRNA_no_shorter_20.fasta"
# merge_to_one(mRNA_base_no_duplicate,general_out/"merge_list.txt",mRNA_base_file_name,cat_err_file)
# extract_mRNA(filtered_pan_protein_file_path,str(mRNA_base_file_name),pan_mRNA_file_name)
# extract_gene(pan_protein_file,str(mRNA_base_file_name),pan_mRNA_file_name)
'''
提取gene

使用gff数据库，唯三需要注意的
70-15不在数据库里面
db文件名是Ina168的需要重命名为ina168因为基因ID存的是ina168，split分解基因ID找gff文件用的是ina168，找不到Ina168
split分解基因ID得到strain_id是WD-3-1需要特别处理
if strain_id=="WD-3-1":
                    strain_id=protein_sequence.id[0:8]
                    gff_protein_id=protein_sequence.id[9:]
'''
contig_gff3_gffutils_db_dir_path=directory_creater(general_out/"contig_gff3_gffutils_db")
# contig_gff3_gffutils_db(gff_path,contig_gff3_gffutils_db_dir_path)
pan_gene_file_path=pan_dir_path/"pan_gene_no_shorter_20.fasta"
# extract_gene_gff(filtered_pan_protein_file_path,contig_gff3_gffutils_db_dir_path,contig_156_path,pan_gene_file_path)

'''
核心下降曲线
包含生成pav 矩阵，excel，pav_orthofinder_file_name，pav_orthofinder.xlsx
'''
orthofinder_assianed_tsv_file_name=orthogroup_result_path/"Orthogroups/Orthogroups.tsv"
orthofinder_unassianed_tsv_file_name=orthogroup_result_path/"Orthogroups/Orthogroups_UnassignedGenes.tsv"

'''
set_minus_orthofinder_result是使用未过滤的pan_protein的结果，
set_minus_orthofinder_result_2是使用过滤掉<20后得到的结果。
'''
set_minus_orthofinder_result=directory_creater(general_out/"set_minus_orthofinder_result_2") 
pav_orthofinder_file_name=set_minus_orthofinder_result/"pav_orthofinder.xlsx"

# R_set_minus_cut_orthofinder(
#     str(orthofinder_assianed_tsv_file_name),
#     str(orthofinder_unassianed_tsv_file_name),
#     str(filtered_pan_id_file_name),
#     str(set_minus_orthofinder_result)+"/"
#     )
# set_minus(str(set_minus_orthofinder_result)+"/")
# drawer_curve(str(set_minus_orthofinder_result)+"/")

'''
使用orthofinder结果的聚类
'''
# color_num=5
color_num=4
pav_with_cluster_out_file_name=set_minus_orthofinder_result/("pav_with_cluster_"+str(color_num)+"_binary.tsv")
clade_set=[str(i) for i in range(1,color_num+1)]
color_set="Set2"
# cluster_clade_file_name=set_minus_orthofinder_result/("strain_clade_category_ID_"+str(color_num)+"_binary.txt")

# cluster(
#     str(pav_orthofinder_file_name),
#     clade_set,
#     color_set,
#     color_num,
#     str(pav_with_cluster_out_file_name),
#     str(cluster_clade_file_name),
#     str(set_minus_orthofinder_result)+"/"
#     )

# write_cluster_result()
# set_minus_with_cluster(
#     str(set_minus_orthofinder_result)+"/",
#     str(cluster_clade_file_name),
#     clade_set,             #决定号码对应的颜色
#     color_set,
#     color_num,
#     clade_set   #决定在核心下降曲线中的顺序
#     )
'''
使用orthofinder结果的pca
'''
# pca_png_result_file_path=set_minus_orthofinder_result/("pca_"+str(color_num)+".png")
# pca_eig_png_file_name=set_minus_orthofinder_result/("pca_eig_"+str(color_num)+".png")
# pca(str(pav_orthofinder_file_name),str(cluster_clade_file_name),str(pca_png_result_file_path),str(pca_eig_png_file_name))

'''
用于rice内部的聚类
'''
color_num_rice=3
clade_set_rice=[str(i) for i in range(1,color_num_rice+1)]
rice_cluster_clade_file_name=set_minus_orthofinder_result/("rice_strain_clade_category_ID_"+str(color_num_rice)+"_binary.txt")
rice_pav_with_cluster_out_file_name=set_minus_orthofinder_result/("rice_pav_with_cluster_"+str(color_num_rice)+"_binary.tsv")
pav_rice_out_excel_name=set_minus_orthofinder_result/"pav_orthofinder_rice.xlsx"
# cluster_rice(
#     str(pav_with_cluster_out_file_name),
#     clade_set_rice,
#     color_set,
#     color_num_rice,
#     str(set_minus_orthofinder_result)+"/",
#     str(rice_pav_with_cluster_out_file_name),
#     str(rice_cluster_clade_file_name),
#     str(pav_rice_out_excel_name)
#     )
'''
用于Augustus_based_orthofinder的聚类
'''
# cluster(
#     "../Pan_genome_data/set_minus_orthofinder_result/pav_orthofinder.xlsx",
#     clade_set,
#     color_set,
#     color_num,
#     str(directory_creater(general_out/"old")/("pav_with_cluster_"+str(color_num)+".tsv")),
#     str(directory_creater(general_out/"old")/("strain_clade_category_ID_"+str(color_num)+".txt")),
#     str(directory_creater(general_out/"old"))+"/"
#     )
# set_minus_with_cluster(
#     "../Pan_genome_data/set_minus_orthofinder_result/",
#     str(directory_creater(general_out/"old")/("strain_clade_category_ID_"+str(color_num)+".txt")),
#     clade_set,
#     color_set,
#     color_num,
#     ["1","2","3","4"]
#     )

'''
基于blastn的pav矩阵
'''
blast_intermediate_out=directory_creater(general_out/"blast_intermediate_out_gene")
pav_blastn=set_minus_orthofinder_result/"pav_blastn_gene.xlsx"
# pav_contig_present_blast_gth_main(contig_156_path,pan_gene_file_path,blast_intermediate_out,pav_blastn)
'''
划分基因类,
基因类在pan_categories_protein中
'''
# pan_categories_protein_dir_path=directory_creater(set_minus_orthofinder_result/"pan_categories_protein")

# R_run_pav_instance_maker(
#     str(pav_orthofinder_file_name),
#     str(set_minus_orthofinder_result)+'/'
#     )
# write_class()
# # heatmap()
# draw_stack()
'''
划分基因类—_mRNA
使用pav_blastn的基因类在pan_categories_mRNA_blastn_dir_path
'''
pan_categories_mRNA_blastn_dir_path=directory_creater(set_minus_orthofinder_result/"pan_categories_mRNA_blastn")
# R_run_pav_instance_maker_blastn(
#     str(pav_blastn),
#     str(set_minus_orthofinder_result)+'/'
# )
# write_class_blastn()
# draw_stack_blastn()

'''
划分基因类—_gene
使用pav_blastn的基因类在pan_categories_gene_blastn_dir_path

属于no present 基因：
MGG_00088，里面有101个n，超过设定的50，被判断为1，计入未出现的基因
考虑可能是因为蛋白没有X，而基因有N
'''
# pan_categories_gene_blastn_dir_path=directory_creater(set_minus_orthofinder_result/"pan_categories_gene_blastn")
# R_run_pav_instance_maker_blastn(
#     str(pav_blastn),
#     str(set_minus_orthofinder_result)+'/'
# )
# write_class_blastn()
# draw_stack_blastn()

'''
基于blastn,pav结果的聚类
'''
# color_num=5
color_num=4
pav_with_cluster_out_file_name_blastn=set_minus_orthofinder_result/("pav_with_cluster_"+str(color_num)+"_blastn.tsv")
cluster_clade_file_name_blastn=set_minus_orthofinder_result/("strain_clade_category_ID_"+str(color_num)+"_blastn.txt")

# cluster_blastn(
#     str(pav_blastn),
#     clade_set,
#     color_set,
#     color_num,
#     str(pav_with_cluster_out_file_name_blastn),
#     str(cluster_clade_file_name_blastn),
#     str(set_minus_orthofinder_result)+"/"
#     )

# write_cluster_result_blastn()
# set_minus_with_cluster(
#     str(set_minus_orthofinder_result)+"/",
#     str(cluster_clade_file_name_blastn),
#     clade_set,             #决定号码对应的颜色
#     color_set,
#     color_num,
#     clade_set,   #决定在核心下降曲线中的顺序
#     "blastn"
#     )
'''
基于pan_gene的pca,没有做
'''
# pca_png_result_file_path=set_minus_orthofinder_result/("pca_"+str(color_num)+".png")
# pca_eig_png_file_name=set_minus_orthofinder_result/("pca_eig_"+str(color_num)+".png")
# pca_blastn(str(pav_orthofinder_file_name),str(cluster_clade_file_name),str(pca_png_result_file_path),str(pca_eig_png_file_name))

'''
upset：各个clade之间基因的差集
提供的结果是交集，交集提供的是核心，即哪些基因在clade 1 是核心，但是在clade 2 就不是了。
要想变为并集，把这里的&改为|即可。
“result=result&pav_logical[[clade]]”
并集提供的是clade整体的情况，哪些基因属于clade 1，但是不在clade 2中。
upset中包含从pav_blastn中得到纯水稻（rice）宿主的R代码
'''

zhong_zhe_merge_df_all=robjects.DataFrame({
      "zhong":robjects.StrVector(["clade_1","clade_2","clade_3","other"]),
      "zhe":robjects.StrVector(["1","4","2","3"])
    })
combination_matrix_distinct_dir_path=directory_creater(general_out/"combination_matrix_distinct_core")
combination_matrix_distinct_all_dir_path=directory_creater(combination_matrix_distinct_dir_path/"all")
clade_combination=robjects.StrVector(["clade_1&clade_2","clade_2&clade_3","clade_1&clade_3","clade_1","clade_2","clade_3","clade_1&clade_2&clade_3",""])
# make_combination_matrix_all(
#     str(pav_blastn),
#     str(cluster_clade_file_name_blastn),
#     zhong_zhe_merge_df_all,
#     "distinct",
#     clade_combination,
#     str(combination_matrix_distinct_all_dir_path)+'/'
#     )

zhong_zhe_merge_df_rice=robjects.DataFrame({
      "zhong":robjects.StrVector(["clade_1","clade_2","clade_3"]),
      "zhe":robjects.StrVector(["1","4","2"])
    })
combination_matrix_distinct_rice_dir_path=directory_creater(combination_matrix_distinct_dir_path/"rice")
rice_pan_gene_file_name=pan_dir_path/"rice_pan_gene.txt"
rice_pav_file_name=set_minus_orthofinder_result/"rice_pav.tsv"
# make_combination_matrix_rice(
#     str(pav_blastn),
#     str(cluster_clade_file_name_blastn),
#     zhong_zhe_merge_df_rice,
#     "distinct",
#     clade_combination,
#     str(combination_matrix_distinct_rice_dir_path)+"/",
#     str(rice_pan_gene_file_name),
#     str(rice_pav_file_name)
# )

'''
go、kegg、pfam的注释
'''
annotation_result_path=general_out/"annotation"
annotation_interproscan_go_pfam_file_name="../Pan_genome_data_2/annotation/pan_protein_no_shorter_20.fasta.gff3"
# annotation_go_pfam(
#     annotation_interproscan_go_pfam_file_name,
#     pan_categories_protein_dir_path,
#     str(annotation_result_path)+"/",
#     str(annotation_result_path)+"/go_annotation.txt",
#     str(annotation_result_path)+"/pfam_annotation.txt"
#     )
kegg_input_path_name="../Pan_genome_data_2/annotation/"
# annotation_kegg(
#     kegg_input_path_name,
#     str(pan_categories_protein_dir_path),
#     str(annotation_result_path)+"/",
#     str(annotation_result_path)+"/kegg_annotation.txt"
#     )
egg_cog_input_file_name="../Pan_genome_data_2/annotation/query_seqs.fa.emapper.annotations.xlsx"
# annotation_eggnog_cog(
#     egg_cog_input_file_name,
#     pan_categories_protein_dir_path,
#     str(annotation_result_path)+"/"
#     )

'''
针对不同clade差异的富集分析
当心：
universe，注释库，和gene_taxonomy的ID是否相符
'''
clade_different_enrichment_path=directory_creater(general_out/"clade_different_enrichment")
# clade_different_go_pfam(
#     annotation_interproscan_go_pfam_file_name,
#     combination_matrix_distinct_rice_dir_path,
#     clade_combination,
#     rice_pan_gene_file_name,
#     str(clade_different_enrichment_path)+"/",
#     str(clade_different_enrichment_path)+"/go_annotation.txt",
#     str(clade_different_enrichment_path)+"/pfam_annotation.txt"
#     )
# clade_different_kegg(
#     kegg_input_path_name,
#     str(combination_matrix_distinct_rice_dir_path),
#     clade_combination,
#     rice_pan_gene_file_name,
#     str(clade_different_enrichment_path)+"/",
#     str(clade_different_enrichment_path)+"/kegg_annotation.txt"
#     )
# clade_different_eggnog_cog(
#     egg_cog_input_file_name,
#     str(combination_matrix_distinct_rice_dir_path),
#     clade_combination,
#     rice_pan_gene_file_name,
#     str(clade_different_enrichment_path)+"/"
#     )

'''
分泌蛋白
'''

Protcomp_file_name="../Pan_genome_data_2/annotation/protcomp.txt"
signalp_file_name="../Pan_genome_data_2/annotation/pan_protein_no_shorter_20.gff3"
tmhmm_file_name="../Pan_genome_data_2/annotation/tmhmm_result.txt"
new_secreted_protein_file_name="../Pan_genome_data_2/annotation/new_secreted_protein.txt"
classic_secreted_protein_file_name="../Pan_genome_data_2/annotation/classic_secreted_protein.txt"
known_secreted_proteins=["MGG_18041", "MGG_13283", "MGG_03685", "MGG_04301", "MGG_12655"]
# secreted_protein(Protcomp_file_name,tmhmm_file_name,signalp_file_name,new_secreted_protein_file_name,classic_secreted_protein_file_name)

# drawer_secreted_protein(
#     str(pav_orthofinder_file_name),
#     new_secreted_protein_file_name,
#     classic_secreted_protein_file_name,
#     known_secreted_proteins,
#     str(annotation_result_path)+"/"
#     )
'''
用于纯水稻宿主完全划分，n是划分份数，s是前后多少份详细表示
'''
# drawer_secreted_protein_rice_division(
#     str(rice_pav_with_cluster_out_file_name),
#     str(new_secreted_protein_file_name),
#     str(classic_secreted_protein_file_name),
#     known_secreted_proteins,
#     str(annotation_result_path)+"/",
#     40,
#     1
# )
# secreted_protein(
#     "../../../../陈瑞齐/secrete_protein/result.txt",
#     "../../../../陈瑞齐/secrete_protein/NP_ill_pro_tmhmm.txt",
#     "../../../../陈瑞齐/secrete_protein/NP_ill_pro_signalp.gff3",
#     "../../../../陈瑞齐/secrete_protein/new",
#     "../../../../陈瑞齐/secrete_protein/classic"
#     )
# secreted_protein(
#     "../../../../陈瑞齐/PB_secrete/molquest_result.txt",
#     "../../../../陈瑞齐/PB_secrete/NP_ill_pro_tmhmm.txt",
#     "../../../../陈瑞齐/PB_secrete/NP_ill_pro_signalp.gff3",
#     "../../../../陈瑞齐/PB_secrete/new",
#     "../../../../陈瑞齐/PB_secrete/classic"
#     )


'''
是否有分泌蛋白富集
'''
# clade_different_secreted_protein(
#     combination_matrix_distinct_rice_dir_path,
#     clade_combination,
#     rice_pan_gene_file_name,
#     "distinct",
#     classic_secreted_protein_file_name,
#     new_secreted_protein_file_name,
#     str(clade_different_enrichment_path)+"/"
# )
ortho_dir_path=directory_creater(general_out/"ortho")
joined_df_name_path=ortho_dir_path/"joined_df.tsv"
# read_file(orthogroup_result_path/"Orthologues"/"Orthologues_70-15",str(joined_df_name_path))

protein_base_file_name=general_out/"protein_base.fasta"
protein_cat_err_file=general_out/"protein_cat_err.txt"
# merge_to_one(protein_base_no_duplicate,general_out/"protein_merge_list.txt",protein_base_file_name,protein_cat_err_file)
# prepare_for_ParaAT(str(joined_df_name_path),mRNA_base_file_name,protein_base_file_name,ortho_dir_path)
# run_paraAT(ortho_dir_path)
# parse_ParaAT_result(ortho_dir_path)

'''
Single_Copy_Orthologue 的联配
mafft
'''
mafft_result_dir_path=directory_creater(general_out/"mafft_result")
mafft_stderr_dir_path=directory_creater(general_out/"mafft_stderr")
Single_Copy_Orthologue_gene_dir_path=directory_creater(general_out/"Single_Copy_Orthologue_gene")
# mafft_run(
#     orthogroup_result_path/"Single_Copy_Orthologue_Sequences",
#     Single_Copy_Orthologue_gene_dir_path,
#     contig_gff3_gffutils_db_dir_path,
#     contig_156_path,
#     mafft_result_dir_path,
#     mafft_stderr_dir_path
# )

'''
mafft那些出现在低于156个菌株的
extract_ortholog_gene，从orthofinder大表格里面提取基因，
先做单拷贝基因自己的blast，确定阈值，
再把有多拷贝的加入进去（参考blast指标选择相似的）
有上下游1k，alignment_extract_orthologs_mafft_slop.py
仅基因，alignment_extract_orthologs_mafft.py


应该以每一个MGG基因作为key，如果MGG位于青蛙的位置，那么会出现一个MGG对应多个其他基因的情况
如果MGG位于人的位置，那么这个MGG有唯一的对应基因
如果MGG基因位于人和狗的位置，那么会出现多个MGG对应一个外类基因，此时，单个的MGG必然会有其单个对应的唯一对应的外类基因，因此需要使用单个的MGG作为key。

只考虑了一个hsp的最大值，没有考虑是否有断
检查刷下来的低于0.95的有没有问题,检查低的是为啥
得到的是不带MGG的基因总数

'''
ortholog_blast_path_name=directory_creater(general_out/"ortholog_blast_mgg_key")
# joined_df_name_path=general_out/"ortholog_blast"/"linshi.txt"
# 需要删掉ortholog_blast中对于已经做过文件读取的注释
# 写
extract_ortholog_gene(
    contig_gff3_gffutils_db_dir_path,
    contig_156_path,
    joined_df_name_path,
    ortholog_blast_path_name
    # orthogroup_sequence_path,
    # pav_orthofinder_file_name,
    # 150,
    # ortholog_blast_path_name
)
'''
找到公共的阈值
每一个函数环节设置检验办法
'''