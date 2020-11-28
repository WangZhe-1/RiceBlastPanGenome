'''
@Author: your name
@Date: 2020-07-18 17:20:03
@LastEditTime: 2020-07-18 23:41:57
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/nucleotide_diversity.py
'''
import gffutils
from pathlib import Path
from gffutils.helpers import asinterval
import pybedtools
from Directory_creater import directory_creater
import subprocess
from pybedtools.featurefuncs import gff2bed

def generate_interval(category_fl):
    for strain_id_raw in category_fl:
        strain_id=strain_id_raw.strip('\n')
        if strain_id[0:3]!="MGG": continue
        yield gff2bed(asinterval(MGG_db[strain_id]),name_field=2)
def vcftool(in_file,out_file,stdout_file):
    call_vcftool=subprocess.Popen(
        [
            "vcftools",
            "--vcf",
            "../Pan_genome_data/call_snp_mummer/merge.vcf",
            "--window-pi",
            "--bed",
            str(in_file),
            "--out",
            str(out_file)
        ],
        stdout=stdout_file.open('w'),
        stderr=subprocess.STDOUT
    )
    call_vcftool.wait()
def calculate_nucleotide_diversity(gene_category_dir_name,output_dir):
    '''
    input 1: gene category dir path
    output 1: nucleotide_diversity_dir
    '''
    global MGG_db
    MGG_db=gffutils.FeatureDB("../Pan_genome_data/ortholog/gffutils_db/MGGdb.db")
    gene_category_dir_path=Path(gene_category_dir_name)
    category_bed_dir_path=directory_creater(output_dir/"category_bed")
    pi_result_dir_path=directory_creater(output_dir/"pi_result")
    pi_std_out_err_dir_path=directory_creater(output_dir/"pi_std_out_err")
    for category_file in gene_category_dir_path.iterdir():
        bed_file_path=category_bed_dir_path/(category_file.stem+".bed")
        if bed_file_path.is_file() is False:
            pybedtools.BedTool(generate_interval(category_file.open())).saveas(bed_file_path)
        vcftool(
            bed_file_path,
            pi_result_dir_path/category_file.stem,
            pi_std_out_err_dir_path/(category_file.stem+"_out_err.txt")
            )

    

