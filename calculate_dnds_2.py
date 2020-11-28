'''
@Author: your name
@Date: 2020-07-20 16:48:36
LastEditTime: 2020-10-18 20:29:02
LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/calculate_dnds.py
'''
from joblib import Parallel,delayed
import re
from tempfile import NamedTemporaryFile
import subprocess
from extract_strain_id import extract_strain_id
from Bio import SeqIO
from rpy2.robjects.packages import importr
from Directory_creater import directory_creater
from pathlib import Path

# def uniformize(in_file,out_file,fun):
#     if out_file.is_file() is False:
#         with open(in_file) as in_fl:
#             with open(out_file,'w') as out_fl:
#                 for seq in SeqIO.parse(in_fl,"fasta"):
#                     seq.id=fun(seq.id)
#                     SeqIO.write(seq,out_fl,"fasta")

# def calculate_dnds(pan_gene,pan_protein,out_dir):
#     uniformize_pan_gene_file_path=out_dir/"uniformize_pan_gene.fasta"
#     uniformize_pan_protein_file_path=out_dir/"uniformize_pan_protein.fasta"
    
#     uniformize(pan_gene,uniformize_pan_gene_file_path,get_id_gene)
#     uniformize(pan_protein,uniformize_pan_protein_file_path,get_id_protein)
def extract_gene(uniformize_id):
    gene_sequence=coding_gene_base[uniformize_id]
    gene_sequence.name=""
    gene_sequence.description=""
    SeqIO.write(gene_sequence,gene_fasta_fl,"fasta")

    protein_sequence=protein_base[uniformize_id]
    protein_sequence.name=""
    protein_sequence.description=""
    SeqIO.write(protein_sequence,protein_fasta_fl,"fasta")
def merge_to_one(fasta_list,id_file,filter_string,list_file,base_name,cat_err_file):
    '''
    input 1: list_head
    input 2: strain list
    input 3: filter_string
    output 1: list file
    output 2: base name
    output 3: cat_err_file
    '''
    strain_95_list=extract_strain_id(id_file)
    GFF_path=Path("../../GFF/")
    with open (list_file,'w+') as list_out:
        for species_id in strain_95_list:
            fasta_file_path=GFF_path/(species_id.strip('\n')+filter_string)
            list_out.write('{}\n'.format(fasta_file_path))
            fasta_list.append(str(fasta_file_path.resolve()))
        
    fasta_cat=subprocess.Popen(
        fasta_list,
        stdout=base_name.open('w'),
        stderr=open(cat_err_file,'w+'),
        universal_newlines=True,
    )
    fasta_cat.wait()
def one2two():
    pass
def one_head(df_row):
    for R_cell in df_row:
        if R_cell[0]=='NA' or (str(R_cell[0])=="NA"):
            continue
        body_list=R_cell[0].split()
        if len(body_list)>1:
            one2two()
        else:
            uniformize_id=body_list[0]
            extract_gene(uniformize_id)
            yield uniformize_id
def two_head(p1,p2):
    pass
def prepare_for_ParaAT(joined_df_file_name,coding_gene_base_file_path,protein_base_file_path,out_dir):
    '''
    input 1: joined_df_file_name
    input 2: coding_gene_base_file_path
    input 3: protein_base_file_path
    output 1: out_dir
    '''
    global coding_gene_base,protein_base,gene_fasta_fl,protein_fasta_fl
    gene_dir_path=directory_creater(out_dir/"nucleotide")
    protein_dir_path=directory_creater(out_dir/"aminoacid")
    fasta_base_dir_path=directory_creater(out_dir/"gene_protein_base")
    homolog_dir_path=directory_creater(out_dir/"homolog")
    coding_gene_base=SeqIO.index(
        str(coding_gene_base_file_path),
        "fasta"
        )
    protein_base=SeqIO.index(
        str(protein_base_file_path),
        "fasta"
        )
    base=importr("base")
    utils=importr("utils")
    ortholog_joined_df=utils.read_table(
        joined_df_file_name,
        sep = "\t",
        header = True,
        **{'stringsAsFactors': False},
        **{'check.names': False}
        )
    ortholog_joined_df_sub=ortholog_joined_df.rx(True,-1)
    for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)):
                df_row=ortholog_joined_df_sub.rx(i, True)
                df_row_iter=iter(df_row)
                head_list=next(df_row_iter)[0].split()
                if len(head_list)==1:
                    gene_fasta=gene_dir_path/(head_list[0]+".fasta")
                    protein_fasta=protein_dir_path/(head_list[0]+".fasta")
                    homolog_file_path=homolog_dir_path/(head_list[0]+".txt")
                    if gene_fasta.is_file() is True:continue
                    with gene_fasta.open('w') as gene_fasta_fl:
                        with protein_fasta.open('w') as protein_fasta_fl:
                            with homolog_file_path.open('w') as homolog_fl:
                                extract_gene(head_list[0])
                                homolog_fl.write(head_list[0]+"\t")
                                for homolog_id in one_head(df_row_iter):
                                    homolog_fl.write(homolog_id+"\t")
                                homolog_fl.write("\n")
                else:
                    two_head(head_list,df_row_iter)
def format_paraAT_parameter(in_dir):
    nucleotide_dir_path=in_dir/"nucleotide"
    homolog_dir_path=in_dir/"homolog"
    protein_dir_path=in_dir/"aminoacid"
    out_dir_path=directory_creater(in_dir/"ParaAT_out")
    for gene_id in nucleotide_dir_path.iterdir():
        f=NamedTemporaryFile('w+t',delete=False)
        f.write('1')
        yield ([
            "/mnt/d/zhes_learning_space/software_in_ubuntu/ParaAT2.0/ParaAT.pl",
            "-h",
            str(homolog_dir_path/(gene_id.stem+".txt")),
            "-n",
            str(nucleotide_dir_path/(gene_id.stem+".fasta")),
            "-a",
            str(protein_dir_path/(gene_id.stem+".fasta")),
            "-p",
            f.name,
            "-m",
            "muscle",
            "-f",
            "axt",
            # "-g",
            "-k",
            "-o",
            str(out_dir_path/gene_id.stem)
        ],
        gene_id.stem
        )
    
def Popen_paraAT(parameter,ParaAT_stdout,ParaAT_stderr):
    proc=subprocess.Popen(
        parameter[0],
        stdout=(ParaAT_stdout/(parameter[1]+".txt")).open('w'),
        stderr=(ParaAT_stderr/(parameter[1]+".txt")).open('w')
    )
    proc.wait()
def run_paraAT(in_dir):
    ParaAT_stdout=directory_creater(in_dir/"ParaAT_stdout")
    ParaAT_stderr=directory_creater(in_dir/"ParaAT_stderr")
    Parallel(n_jobs=12)(delayed(Popen_paraAT)(i,ParaAT_stdout,ParaAT_stderr) for i in format_paraAT_parameter(in_dir))
def parse_ParaAT_result(paraAT_result_dir_path):
    '''
    input 1:paraAT_result_dir_path
    output 2: parsed paraAT result:
    like:
    Sequence\tMethod\tKa\tKs\tKa/Ks\n
    '''
    in_path=paraAT_result_dir_path/"ParaAT_out"
    out_file=paraAT_result_dir_path/"parsed_ParaAT_result.tsv"
    with open(out_file,'w') as out_fl:
        with open(paraAT_result_dir_path/"grep_err.txt",'w') as grep_err_fl:
            out_fl.write("Sequence\tMethod\tKa\tKs\tKa/Ks\n")
            for kaks_result_file in in_path.rglob("*.kaks"):
                cut=subprocess.Popen(
                    [
                        "cut",
                        "-f",
                        "1,2,3,4,5",
                        str(kaks_result_file)
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                cut_output, cut_errors=cut.communicate()
                if cut.returncode != 0:
                    print ("cut failed for %s: %s" % (kaks_result_file, cut_errors))
                else:
                    grep=subprocess.Popen(
                        [
                            "grep",
                            "-v",
                            "Sequence"
                        ],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                    grep_out,gerp_err=grep.communicate(input=cut_output)
                    if grep.returncode!=0:
                        grep_err_fl.write(str(kaks_result_file)+'\t'+gerp_err+'\n')
                    else:
                        grep_out_sub=re.sub(".*(MGG_.+T0).*?\\t",lambda m: m.group(1)+"\t",grep_out)
                        out_fl.write(grep_out_sub)




