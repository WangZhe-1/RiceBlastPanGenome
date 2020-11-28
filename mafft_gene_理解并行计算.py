'''
Author: your name
Date: 2020-10-18 11:06:40
LastEditTime: 2020-10-18 18:58:46
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/mafft.py
'''
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet
from Bio import SeqIO
import gffutils
import subprocess
from joblib import Parallel,delayed
def extract_gene_gff_mafft(protein_file_path,contig_path,out_gene_file):
    with out_gene_file.open("w") as out_fl:
        for protein_sequence in SeqIO.parse(protein_file_path,"fasta"):
            protein_id_list=protein_sequence.id.split("_",1)
            strain_id=protein_id_list[0]
            gff_protein_id=protein_id_list[1]
            if strain_id=="MGG":
                SeqIO.write(globals().get("MGG_db")[protein_sequence.id[0:9]],out_fl,"fasta")
            else:
                if strain_id=="WD-3-1":
                    strain_id=protein_sequence.id[0:8]
                    gff_protein_id=protein_sequence.id[9:]
                gene_plain_sequences=globals().get(strain_id+"_db")[gff_protein_id].sequence(
                    str(contig_path/(strain_id+".fasta")),
                    use_strand=True
                    )
                record = SeqRecord(
                    Seq(gene_plain_sequences,Bio.Alphabet.IUPAC.unambiguous_dna),
                    id=protein_sequence.id,
                    description=""
                    )
                SeqIO.write(record, out_fl,"fasta")
def mafft_popen(Single_Copy_Orthologue_file_path,gene_dir_path,strain_db_dir_path,contig_path,out_dir_path,out_stderr_dir_path):
    gene_file_path=gene_dir_path/(Single_Copy_Orthologue_file_path.stem+".fasta")
    if gene_file_path.exists() is False:
        extract_gene_gff_mafft(Single_Copy_Orthologue_file_path,contig_path,gene_file_path)
    ma_popen=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/miniconda3/bin/mafft",
            "--genafpair",
            "--maxiterate",
            "16",
            "--reorder",
            gene_file_path.resolve()
        ],
        stdout=(out_dir_path/Single_Copy_Orthologue_file_path.stem).open("w"),
        stderr=(out_stderr_dir_path/(Single_Copy_Orthologue_file_path.stem+"_stderr.txt")).open("w")
    )
    ma_popen.wait()
def mafft_run(Single_Copy_Orthologue_dir_path,gene_file_path,strain_db_dir_path,contig_path,mafft_out_dir_path,mafft_stderr_dir_path):
    '''
    input 1: Single_Copy_Orthologue_dir_path
    input 2: strain_db_dir_path
    input 3: contig_path
    output 1: gene_file_path(the output of extract_gene_gff function)
    output 1: mafft_out_dir_path
    outpur 2: mafft_stderr_dir_path
    '''
    global_names=globals()
    for db_file in strain_db_dir_path.iterdir():
        global_names[db_file.stem.strip()+"_db"]=gffutils.FeatureDB(db_file)
    global_names["MGG_db"]=SeqIO.index("../../70-15_refference_genome/magnaporthe_oryzae_70-15_8_genes.fasta","fasta")
    Parallel(n_jobs=12)(
        delayed(mafft_popen)(
            Single_Copy_Orthologue_file_path,
            gene_file_path,
            strain_db_dir_path,
            contig_path,
            mafft_out_dir_path,
            mafft_stderr_dir_path
            ) 
            for Single_Copy_Orthologue_file_path in Single_Copy_Orthologue_dir_path.iterdir()
        )