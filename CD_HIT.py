'''
@Author: your name
@Date: 2020-08-01 15:24:00
@LastEditTime: 2020-08-04 13:24:22
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/CD-HIT.py
'''
import subprocess
from Bio import SeqIO
def cd_hit(input_gene_dir_path,in_protein_dir_path,output_gene_dir_path,output_protein_dir_path,std_out_dir,std_err_dir):
    for input_gene_file in input_gene_dir_path.iterdir():
        if input_gene_file.stem == "70-15": continue
        output_gene_file_path=output_gene_dir_path/(input_gene_file.stem+".fasta")
        std_out=std_out_dir/(input_gene_file.stem+"_out.txt")
        std_err=std_err_dir/(input_gene_file.stem+"_err.txt")
        call_cd_hit=subprocess.Popen(
            [
                "/mnt/d/zhes_learning_space/software_in_ubuntu/cd-hit-v4.8.1-2019-0228/cd-hit-est",
                "-i",
                str(input_gene_file),
                "-o",
                str(output_gene_file_path),
                "-c",
                "1",
                "-aS",
                "1",
                "-n",
                "10",
                "-T",
                "12",
                "-d",
                "0"
            ],
            stderr=open(std_err,'w'),
            stdout=open(std_out,'w')
        )
        call_cd_hit.wait()
        if call_cd_hit.returncode==0:
            out_protein_file_path=output_protein_dir_path/(input_gene_file.stem.strip()+".fasta")
            with output_gene_file_path.open() as output_gene_fl:
                input_file_index=SeqIO.index(str(in_protein_dir_path/(input_gene_file.stem+".fasta")),'fasta')
                with out_protein_file_path.open('w') as out_protein_fl:
                    for seq_record in SeqIO.parse(output_gene_fl,"fasta"):
                        SeqIO.write(input_file_index[seq_record.id],out_protein_fl,"fasta")
        else:
            print(str(input_gene_file))


# if __name__ == "__main__":
#     cd_hit("../Pan_genome_data/TW-6-2-2-B-1_gene_5.fasta","../Pan_genome_data/TW-6-2-2-B-1_gene_5_cd_hit_1","../Pan_genome_data/TW-6-2-2-B-1_gene_5_cd_hit_std_err_out")
    # cd_hit_protein("../Pan_genome_data/TW-6-2-2-B-1_protein_verify.fasta","../Pan_genome_data/TW-6-2-2-B-1_protein_verify_1","../Pan_genome_data/TW-6-2-2-B-1_gene_5_cd_hit_std_err_out")

# from multiprocessing.dummy import Pool # use threads
# from subprocess import check_output

# def md5sum(filename):
#     try:
#         return check_output(["md5sum", filename]), None
#     except Exception as e:
#         return None, e

# if __name__ == "__main__":
#     p = Pool(number_of_processes) # specify number of concurrent processes
#     with open("md5sums.txt", "wb") as logfile:
#         for output, error in p.imap(md5sum, filenames): # provide filenames
#             if error is None:
#                logfile.write(output)