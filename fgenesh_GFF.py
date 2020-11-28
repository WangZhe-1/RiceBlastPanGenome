'''
@Author: your name
@Date: 2020-08-04 17:08:28
@LastEditTime: 2020-08-04 17:33:00
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/fgenesh_GFF.py
'''

def writen_name(in_dir_path,out_list_file_name):
    with out_list_file_name.open('w') as out_fl:
        for contig_file in in_dir_path.iterdir():
            out_fl.write(contig_file.stem+"\n")
