'''
Author: your name
Date: 2020-09-28 15:50:28
LastEditTime: 2020-09-28 22:36:09
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /Pan_genome_github/fgenesh_gff_to_targer_dir.py
'''
import re
import subprocess
rgx_strain_id=re.compile("oryzae\s(?:isolate|strain)*(.+?)\s")
rgx_fr13=re.compile
def call_copy(from_file,to_dir,stdout_err_file_path):
    copy_popen=subprocess.Popen(
        [
            "cp",
            from_file,
            to_dir
        ],
        stdout=stdout_err_file_path.open('w'),
        stderr=subprocess.STDOUT
    )
    copy_popen.wait()
def fgenesh_gff_to_targer_dir(fgenesh_gff_dir_path,target_dir_path,gff_cp_stdout_err_dir_path):
    for file_path in fgenesh_gff_dir_path.rglob('result.txt'):
        if file_path.is_file():
            strain_id=None
            flag=True
            with file_path.open() as in_fl:
                for i,line in enumerate(in_fl):
                    if i==2:
                        search_result=rgx_strain_id.search(line)
                        if search_result is not None:
                            strain_id=search_result.group(1)
                            if strain_id=="genome":
                                strain_id="FR13"
                            elif strain_id=="70-15":
                                flag=False
                                break
                        elif re.search("Ina168",line) is not None:
                            strain_id="ina168"
                    elif i>2:
                        break
            if not flag:continue
            call_copy(file_path,target_dir_path/(strain_id+".gff"),gff_cp_stdout_err_dir_path/(strain_id+".txt"))
