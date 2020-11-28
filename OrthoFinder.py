'''
@Author: your name
@Date: 2020-08-03 20:07:54
@LastEditTime: 2020-08-03 20:07:55
@LastEditors: your name
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/OrthoFinder.py
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

# protein_path="../Pan_genome_data/b_orthofinder_input_protein/"
protein_path="../Pan_genome_data_2/protein_base_no_duplicate/"

def start_from_protein():
    OrthoFinder=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/OrthoFinder/orthofinder",
            "-f",
            protein_path
        ],
        stderr=open("../Pan_genome_data_2/orthofinder_err.txt","w+"),
        stdout=open("../Pan_genome_data_2/orthofinder_out.txt","w+")
    )
    OrthoFinder.wait()

def remove_add_70_15():
    OrthoFinder_70_15=subprocess.Popen(
        [
            "orthofinder",
            "-b",
            "/mnt/c/WorkingDirectory",
            "-f",
            "../../wangzhe2/OrthoFinder_result/new_70_15/"
        ],
        stderr=open("../../wangzhe2/OrthoFinder_result/orthofinder_err.txt","w+"),
        stdout=open("../../wangzhe2/OrthoFinder_result/orthofinder_out.txt","w+")
    )
    OrthoFinder_70_15.wait()
start_from_protein()
# remove_add_70_15()
