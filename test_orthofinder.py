'''
@Author: your name
@Date: 2020-08-04 13:11:20
@LastEditTime: 2020-08-04 13:54:03
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/test_orthofinder.py
'''
qHasAA=False
mLinesToCheck=100
with open("../Pan_genome_data_2/protein_base_no_duplicate/ 13FM-16-1.fasta", 'r') as fastaFile:
    for iLine, line in enumerate(fastaFile):
        if line.isspace(): continue
        if len(line) > 0 and line[0] == ">":
            pass
        else:
            line = line.upper()    # allow lowercase letters in sequences
            if not qHasAA and (iLine < mLinesToCheck):
                qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
    if not qHasAA:
        qOk = False
        print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences." % "13")