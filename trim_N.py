#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

def sequence_cleaner_detect(in_file,out_file,out_ID,letter):
    '''
    1:in_file pan_protein/pan_gene
    2: out_file is a table
        ID N_count seq_len per
    3: out_ID is a list
        contain protein/gene which have "x"/"n"
    4:letter,protein is x,gene is n
    '''

    # Using the Biopython fasta parse we can read our fasta input
    with open(out_file,'w+') as fl_out:
        with open (out_ID,'w+') as out_ID_fl:
            fl_out.write('ID\tN_count\tseq_len\tper\n')
            for seq_record in SeqIO.parse(in_file, "fasta"):
                # Take the current sequence
                sequence = str(seq_record.seq).upper()
                # Check if the current sequence is according to the user parameters
                N_count=sequence.count(letter)
                if N_count > 0:
                    out_ID_fl.write('{}\n'.format(seq_record.id))
                seq_len=len(sequence)
                per=(float(N_count)/float(seq_len))*100
                fl_out.write('{}\t{}\t{}\t{}\n'.format(seq_record.id,N_count,seq_len,per))
