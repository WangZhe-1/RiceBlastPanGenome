#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
pan_protein=Path('./pan_protein.fasta')
out_file='./'


interpro_cmd=subprocess.Popen(
    [
        "/media/disk1/Biosoft/interproscan-5.44-79.0/interproscan.sh",
        "-i",
        pan_protein,
        "-f",
        "GFF3",
        "-d",
        out_file,
        "-appl",
        "Pfam",
        "-goterms",
        "-pa",
        "-iprlookup",
        "-cpu",
        "20"
    ],
    stdout=open('./inperpro_out','w+'),
    stderr=open('./interpro_err',"w+"),
    universal_newlines=True,
    encoding='UTF-8',
)
interpro_cmd.wait()
