#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import re

protein_path=Path("../GFF/")
#protein_path=Path("../../GFF/")
rgx_strain_id=re.compile("(.+)_p")

for protein_file in protein_path.glob("*_protein.fasta"):
    # if protein_file.stem<"BTGP6F_protein":
    #     continue    
    strain_id=rgx_strain_id.search(protein_file.stem).group(1)
    busco_cmd=subprocess.Popen(
        [
            "busco",
            "-m",
            "protein",
            "-i",
            str(protein_file),
            "-o",
            "{}".format(strain_id),
            "-l",
            #"/mnt/d/zhes_learning_space/software_in_ubuntu/sordariomycetes_odb10",
            "/media/disk1/Database/busco_downloads/lineages/sordariomycetes_odb10",
            "-c",
            "10",
            "-f",
            "--offline"            
        ],
        stderr=open("../qc_contig_by_busco_err/{}.txt".format(strain_id),"w+"),
        stdout=open("../qc_contig_by_busco_out/{}.txt".format(strain_id),"w+")
    )
    busco_cmd.wait()
