#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
rgx_percent=re.compile("C:(.+?)%")
rgx_strain=re.compile("10\.(.+)$")
busco_path=Path("/media/disk1/wangzhe/busco")
out_file="../read_busco.txt"
with open(out_file,"w+") as out_fl:    
    for busco_file in busco_path.glob("*/*.txt"):
        with busco_file.open() as busco_fl:
            for line in busco_fl:
                result_percent=rgx_percent.search(line)
                if result_percent is not None:
                    out_fl.write(
                        "{}\t{}\n".format(
                            rgx_strain.search(busco_file.stem).group(1),
                            result_percent.group(1)
                        )
                    )
                    break
