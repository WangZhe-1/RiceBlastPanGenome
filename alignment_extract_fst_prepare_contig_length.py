from pathlib import Path
import pysam
from extract_strain_id import extract_strain_id
import json
def contig_length2json(contig_path_name,id_file,json_dir_path_name):
    '''
    input 1: contig_path_name
    input 2: id_file
    output 1: json_dir_path_name
    '''
    species_95_list=extract_strain_id(id_file)
    species_95_list.append("70-15")
    species_95_list.append("ina168")
    species_95_list.remove("magnaporthe_oryzae_70-15_8_proteins_T0")
    contig_path=Path(contig_path_name)
    for strain_id in species_95_list:
        for strain_file in contig_path.glob(strain_id+".fasta"):
            with pysam.FastaFile(strain_file) as contig_length_it:
                contig_length_dic={}
                for contig_id,contig_length in zip(contig_length_it.references,contig_length_it.lengths):
                    contig_length_dic[contig_id]=(0, int(contig_length))
            json_file_name=json_dir_path_name/(strain_file.stem+".json")
            with open(json_file_name,'w+') as json_fl:
                json.dump(contig_length_dic,json_fl)



