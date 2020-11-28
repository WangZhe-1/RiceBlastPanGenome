import subprocess
from pathlib import Path
from Directory_creater import directory_creater
from extract_strain_id import extract_strain_id
def nucmer(in_file_path,out_file_path,stdout_file_path,stderr_file_path):
    call_mummer=subprocess.Popen(
        [
            "nucmer",
            "--delta",
            str(out_file_path),
            "--maxmatch",
            "-c",
            "90",
            "-l",
            "40",
            "--threads=12",
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta",
            str(in_file_path)
        ],
        stdout=stdout_file_path.open('w'),
        stderr=stderr_file_path.open('w')
    )
    call_mummer.wait()
def call_filter(in_file,out_file,stderr_file):
    call_filter=subprocess.Popen(
        [
            "delta-filter",
            "-i",
            "89",
            "-l",
            "1000",
            "-1",
            str(in_file)
        ],
        stdout=out_file.open('w'),
        stderr=stderr_file.open('w')
    )
    call_filter.wait()
def snp(in_file,out_file,stderr_file):
    show_snp=subprocess.Popen(
        [
            "show-snps",
            "-Clr",
            "-x",
            "1",
            "-T",
            str(in_file)
        ],
        stdout=out_file.open('w'),
        stderr=stderr_file.open('w')
    )
    show_snp.wait()
def To_vcf(in_file,out_file,stderr,bgzip_out_err_file_path,index_out_err_file_path):
    '''
    include bgzip,index
    '''
    convert_to_vcf=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/zhe_MUMmerSNPs2VCF.py",
            str(in_file),
            str(out_file)
        ],
        stdout=stderr.open('w'),
        stderr=subprocess.STDOUT
    )
    convert_to_vcf.wait()
    bgzip=subprocess.Popen(
        [
            "bgzip",
            str(out_file)
        ],
        stdout=bgzip_out_err_file_path.open('w'),
        stderr=subprocess.STDOUT
    )
    bgzip.wait()
    index=subprocess.Popen(
        [
            "bcftools",
            "index",
            str(out_file)+'.gz'
        ],
        stdout=index_out_err_file_path.open('w'),
        stderr=subprocess.STDOUT
    )
    index.wait()
def call_snp_mummer(id_file,contig_dir_name,out_dir):
    '''
    input 1: strain 95 file
    input 2: contig_dir_name
    output 1: out dir
    '''
    strain_95_list=extract_strain_id(id_file)
    strain_95_list.append("ina168")
    strain_95_list.remove("magnaporthe_oryzae_70-15_8_proteins_T0")

    contig_dir_path=Path(contig_dir_name)
    intermediate_dir_path=directory_creater(out_dir/"intermediate_files")
    nucmer_stdout_dir_path=directory_creater(intermediate_dir_path/"nucmer_std_out")
    nucmer_stderr_dir_path=directory_creater(intermediate_dir_path/"nucmer_std_err")
    delta_dir_path=directory_creater(intermediate_dir_path/"delta")
    delta_filter_err_dir_path=directory_creater(intermediate_dir_path/"filter_err")
    delta_filter_file_dir_path=directory_creater(intermediate_dir_path/"filter")
    delta_snp_dir_path=directory_creater(intermediate_dir_path/"snp")
    delta_snp_err_path=directory_creater(intermediate_dir_path/"snp_err")
    delta_vcf_file_path=directory_creater(intermediate_dir_path/"vcf")
    delta_vcf_err_path=directory_creater(intermediate_dir_path/"To_vcf_err")
    merge_std_out_err_dir_path=directory_creater(intermediate_dir_path/"merge_std_out_err")
    bgzip_out_err_dir_path=directory_creater(intermediate_dir_path/"bgzip")
    index_out_err_dir_path=directory_creater(intermediate_dir_path/"index")
    global merge_list
    merge_list=["bcftools","merge"]
    for strain_95_id in strain_95_list:
        delta_file_path=delta_dir_path/(strain_95_id+".delta")
        if delta_file_path.is_file() is False:
            nucmer(
                str(contig_dir_path/(strain_95_id+".fasta")),
                str(delta_file_path),
                nucmer_stdout_dir_path/(strain_95_id+'_stdout.txt'),
                nucmer_stderr_dir_path/(strain_95_id+'_stderr.txt')
            )
        delta_filter_file_path=delta_filter_file_dir_path/(strain_95_id+"_filter.txt")
        if delta_filter_file_path.is_file() is False:
            call_filter(
                delta_file_path,
                delta_filter_file_path,
                delta_filter_err_dir_path/(strain_95_id+"_filter_err.txt")
            )
        delta_snp_file_path=delta_snp_dir_path/(strain_95_id+"_snp.txt")
        if delta_snp_file_path.is_file() is False:
            snp(
                delta_filter_file_path,
                delta_snp_file_path,
                delta_snp_err_path/(strain_95_id+"_snp_err.txt")
            )
        delta_snp2vcf_path=delta_vcf_file_path/(strain_95_id+".vcf")
        delta_snp2vcf_gz_path=delta_vcf_file_path/(strain_95_id+".vcf.gz")
        if delta_snp2vcf_gz_path.is_file() is False:
            To_vcf(
                delta_snp_file_path,
                delta_snp2vcf_path,
                delta_vcf_err_path/(strain_95_id+"_to_vcf_err.txt"),
                bgzip_out_err_dir_path/(strain_95_id+"_bgzip_out_err.txt"),
                index_out_err_dir_path/(strain_95_id+"_index_out_err.txt")
            )
        merge_list.append(str(delta_snp2vcf_path)+".gz")
    merge_list.extend(["-O","v","-o",str(out_dir/"merge.vcf")])
    call_merge=subprocess.Popen(
        merge_list,
        stdout=(merge_std_out_err_dir_path/"merge_out.txt").open('w'),
        stderr=subprocess.STDOUT
    )
    call_merge.wait()
# def merge_vcf(snp_mummer_result_path):
#     '''
#     input 1:snp_mummer_result_path
#     '''
#     vcf_dir_path=snp_mummer_result_path/"intermediate_files"/"vcf"
#     for vcf_file in vcf_dir_path.glob("*.gz"):
#         merge_list.append()
