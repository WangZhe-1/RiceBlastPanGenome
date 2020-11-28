import subprocess
import re
from Directory_creater import directory_creater
from extract_strain_id import extract_strain_id
from pathlib import Path

def replace_header(sam_file_name,out_sam_file):
    with sam_file_name.open('r') as in_fl:
        with out_sam_file.open('w') as out_fl:
            for line in in_fl:
                if line[0]=="@":
                    line=re.sub("@PG.+","@PG\tID:nucmer\tPN:nucmer\tVN:4.0\tCL:\"nucmer\"",line)
                    line=re.sub(" +","\t",line)
                    line=re.sub("VN1.0","VN:1.0",line)
                    
                out_fl.write(line)
def mummer_call(strain_id,out_file):
    call_mummer=subprocess.Popen(
        [
            "nucmer",
            "--sam-long",
            str(sam_out_dir/(strain_id+".sam")),
            "--threads=12",
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta",
            str(contig_dir_path/(strain_id+".fasta"))
        ],
        stdout=open(str(mummer_stdout_dir/(strain_id+"_out.txt")),'w'),
        stderr=open(str(mummer_stderr_dir/(strain_id+"_err.txt")),'w')
    )
    call_mummer.wait()
def sam2bam(in_file,out_file,err_file):
    call_sam2bam=subprocess.Popen(
        [
            "samtools",
            "view",
            "-t",
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta.fai",
            "-S",
            "-b",
            str(in_file)
        ],
        stdout=open(out_file,'w+'),
        stderr=open(err_file,'w+')
    )
    call_sam2bam.wait()
def add_RG(in_file,out_file,stdout_file_path,stderr_file_path):
    call_add_RG=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
            "AddOrReplaceReadGroups",
            "-I",
            str(in_file),
            "-O",
            str(out_file),
            "-LB",
            "lib1",
            "-PU",
            "unit1",
            "-PL",
            "UNKNOWN",
            "-SM",
            in_file.stem
        ],
        stdout=stdout_file_path.open('w'),
        stderr=stderr_file_path.open('w')
    )
    call_add_RG.wait()
def bam_sort(in_file,out_file,stdout_file_path,err_file):
    call_sam_sort=subprocess.Popen(
        [
            "samtools",
            "sort",
            "-O",
            "bam",
            "-o",
            str(out_file),
            str(in_file)
        ],
        stdout=stdout_file_path.open('w'),
        stderr=err_file.open('w')
    )
    call_sam_sort.wait()
def bam_index(in_file,stdout_file_path,err_file):
    call_sam_index=subprocess.Popen(
        [
            "samtools",
            "index",
            str(in_file)
        ],
        stdout=stdout_file_path.open('w'),
        stderr=err_file.open('w')
    )
    call_sam_index.wait()
def haplotypecaller(in_file,out_file,stdout_file_path,stderr_file_path):
    call_haplotypecaller=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
            "--java-options",
            "-Xmx20g",
            "HaplotypeCaller",
            # "--disable-read-filter",
            # "WellformedReadFilter",
            # # "{WellformedReadFilter,MappingQualityReadFilter}",
            # "--disable-read-filter",
            # "MappingQualityReadFilter",
            "--disable-tool-default-read-filters",
            "-R",
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta",
            "-I",
            str(in_file),
            "-O",
            str(out_file),
            "-ERC",
            "GVCF"
        ],
        stdout=stdout_file_path.open('w'),
        stderr=stderr_file_path.open('w')
    )
    call_haplotypecaller.wait()
def call_snp_HaplotypeCaller(id_file,contig_dir_name,out_dir):
    '''
    input 1: strain 95 file
    input 2: contig_dir_name
    output 1: out dir
    '''
    global mummer_stderr_dir,mummer_stdout_dir,contig_dir_path,sam_out_dir
    sam_out_dir=directory_creater(out_dir/"mummer_sam_files")
    raw_bam_out_dir=directory_creater(out_dir/"raw_bam_files")
    RG_bam_dir=directory_creater(out_dir/"RG_bam_files")
    intermediate_dir=directory_creater(out_dir/"intermediate_files")
    mummer_stderr_dir=directory_creater(intermediate_dir/"mummer_err")
    mummer_stdout_dir=directory_creater(intermediate_dir/"mummer_out")
    samview_stderr_dir=directory_creater(intermediate_dir/"sam_view_err")
    add_RG_stdout_dir=directory_creater(intermediate_dir/"add_RG_out")
    add_RG_stderr_dir=directory_creater(intermediate_dir/"add_RG_err")
    strain_95_list=extract_strain_id(id_file)
    strain_95_list.append("ina168")
    strain_95_list.remove("magnaporthe_oryzae_70-15_8_proteins_T0")
    contig_dir_path=Path(contig_dir_name)
    replace_header_sam_dir=directory_creater(intermediate_dir/"replace_header_sam")
    bam_sort_dir_path=directory_creater(intermediate_dir/"bam_sort")
    bam_sort_stdout_dir_path=directory_creater(intermediate_dir/"bam_sort_out")
    bam_sort_stderr_dir_path=directory_creater(intermediate_dir/"bam_sort_err")
    bam_index_stdout_dir_path=directory_creater(intermediate_dir/"bam_index_out")
    bam_index_stderr_dir_path=directory_creater(intermediate_dir/"bam_index_err")
    gvcf_dir=directory_creater(out_dir/"gvcf")
    haplotypecaller_stdout_dir_path=directory_creater(intermediate_dir/"haplotypecaller_stdout")
    haplotypecaller_stderr_dir_path=directory_creater(intermediate_dir/"haplotypecaller_stderr")
    for strain_95_id in strain_95_list:
        mummer_sam_out=sam_out_dir/(strain_95_id+".sam")
        if mummer_sam_out.is_file() is False:
            mummer_call(
                strain_95_id,
                str(mummer_sam_out)
                )
        replace_header_sam_file_path=replace_header_sam_dir/(strain_95_id+".sam")
        if replace_header_sam_file_path.is_file() is False:
            replace_header(mummer_sam_out,replace_header_sam_file_path)
        raw_bam_file_path=raw_bam_out_dir/(strain_95_id+".bam")
        if raw_bam_file_path.is_file() is False:
            sam2bam(replace_header_sam_file_path,raw_bam_file_path,samview_stderr_dir/(strain_95_id+"_samview_err.txt"))
        RG_bam_file_path=RG_bam_dir/(strain_95_id+".bam")
        if RG_bam_file_path.is_file() is False:
            add_RG(
                raw_bam_file_path,
                RG_bam_file_path,
                add_RG_stdout_dir/(strain_95_id+"_out.txt"),
                add_RG_stderr_dir/(strain_95_id+"_err.txt")
                )
        bam_RG_sort_file_path=bam_sort_dir_path/(strain_95_id+"_sort.bam")
        if bam_RG_sort_file_path.is_file() is False:
            bam_sort(
                RG_bam_file_path,
                bam_RG_sort_file_path,
                bam_sort_stdout_dir_path/(strain_95_id+"_out.txt"),
                bam_sort_stderr_dir_path/(strain_95_id+"_err.txt")
                )
            bam_index(
                bam_RG_sort_file_path,
                bam_index_stdout_dir_path/(strain_95_id+"_out.txt"),
                bam_index_stderr_dir_path/(strain_95_id+"_err.txt")
            )
        gvcf_file_path=gvcf_dir/(strain_95_id+".vcf")
        if gvcf_file_path.is_file() is False:
            haplotypecaller(
                bam_RG_sort_file_path,
                gvcf_file_path,
                haplotypecaller_stdout_dir_path/(strain_95_id+"_out.txt"),
                haplotypecaller_stderr_dir_path/(strain_95_id+"_err.txt")
                )
def GenomicsDBImport(snp_dir):
    '''
    input 1:snp_dir generate by call snp
    '''
    gvcf_dir_path=snp_dir/"gvcf"
    vcfdb_dir=snp_dir/"vcfdb"
    std_out_err=directory_creater(snp_dir/"intermediate_files"/"GenomicsDBImport_std_out_err")
    call_list=[
        "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
        "--java-options",
        "-Xmx20g -Xms20g",
        "GenomicsDBImport",
        "--genomicsdb-workspace-path",
        str(vcfdb_dir)
    ]
    for vcf_file in gvcf_dir_path.glob("*.vcf"):
        call_list.append("-V")
        call_list.append(str(vcf_file))
    call_GenomicsDBImport=subprocess.Popen(
        call_list,
        stdout=(std_out_err/"GenomicsDBImport_stdout.txt").open('w'),
        stderr=(std_out_err/"GenomicsDBImport_stderr.txt").open('w')
    )
    call_GenomicsDBImport.wait()
def CombineGVCFs(snp_dir):
    '''
    input 1:snp_dir generate by call snp
    '''
    vcfdb_dir=directory_creater(snp_dir/"CombineGVCFs_vcf")
    std_out_err=directory_creater(snp_dir/"intermediate_files"/"CombineGVCFs_std_out_err")
    call_list=[
        "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
        "--java-options",
        "-Xmx20g -Xms20g",
        "CombineGVCFs",
        "-R",
        "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta",
        "-O",
        str(vcfdb_dir/"pan_genome_CombineGVCFs.vcf")
    ]
    for vcf_file in gvcf_dir_path.glob("*.vcf"):
        call_list.append("-V")
        call_list.append(str(vcf_file))
    call_CombineGVCFs=subprocess.Popen(
        call_list,
        stdout=(std_out_err/"CombineGVCFs_stdout.txt").open('w'),
        stderr=(std_out_err/"CombineGVCFs_stderr.txt").open('w')
    )
    call_CombineGVCFs.wait()
def GenotypeGVCFs(snp_dir):
    input_vcf_file=snp_dir/"CombineGVCFs_vcf"/"pan_genome_CombineGVCFs.vcf"
    output_vcf_file=snp_dir/"call_snp_gatk.vcf"
    std_out_err_dir=directory_creater(snp_dir/"intermediate_files"/"GenotypeGVCFs_out_err")
    call_GenotypeGVCFs=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
            "--java-options",
            "-Xmx20g -Xms20g",
            "GenotypeGVCFs",
            "-R",
            "/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta",
            "-V",
            str(input_vcf_file),
            "-O",
            str(output_vcf_file)
        ],
        stdout=(std_out_err_dir/"GenotypeGVCFs_out.txt").open('w'),
        stderr=(std_out_err_dir/"GenotypeGVCFs_err.txt").open('w')
    )
    call_GenotypeGVCFs.wait()
def print_WellformedReadFilter(in_file,out_file,filter):
    call_print=subprocess.Popen(
        [
            "/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk",
            "PrintReads",
            "--disable-tool-default-read-filters",
            "--read-filter",
            # "--disable-read-filter",
            filter,
            "-I",
            str(in_file),
            "-O",
            str(out_file)
        ],
        stdout=open("../Pan_genome_data/print_stdout"+'_'+filter+".txt",'w'),
        stderr=open("../Pan_genome_data/print_stderr"+'_'+filter+".txt",'w')
    )
    call_print.wait()
if __name__ == '__main__':
    # replace_header(Path("../Pan_genome_data/out_sam_long_13FM-16-1.sam"),Path("../Pan_genome_data/out_sam_py.sam"))
    # sam2bam("../Pan_genome_data/out_sam_py.sam","../Pan_genome_data/out_sam_py.bam","../Pan_genome_data/out_sam2bam_err.txt")
    # add_RG(Path("../Pan_genome_data/out_sam_py.bam"),"../Pan_genome_data/out_sam_py_RG.bam",Path("../Pan_genome_data/RG_out.txt"),Path("../Pan_genome_data/RG_err.txt"))
    # bam_sort("../Pan_genome_data/out_sam_py_RG.bam","../Pan_genome_data/out_sam_py_RG_sort.bam",Path("../Pan_genome_data/sort_out.txt"),Path("../Pan_genome_data/sort_err.txt"))
    # bam_index("../Pan_genome_data/out_sam_py_RG_sort.bam",Path("../Pan_genome_data/index_out.txt"),Path("../Pan_genome_data/index_err.txt"))
    # haplotypecaller("../Pan_genome_data/out_sam_py_RG_sort.bam","../Pan_genome_data/out_sam_py_RG_sort.vcf",Path("../Pan_genome_data/HaplotypeCaller_out.txt"),Path("../Pan_genome_data/HaplotypeCaller_err.txt"))
    filter_list=[
        "ValidAlignmentStartReadFilter",
        "ValidAlignmentEndReadFilter",
        "AlignmentAgreesWithHeaderReadFilter",
        "HasReadGroupReadFilter",
        "MatchingBasesAndQualsReadFilter",
        "ReadLengthEqualsCigarLengthReadFilter",
        "SeqIsStoredReadFilter",
        "CigarContainsNoNOperator"
    ]
    # for filter_cat in filter_list:
    #     # print_WellformedReadFilter("../Pan_genome_data/call_snp/intermediate_files/bam_sort/13FM-16-1_sort.bam","../Pan_genome_data/"+filter_cat+".sam",filter_cat)
    #     print_WellformedReadFilter("../Pan_genome_data/out_sam_py_RG_sort.bam","../Pan_genome_data/"+filter_cat+".sam",filter_cat)
    print_WellformedReadFilter("../Pan_genome_data/out_sam_py_RG_sort.bam","../Pan_genome_data/MappingQualityReadFilter.sam","MappingQualityReadFilter")