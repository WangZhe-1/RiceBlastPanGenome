'''
@Author: your name
@Date: 2020-07-26 13:00:13
@LastEditTime: 2020-07-26 14:07:20
@LastEditors: Please set LastEditors
@Description: In User Settings Edit
@FilePath: /Pan_genome_github/linsdu.py
'''
nucmer --sam-short=./out_sam_short_13FM-16-1.sam --threads=12 /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/contig/13FM-16-1.fasta
nucmer --delta=../Pan_genome_data/test.delta --threads=12 /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/contig/13FM-16-1.fasta
show-snps -T test_only.filter
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/my-mummer-2-vcf.py --snps test_only_snp --input-header --reference /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta > test_only_vcf
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R 
delta-filter -i 89 -l 1000 -1 test.delta > test_89_1000.filter
delta-filter -1 test.delta > test_only.filter
show-coords -c test_89_1000.filter > test_89_1000.filter.coords
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_89_1000.filter.coords -o test_89_1000.plot -m 1000 -q 300000 -s -t -l 

#没有maxmatch，-i 79
delta-filter -i 79 -l 1000 -1 test.delta > test_79_1000.filter
show-coords -c test_79_1000.filter > test_79_1000.filter.coords
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_79_1000.filter.coords -o test_79_1000.plot -m 1000 -q 300000 -s -t -l 
#画图没有t
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_79_1000.filter.coords -o test_79_1000.plot_t -m 1000 -q 300000 -s -l 

nucmer --delta=../Pan_genome_data/test_maxmatch.delta --maxmatch --threads=12 /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/contig/13FM-16-1.fasta
delta-filter -i 89 -l 1000 -1 test_maxmatch.delta > test_89_1000_maxmatch.filter
delta-filter -1 test_maxmatch.delta > test_only_maxmatch.filter
show-coords -c test_89_1000_maxmatch.filter > test_89_1000_maxmatch.filter.coords
show-coords -c test_only_maxmatch.filter > test_only_maxmatch.filter.filter.coords
show-coords -c test_maxmatch.delta > test_maxmatch.delta.coords
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_89_1000_maxmatch.filter.coords -o test_89_1000_maxmatch.plot -m 1000 -q 300000 -s -t -l 
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_only_maxmatch.filter -o test_only_maxmatch.plot -m 1000 -q 300000 -s -t -l 
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_maxmatch.delta.coords -o test_all_maxmatch.plot -m 1000 -q 300000 -s -t -l

#maxmatch+–c 90 –l 40
nucmer --delta=../Pan_genome_data/test_maxmatch–c90–l40.delta --maxmatch –c 90 –l 40 --threads=12 /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/contig/13FM-16-1.fasta
delta-filter -i 89 -l 1000 -1 test_maxmatch–c90–l40.delta > test_89_1000_maxmatch–c90–l40.filter
show-coords -c test_89_1000_maxmatch–c90–l40.filter > test_89_1000_maxmatch–c90–l40.filter.coords
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/mummerCoordsDotPlotly.R -i test_89_1000_maxmatch–c90–l40.filter.coords -o test_89_1000_maxmatch–c90–l40.plot -m 1000 -q 300000 -s -t -l 
show-snps -T test_89_1000_maxmatch–c90–l40.filter>test_89_1000_maxmatch–c90–l40.filter.snp
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/my-mummer-2-vcf.py --snps test_89_1000_maxmatch–c90–l40.filter.snp --input-header --reference /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta > test_89_1000_maxmatch–c90–l40.filter.vcf
#header
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/my-mummer-2-vcf.py --snps test_89_1000_maxmatch–c90–l40.filter.snp --output-header --reference /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta > test_89_1000_maxmatch–c90–l40_header.filter.vcf
/mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_github/my-mummer-2-vcf.py --snps test_89_1000_maxmatch–c90–l40.filter.snp --input-header --output-header --reference /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta > test_89_1000_maxmatch–c90–l40_header.filter.vcf

#vcftools部分
vcftools --vcf test_snp_10000 --TajimaD 3000 --out test_snp_10000_TajimaD_out
vcftools --vcf test_snp_10000 --window-pi 3000 --out test_snp_10000_pi
修改header，
补加参考基因组信息
samtools view -t /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta.fai -S -b /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_short_amend_13FM-16-1.sam> /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_short_13FM-16-1.bam 
gatk AddOrReplaceReadGroups
/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk AddOrReplaceReadGroups -I /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_short_13FM-16-1.bam -O /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_RG.bam -LB lib1 -PU unit1 -PL UNKNOWN -SM pan_data
/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk --java-options "-Xmx20g" HaplotypeCaller -R /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta -I /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_RG.bam -O output.g.vcf -ERC GVCF

samtools sort -O bam -o out_sam_RG.sorted.bam out_sam_RG.bam
samtools index /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_RG.sorted.bam
/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk --java-options "-Xmx20g" HaplotypeCaller -R /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/70-15_refference_genome/70-15_supercontigs.fasta -I /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_RG.sorted.bam -O output.g.vcf -ERC GVCF

/mnt/d/zhes_learning_space/software_in_ubuntu/gatk-4.1.8.0/gatk AddOrReplaceReadGroups -I /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_py.bam -O /mnt/d/zhes_learning_space/the_assignment/pan_genome/Mag_genomes/wangzhe2/Pan_genome_data/out_sam_py_out.bam -LB lib1 -PU unit1 -PL UNKNOWN -SM pan_data
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Homo_sapiens_assembly38.fasta \
   -V gendb://my_database \
   -O output.vcf.gz \
   --tmp-dir=/path/to/large/tmp
vcf_file=call_snp/"pan_genome.vcf"
gatk CombineGVCFs \
   -R reference.fasta \
   --variant sample1.g.vcf.gz \
   --variant sample2.g.vcf.gz \
   -O cohort.g.vcf.gz
vcftools --vcf call_snp_mummer/merge.vcf --window-pi --bed nucleotide_diversity/core_100.bed --out nucleotide_diversity/core_100
/mnt/d/zhes_learning_space/software_in_ubuntu/ParaAT2.0/ParaAT.pl -h ../Pan_genome_data/prepare_for_ParaAT/homolog/MGG_01191T0.txt -n ../Pan_genome_data/prepare_for_ParaAT/nucleotide/MGG_01191T0.fasta -a ../Pan_genome_data/prepare_for_ParaAT/aminoacid/MGG_01191T0.fasta -p ../Pan_genome_data/prepare_for_ParaAT/th.txt -m muscle -f axt -g -k -o ../Pan_genome_data/result_dir_ziji
"../Pan_genome_data/prepare_for_ParaAT/homolog/MGG_01191T0.txt"
"../Pan_genome_data/prepare_for_ParaAT/nucleotide/MGG_01191T0.fasta"
"../Pan_genome_data/prepare_for_ParaAT/aminoacid/MGG_01191T0.fasta"
"../Pan_genome_data/prepare_for_ParaAT/th.txt"
/mnt/d/zhes_learning_space/software_in_ubuntu/ParaAT2.0/ParaAT.pl -h ../Pan_genome_data/prepare_for_ParaAT/homolog/MGG_02047T0.txt -n ../Pan_genome_data/prepare_for_ParaAT/nucleotide/MGG_02047T0.fasta -a ../Pan_genome_data/prepare_for_ParaAT/aminoacid/MGG_02047T0.fasta -p ../Pan_genome_data/prepare_for_ParaAT/th.txt -m muscle -f axt -g -k -o ../Pan_genome_data/result_dir_ziji_duan
/mnt/d/zhes_learning_space/software_in_ubuntu/ParaAT2.0/ParaAT.pl -h ../Pan_genome_data/prepare_for_ParaAT/homolog/MGG_00040T0.txt -n ../Pan_genome_data/prepare_for_ParaAT/nucleotide/MGG_00040T0.fasta -a ../Pan_genome_data/prepare_for_ParaAT/aminoacid/MGG_00040T0.fasta -p ../Pan_genome_data/prepare_for_ParaAT/th.txt -m muscle -f axt -g -k -t -o ../Pan_genome_data/result_dir_MGG_00040T0_t
/mnt/d/zhes_learning_space/software_in_ubuntu/ParaAT2.0/ParaAT.pl -h ../Pan_genome_data/prepare_for_ParaAT/homolog/MGG_00040T0.txt -n ../Pan_genome_data/prepare_for_ParaAT/nucleotide/MGG_00040T0.fasta -a ../Pan_genome_data/prepare_for_ParaAT/aminoacid/MGG_00040T0.fasta -p ../Pan_genome_data/prepare_for_ParaAT/th.txt -m muscle -f axt -k -o ../Pan_genome_data/result_dir_MGG_00040T0_with_gap
"../Pan_genome_data/annotation/secreted_proteins/new_tpye_secreted_protein.txt"
"../Pan_genome_data/annotation/secreted_proteins/classic_secreted_protein.txt"
"../Pan_genome_data/pan_categories_protein/core.txt"
"../Pan_genome_data/R_result/strain_clade_category_ID.txt"
"../Pan_genome_data/avr_2/HZVJZGPH01R-Alignment-Descriptions.csv"