from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastformatterCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

def blastdb (in_file,db_file):
    make_db_cmd=NcbimakeblastdbCommandline(
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.0+/bin/makeblastdb',
        dbtype='nucl',
        input_file=in_file,
        out=db_file
    )
    make_db_cmd()
def blast(db_file,query_file,blast_out_asn_file,blast_out_xml_file,blast_out_txt_file):
    blast_asn_cmd=NcbiblastnCommandline(
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.0+/bin/blastn',
        query=query_file,
        db=db_file,
        outfmt=11,
        out=blast_out_asn_file
    )
    blast_asn_cmd()
    blast_xml_cmd=NcbiblastformatterCommandline(
        archive=blast_out_asn_file,
        outfmt=5,
        out=blast_out_xml_file,
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.0+/bin/blast_formatter'
    )
    blast_xml_cmd()
    blast_txt_cmd=NcbiblastformatterCommandline(
        archive=blast_out_asn_file,
        outfmt=7,
        out=blast_out_txt_file,
        cmd='/mnt/d/zhes_learning_space/software_in_ubuntu/ncbi-blast-2.10.0+/bin/blast_formatter'
    )
    blast_txt_cmd()
def read_blast_result(blast_out_xml_file,query_file,MGG_unpresent_Augustus):
    hits = []
    with open(MGG_unpresent_Augustus,'w+') as out_fl:
        for record in NCBIXML.parse(open(blast_out_xml_file)):
            if record.alignments:
                hits.append(record.query.split()[0])
            else:
                out_fl.write("{}\n".format(record.query.split()[0]))
    query_dict = SeqIO.index(query_file, "fasta")
    set_1=set(query_dict.keys())
    set_2=set(hits)
    misses = set_1 - set_2
    print(len(misses))
def MGG_unpresent_Augustus_locate_in_pav_orthofinder(
    MGG_unpresent_Augustus_list_file_name,
    pav_orthofinder_file_name,
    gene_protein_mapping_table_file_name,
    pav_MGG_unpresent_Augustus_file_name,
    orthofinder_unassianed_tsv_file_name,
    MGG_unpresent_Augustus_unassianed_list_file_name
    ):
    R_code='''
    MGG_unpresent_Augustus_locate_in_pav_orthofinder=function(MGG_unpresent_Augustus_list_file_name,pav_orthofinder_file_name,gene_protein_mapping_table_file_name,pav_MGG_unpresent_Augustus_file_name,orthofinder_unassianed_tsv_file_name,MGG_unpresent_Augustus_unassianed_list_file_name){
        require(readxl)
        require(WriteXLS)
        require(dplyr)
        MGG_unpresent_Augustus_list=read.table(MGG_unpresent_Augustus_list_file_name,stringsAsFactors = F)
        pav_orthofinde=read_xlsx(pav_orthofinder_file_name)
        pav_orthofinde=pav_orthofinde[,-1]

        gene_protein_mapping_table=read.table(gene_protein_mapping_table_file_name)

        pav_orthofinde_protein=merge(MGG_unpresent_Augustus_list,gene_protein_mapping_table,by.x = 1,by.y = 1)
        pav_MGG_unpresent_Augustus =pav_orthofinde %>% 
        filter(protein_id %in% pav_orthofinde_protein$V2)
        pav_MGG_unpresent_Augustus=pav_MGG_unpresent_Augustus[,c(158,1:157)]
        WriteXLS::WriteXLS(pav_MGG_unpresent_Augustus,pav_MGG_unpresent_Augustus_file_name)

        orthofinder_unassianed=read.table(orthofinder_unassianed_tsv_file_name,sep = "\t",header = T,check.names = F)
        MGG_unpresent_Augustus_unassianed_list=intersect(pav_MGG_unpresent_Augustus$protein_id,orthofinder_unassianed$`70-15_protein`)
        write.table(MGG_unpresent_Augustus_unassianed_list,MGG_unpresent_Augustus_unassianed_list_file_name,quote = F,row.names = F,col.names = F)
}
    '''
    R_MGG_unpresent_Augustus_locate_in_pav_orthofinder = SignatureTranslatedAnonymousPackage(R_code, "R_MGG_unpresent_Augustus_locate_in_pav_orthofinder")
    R_MGG_unpresent_Augustus_locate_in_pav_orthofinder.MGG_unpresent_Augustus_locate_in_pav_orthofinder(MGG_unpresent_Augustus_list_file_name,pav_orthofinder_file_name,gene_protein_mapping_table_file_name,pav_MGG_unpresent_Augustus_file_name,orthofinder_unassianed_tsv_file_name,MGG_unpresent_Augustus_unassianed_list_file_name)
def remove_MGG_unpresent_Augustus(pav_orthofinder,MGG_unpresent_Augustus_unassianed_list_file_name,pav_orthofinder_1574):
    '''
    用于从pav_orthofinder中删掉1574个出现在unassigned gene中的
    input 1: pav_orthofinder
    input 2: MGG_unpresent_Augustus_unassianed_list
    output 1: pav_orthofinder_1574
    '''
    R_code_remove_MGG_unpresent_Augustus='''
    R_remove_MGG_unpresent_Augustus=function(
    pav_orthofinder_file_name,
    MGG_unpresent_Augustus_unassianed_list_file_name,
    pav_orthofinder_1574_file_name
    ){
        require(readxl)
        require(WriteXLS)
        require(dplyr)
        MGG_unpresent_Augustus_unassianed_list=read.table(MGG_unpresent_Augustus_unassianed_list_file_name)
        pav_orthofinder=read_xlsx(pav_orthofinder_file_name)
        pav_orthofinder_1574=pav_orthofinder %>% 
            filter(!(protein_id %in% MGG_unpresent_Augustus_unassianed_list$V1))
        pav_orthofinder_1574=pav_orthofinder_1574[,-1]
        pav_orthofinder_1574=pav_orthofinder_1574[,c(1,158,2:157)]
        WriteXLS::WriteXLS(pav_orthofinder_1574,pav_orthofinder_1574_file_name)
    }
    '''
    R_remove_MGG_unpresent_Augustus = SignatureTranslatedAnonymousPackage(R_code_remove_MGG_unpresent_Augustus, "R_remove_MGG_unpresent_Augustus")
    R_remove_MGG_unpresent_Augustus.R_remove_MGG_unpresent_Augustus(str(pav_orthofinder),str(MGG_unpresent_Augustus_unassianed_list_file_name),str(pav_orthofinder_1574))