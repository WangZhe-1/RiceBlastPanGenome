from rpy2.robjects.packages import importr
utils=importr("utils")
base=importr("base")
ortholog_joined_df=utils.read_table(
        "../Pan_genome_data/ortholog/joined_df.tsv",
        sep = "\t",
        header = True,
        **{'stringsAsFactors': False},
        **{'check.names': False}
        )
for i in range(1,(int(base.nrow(ortholog_joined_df)[0])+1)):
                df_row=ortholog_joined_df.rx(i, True)
                df_row_iter=iter(df_row)
                head_list=next(df_row_iter)
                for j in df_row_iter:
                    print(j[0])
                    # if len(head_list==1):
                    #     one2one(head_list[0],df_row)
                    # else:
                    #     two_head(head_list,df_row)