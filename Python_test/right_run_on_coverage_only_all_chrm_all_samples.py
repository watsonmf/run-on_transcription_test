import pandas as pd
import string
import numpy as np
import glob
import os
import csv
#Check files to make sure that the reads are going in the right direction. (pos or min)

def header_int_to_char(df):
    """this function will change all numbers to letters for pandas"""
    import string
    x = 0  #rename int headers to letters
    for num in (list(df)):
        df.rename(columns={int(num) : str(string.ascii_uppercase[x])}, inplace=True)
        x = x + 1
def get_cov_perc(df, chr, start, end):
    """this function will take start and end and return the percent coverage over the area.
    for a BedGraph file for only one chromosome"""
    rows_df = df[ (df["end"] > start) & (df["start"] < end) ]
    if len(rows_df) == 0:
        return 0.0
    else:
        start_no_cov_length = 0
        if start > rows_df["start"].iloc[0]:
            start_no_cov_length = start - rows_df["start"].iloc[0]

        end_no_cov_length = 0
        if end < rows_df["end"].iloc[-1]:
            end_no_cov_length = rows_df["end"].iloc[-1] - end

        cov_sum = rows_df["length"].sum() - (start_no_cov_length + end_no_cov_length)    #add new start and end distance and all rows inbetween
        cov_length = end - start
        if cov_length == 0:
            cov_per = 0.0
            return cov_per
        else:
            cov_per = cov_sum / cov_length
            return cov_per
def get_chrm_dic(gtf, plu_or_min, type, gene_biotype):
    """This function will take a gtf and return strand specific dictionary of different chrm"""
    my_col = ["chr","source","type","start","end","dot","strand","dot2","gene_id"]
    df_gtf = pd.read_csv(gtf, sep="\t",header=None,names=my_col)
    del df_gtf["source"]
    del df_gtf["dot"]
    del df_gtf["dot2"]
    df_gtf["length"] = df_gtf["end"] - df_gtf["start"]
    df_gtf = df_gtf[~df_gtf.chr.str.contains("MtDNA")] # Filter out bogus tRNA
    #Make dictionary of exons
    names=df_gtf['type'].unique().tolist()
    d_type = {type : df_gtf.loc[df_gtf.type==type] for type in names}
    df_gene = d_type[type] #Input for type
    print("Length of gene only gtf is", len(df_gene))
    #Think about not doing this next step because I can now pick up between genes. (To find ncRNA)
    ##split the last col of gtf file into a new df seperated by ;  merge back to gtf
    df_gene_id = df_gene["gene_id"].str.split(';',expand=True)
    header_int_to_char(df_gene_id)
    df_merge = pd.merge(
        df_gene,
        df_gene_id,
        how="left",
        left_index = True,
        right_index = True,
    )
    #Grab all genes that say protein coding
    print(df_merge["D"])
    names=df_merge['D'].unique().tolist()
    print(names)
    d_prot = {prot : df_merge.loc[df_merge.D==prot] for prot in names}
    df_prot = d_prot[" gene_biotype \"protein_coding\""]  #Input for gene_biotype
    print("Length of gene only everything that is protein coding is", len(df_prot))
    del df_prot["E"]
    del df_prot["F"]
    del df_prot["G"]
    del df_prot["C"]
    del df_prot["A"]
    del df_prot["B"]
    del df_prot["D"]

    #Make positive and negative strand dictionary of dataframe
    names=df_prot['strand'].unique().tolist()
    d_strand = {strand : df_prot.loc[df_prot.strand==strand] for strand in names}
    df_plu_or_min = d_strand[plu_or_min]
    print("Lenth of min gtf is ", len(df_plu_or_min))
    #Make df of chrom
    names=df_plu_or_min['chr'].unique().tolist()
    d_gtf_chr = {chrom : df_plu_or_min.loc[df_plu_or_min.chr==chrom] for chrom in names}
    return d_gtf_chr
def get_right_run_on(df):
        """This function will take a dataframe and return a bed file of right run on and gtf of all genes that have run on"""
        print("Length of min gene only everything that is protein coding chrm", chrom," is", len(df))
        df = df.copy()
        df = df.reset_index()
        #Drop first and last row (mess with boolean filtering)
        df.fillna(value=0.0,inplace=True)    #fill with na
        print("Length of df with fill na values", len(df))
        #Make dataframes of genes that overlap, and are between other genes
        df_inside_gene = df[ ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )]         #this finds all genes inside bigger genes
        print("Number of genes inside other genes", len(df_inside_gene))
        df_left_overlap = df[ df["start"] < df["end"].shift(1) ]
        print("Number of genes with left overlap", len(df_left_overlap))
        df_right_overlap = df[ df["start"].shift(-1) < df["end"] ].copy()
        print("Number of genes with right overlap", len(df_right_overlap))
        df_both_end_overlap = df[ ( df["start"] < df["end"].shift(1) ) & ( df["start"].shift(-1) < df["end"] )]
        print("Number of genes with both end overlap", len(df_both_end_overlap))

        #For positive (take out genes inside of other genes and that have right overlap)
        df["inside_gene"] = ( df["start"].shift(1) < df["start"] ) & ( df["end"] < df["end"].shift(1) )
        df = df[ ~df.inside_gene ]
        print("Number of genes without inside genes", len(df))

        df["right_int_length"] = df["start"].shift(-1) - df["end"]       #Calculate 3p int after you get rid of inside genes
        df["right_overlap"] = df["start"].shift(-1) < df["end"]       #Get rid of genes that have another gene in 3' end. (Keep while calculating 3p int length because 5' end is important)
        df = df[ ~df.right_overlap ]       #drop the values
        print("Number of genes without right overlap", len(df))

        #Drop some of the columns not needed
        del df["inside_gene"]
        del df["right_overlap"]

        # Run on Loop (with coverage)
        df["run_on_length"] = 0        #Initialize run on length column
        df_final = pd.DataFrame()      #Initialize final dataframe
        window_size = 100

        #For everything less than window size, give full length if over 95 cov. else 0
        df_LT_w = df [ df['right_int_length'] <= window_size ].copy()    #Less than window size
        print("Number of genes under intergenic ", window_size, "window size = ", len(df_LT_w))
        df_LT_w['cov_perc'] = df_LT_w.apply(lambda row: get_cov_perc(d_bed_chr[row['chr']], row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
        df_LT_w["run_on_length"] = np.where(df_LT_w['cov_perc'] > 0.95 , df_LT_w['right_int_length'], 0 )
        df_final = df_final.append(df_LT_w)

        #For genes over window size
        df.fillna(value=100000,inplace=True) #Give the last gene a 100kb dummy length (Or put chrom sizes in later)
        df = df [ df['right_int_length'] > window_size ]
        print("Number of genes over intergenic ", window_size, "window size = ", len(df))
        counter = window_size
        while len(df) > 0:
            df['cov_perc'] = df.apply(lambda row: get_cov_perc(d_bed_chr[row['chr']], row['chr'], row['end'], row['end'] + window_size), axis=1)
            df_under_95 = df [ df['cov_perc'] < 0.95  ]
            if len(df_under_95) > 0:
                df_under_95 = df [ df['cov_perc'] < 0.95  ].copy()
                df_under_95["run_on_length"] = df["run_on_length"] + window_size - counter
                df_final = df_final.append(df_under_95)
            df = df [ df['cov_perc'] > 0.95  ]
            window_size = window_size + counter

        df_final = df_final.append(df_right_overlap)   #Add back genes in that had right overlap
        df = df_final                                  #Change name back to df_final for easier coding
        df.fillna(value=0.0,inplace=True)
        df = df.sort_values("start")

        print("Number of genes after calculating run_on ", len(df), "sanity check that all genes have been put back")

        df["start_right_run_on_meta"] = df["end"]                                      #start coordinate of run on
        df["end_right_run_on_meta"] = df["start_right_run_on_meta"] + df["run_on_length"]   #end coordinate of run on
        df["run_on_into_downstream_gene"] = ( df["run_on_length"] >= df['right_int_length'] ) & ( df["run_on_length"] > 0.0 )
        df["end_right_run_on_local"] = np.where(df["run_on_into_downstream_gene"] == True , df['start'].shift(-1), df["end_right_run_on_meta"] ) # run on length when stops at next gene
        df["run_on_from_upstream_gene"] = np.where(df["run_on_into_downstream_gene"].shift(1) == True , True, False)                   # run on from upstream gene
        return df

#Pos or min gtf genes
gtf_plu_or_min = "+"
bedgraph_plu_or_min = "pos"
sample = "HS"

path = "../0.0_files/" + sample + "_files/"
filelist_with_path = glob.glob(path + "*" + bedgraph_plu_or_min + '.BedGraph')
filelist = [os.path.basename(x) for x in glob.glob(path + "*" + bedgraph_plu_or_min + '.BedGraph')]



#gtf of min strand
gtf = "../0.0_files/WS258_named_all.gtf"
#gtf = "../0.0_files/WS258.canonical_geneset.gtf"
d_gtf_chr = get_chrm_dic(gtf, gtf_plu_or_min, "gene", " gene_biotype \"protein_coding\"")

directory = "./gtf/"
if not os.path.exists(directory):
    os.makedirs(directory)

directory = "./run_on_contigs/"
if not os.path.exists(directory):
    os.makedirs(directory)

directory = "./df/"
if not os.path.exists(directory):
    os.makedirs(directory)


df_final_all_chrm_all_files = pd.DataFrame()
for f1 in filelist_with_path:
    print("Working on file", f1, "*"*20)
    #Make bedfile dictionary
    my_col = ["chr", "start", "end","count"]
    df_bed = pd.read_csv(f1, sep="\t", header=None, names=my_col)
    df_bed["length"] = df_bed["end"] - df_bed["start"]
    df_bed = df_bed[ df_bed["count"] > 1 ]
    names=df_bed['chr'].unique().tolist()
    d_bed_chr = {chrom : df_bed.loc[df_bed.chr==chrom] for chrom in names}

    df_final_all_chrm = pd.DataFrame() #Initialize final dataframe
    #Do in a loop to get all chrm
    for chrom in d_gtf_chr.keys():
        df = d_gtf_chr[chrom]
        df = get_right_run_on(df)
        #print(df)
        df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df)


df_final_all_chrm_all_files = df_final_all_chrm_all_files.append(df_final_all_chrm)
print(len((df_final_all_chrm_all_files)))
df_final_all_chrm_all_files = df_final_all_chrm_all_files.sort_values('run_on_length', ascending=False).drop_duplicates("gene_id")
print(len(df_final_all_chrm_all_files))

df = df_final_all_chrm_all_files
#Turn floats to ints for bed files later
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
df["length"] = df["length"].astype(int)
df["start_right_run_on_meta"] = df["start_right_run_on_meta"].astype(int)
df["end_right_run_on_meta"] = df["end_right_run_on_meta"].astype(int)
df["end_right_run_on_local"] = df["end_right_run_on_local"].astype(int)
df_all_run_on = df [ df["run_on_length"] > 100]


df_all_run_on.to_csv("./df/" + sample + "_bedgraph_" + bedgraph_plu_or_min + "_gtf_" + gtf_plu_or_min + "_right_cov_local_df.txt", sep="\t", index=None, quoting=csv.QUOTE_NONE)


#make bed file to look at in IGV
df_all_run_on_bed = df_all_run_on[ [ "chr","start_right_run_on_meta","end_right_run_on_local","gene_id"] ].copy()
df_all_run_on_bed.to_csv("./run_on_contigs/" + sample + "_bedgraph_" + bedgraph_plu_or_min + "_gtf_" + gtf_plu_or_min + "_right_cov_local_contigs.bed", sep="\t", header=None, index=None,quoting=csv.QUOTE_NONE)

df_all_gtf_out = df_all_run_on[ ["chr"] ].copy()   #,"source","type","start","end","dot","strand","dot2","gene_id"]].copy()
df_all_gtf_out["source"] = "Wormbase"
df_all_gtf_out["type"] = "exon"
df_all_gtf_out["start_right_run_on_meta"] = df_all_run_on["start_right_run_on_meta"]
df_all_gtf_out["end_right_run_on_meta"] = df_all_run_on["end_right_run_on_meta"]
df_all_gtf_out["dot"] = "."
df_all_gtf_out["strand"] = gtf_plu_or_min
df_all_gtf_out["dot2"] = "."
df_all_gtf_out["gene_id"] = df_all_run_on["gene_id"]
df_all_gtf_out.to_csv("./gtf/" + sample + "_bedgraph_" + bedgraph_plu_or_min + "_gtf_" + gtf_plu_or_min + "_right_cov_out.gtf", sep="\t", header=None, index=None,quoting=csv.QUOTE_NONE)

print("Number of all genes that have run on ", len(df_all_run_on))









