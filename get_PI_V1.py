import pysam
import pandas as pd
import os

import argparse
##拆文件
parser = argparse.ArgumentParser(description="reads length ladder")

parser.add_argument("-input-BMA", "--Input_BAM_file",
                    help="Input the BAM file(such as pol2)",required=True)
parser.add_argument("-o", "--output_PI_txt",
                    help="Output the PI txt filename", type=str)
parser.add_argument("-input_TSS_bed", "--Input_TSS_bed_file",
                    help="Input the TSS bed file",required=True)
parser.add_argument("-input_GeneBody_bed", "--Input_GeneBody_bed_file",
                    help="Input the GeneBody bed file",required=True)
parser.add_argument("-input_GenePred_file", "--Input_GenePred_file",
                    help="Input the GenePred bed file",required=True)

ARGS = parser.parse_args()

inputBam = ARGS.Input_BAM_file

TSS_bed_file = ARGS.Input_TSS_bed_file
GeneBody_bed_file = ARGS.Input_GeneBody_bed_file
GenePred_file = ARGS.Input_GenePred_file

outPI = ARGS.output_PI_txt


data = pysam.AlignmentFile(inputBam,"r")

dict_tss = {}
dict_GeneBody = {}

def get_region_Read_counts(samfile,chrosome,region_start,region_end):
    region_start = int(region_start)
    region_end = int(region_end)
    region_length = abs(region_start-region_end)
    region_readCount = samfile.count(contig=chrosome, start=region_start, stop=region_end)
    region_readCountRatio = region_readCount/region_length
    return region_readCount,region_readCountRatio

dict_TSS_region = {}
dict_GeneBody_region = {}
dict_NM_GeneName = {}

with open(GenePred_file,"r") as f:
    next(f)
    for lines in f:
        line = lines.strip().split()
        NM = line[0]
        geneName = line[10]
        exonCount = line[7]
        txStart = int(line[3])
        txEnd = int(line[4])
        geneLength = abs(txStart - txEnd)
        strand = line[2]
        dict_NM_GeneName[NM] = [geneName,str(exonCount),str(geneLength),strand]

df_NM_GeneName = pd.DataFrame.from_dict(dict_NM_GeneName,orient='index',columns=["GeneName","exonCount","geneLength","strand"])

with open(TSS_bed_file,"r") as f:
    for lines in f:
        line = lines.strip().split()
        Chr = line[0]
        start = int(line[1]) - 1
        end = int(line[2]) - 1 
        NM = line[3]
        dict_TSS_region[NM]=[Chr,str(start),str(end)]
        Tss_readCount,Tss_readCountRatio = get_region_Read_counts(data,Chr,start,end)
        dict_tss[NM] = Tss_readCountRatio
        # print(Chr,start,end)

df_TSS_region = pd.DataFrame.from_dict(dict_TSS_region,orient='index',columns=["chr_TSS","chr_tss_start","chr_tss_end"])
# print(dict_TSS_region)

with open(GeneBody_bed_file,"r") as f:
    for lines in f:
        line = lines.strip().split()
        Chr = line[0]
        start = int(line[1]) - 1
        end = int(line[2]) - 1 
        NM = line[3]
        dict_GeneBody_region[NM]=[Chr,str(start),str(end)]
        try:
            GeneBody_readCount,GeneBody_readCountRatio = get_region_Read_counts(data,Chr,start,end)
            dict_GeneBody[NM] = GeneBody_readCountRatio
        except:
            pass

df_GeneBody_region = pd.DataFrame.from_dict(dict_GeneBody_region,orient='index',columns=["chr_GB","chr_GB_start","chr_GB_end"])
Data_info1 = pd.merge(df_NM_GeneName,df_TSS_region,left_index=True,right_index=True)
Data_info = pd.merge(Data_info1,df_GeneBody_region,left_index=True,right_index=True)


##将dict转成DataFrame
df_TSS = pd.DataFrame.from_dict(dict_tss,orient='index',columns=["TSS_ReadCountRatio"])
df_GeneBody = pd.DataFrame.from_dict(dict_GeneBody,orient='index',columns=["GeneBody_ReadCountRatio"])
##两个DataFrame根据index合并
Data_TSS_GeneBody = pd.merge(df_TSS,df_GeneBody,left_index=True,right_index=True)
Data_TSS_GeneBody["PI"] = Data_TSS_GeneBody["TSS_ReadCountRatio"]/1.0/Data_TSS_GeneBody["GeneBody_ReadCountRatio"]
##将information信息与PI值合并在一起
df_data = pd.merge(Data_info,Data_TSS_GeneBody,left_index=True,right_index=True)

# print(df_data)

df_data.to_csv(outPI,sep="\t")

data.close()



