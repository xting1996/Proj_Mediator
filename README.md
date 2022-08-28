# Proj_Mediator
Codes of Mediator Project


#### 1. Pausing Index count
```
GeneBody=/my_path/Gene_BODY_1218_genelengthEND.V1.bed
GeneTSS=/my_path/Gene_TSS_1218_genelengthEND.V1.bed
GenePred=/my_path/hg19_genePred_1216
getPI=/my_path/get_PI_V1.py

python ${getPI} -input_TSS_bed ${GeneTSS} -input_GeneBody_bed ${GeneBody} -input_GenePred_file ${GenePred} -input-BMA 5_pol2_FKDL210360537-1a.date1017.sort.rmdup.q30.bam -o 5_pol2_FKDL210360537-1a.date1017.sort.rmdup.q30.PI.txt
```
