#!/bin/bash/

for fq in fastq_trimmed/*_R1_001_val_1.fq.gz
do 
    sample=`basename $fq _R1_001_val_1.fq.gz`
    r1=$sample"_R1_001_val_1.fq.gz"
    r2=$sample"_R2_001_val_2.fq.gz"
    salmon quant -i ref_transcriptome/salmon_index -l A -1 fastq_trimmed/$r1 -2 fastq_trimmed/$r2 --seqBias --gcBias --posBias -o tx_quantification/$sample
done