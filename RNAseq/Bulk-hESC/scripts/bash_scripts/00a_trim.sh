# Greg's trim code
find fastq_raw -name "*_R1_001.fastq.gz" | parallel -j 4 trim_galore --paired -o fastq_trimmed/ {} '{= s/R1/R2/ =}'