# check trimmed reads
mkdir fastqc_trimmed
fastqc fastq_trimmed/*.fq.gz --outdir fastqc_trimmed
multiqc fastqc_trimmed/ -n multiqc_trimmed