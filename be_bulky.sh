#! /bin/bash

myreads=$1;
outputs=$2;
genome=$3;


######## QUALITY_CHECK #########


mkdir -p $outputs/qcoutput;

for fq in $myreads/*.fastq;
do
fastqc -t 12 -o $outputs/qcoutput $myreads/$fq  &> fastqc.log;
done;

######### ALIGNMENTS ##########

mkdir -p $outputs/alignment;
mkdir -p $outputs/alignment/logs;

for file in $myreads/*.fastq;
do
 filename=$(basename $file .fastq);
 hisat2 -x $genome/index/hisat_hg38 -U $file -S $outputs/alignment/${filename}.sam flag --mm &> $outputs/alignment/logs/${filename}.log;

samtools view -@ 12 -bo $outputs/alignment/${filename}.bam $outputs/alignment/${filename}.sam &> $outputs/alignment/logs/${filename}_samtools.log;

done;


######## FEATURE_COUNTS ########


mkdir -p $outputs/featurecounts;

featureCounts -g gene -t exon -a $genome/genome/*.gff -o $outputs/featurecounts/genematrix.txt $outputs/alignment/*.bam &> $outputs/featurecounts/featurecounts.log;

cut -f 1,7- $outputs/featurecounts/genematrix.txt | tail +2 > $outputs/featurecounts/genematrix_rawcounts.txt;



############ DESEQ2 ############

Rscript deseq.R $sampleMap $outputs/featurecounts/genematrix_rawcounts.txt;
