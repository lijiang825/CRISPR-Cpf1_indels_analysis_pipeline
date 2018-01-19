#!/bin/bash

if [ "$domain"=="Orchestra" ]; 
then
    module load seq/fastx/0.0.13
    module load seq/bwa/0.7.8
    module load seq/samtools/1.3
    module load seq/emboss/6.6.0
    module load dev/python/2.7.10
    module load dev/java/jdk1.8
fi

# Full path to GFP index.
INDEX_DIR="/home/lj63/data/cpf1_mutation/20171006_Cpf1_GFP_DR_mutant_R49_U937/ref_index/Cpf1_DR_mutant_GFP.fas"

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR="/home/lj63/data/cpf1_mutation/20171006_Cpf1_GFP_DR_mutant_R49_U937/fastq"

# Full path to output folder.
OUTPUT_DIR="/home/lj63/data/cpf1_mutation/20171006_Cpf1_GFP_DR_mutant_R49_U937/analysis/"

# BWA pipeline 
for file in $FASTQ_DIR/*.fastq
do
   NAME=$(basename $file .fastq)
   echo $NAME
   fastq_quality_filter -v -q 25 -p 95 -Q33 -i $file -o $OUTPUT_DIR/$NAME.filtered.fastq
   bwa mem -O 10 -E 1 -t 4 -R '@RG\tID:20171006\tPL:illumina\tPU:unit1\tLB:${NAME}\tSM:${NAME}' $INDEX_DIR $FASTQ_DIR/$NAME.fastq > $OUTPUT_DIR/$NAME.sam    

   grep -v "@" $OUTPUT_DIR/$NAME.sam | awk -F"\t" 'BEGIN{print "flag\toccurrences"} {a[$2]++} END{for(i in a)print i"\t"a[i]}' > $OUTPUT_DIR/$NAME.sum

   samtools view -Shu $OUTPUT_DIR/$NAME.sam | samtools sort - -o $OUTPUT_DIR/$NAME.bam 
   samtools index $OUTPUT_DIR/$NAME.bam
   
done

for file in $FASTQ_DIR/*.fastq
do
    NAME=$(basename $file .fastq)
    echo $NAME   
    python VariantsCalling.py $OUTPUT_DIR/$NAME.sam 
done



