#!/bin/sh
#The Vanilla (Exact UMI Collapsing) SCUMI Pipeline
#Taihsien Ouyang, Yixuan Guo, Darwin Shen, Ranran Hu

## Prerequisites ##
#SRAToolkit, SAMTools, BEDTools, umitools, CTK, CZPlib, R
#Aligner can be Bowtie or BWA

## PATHS ##

ALIGNER_PATH=/bowtie/path
SRATOOLKIT_PATH=/sratoolkit/path
HOMER_PATH=/homer/path
CTK_PATH=/CTK/path

WORKING_PATH=/working/path
SRA_FILES=/working/path/*.sra

## CONFIGURATIONS ##
#Options of reference genome: ERCC92, GRCh38, GRCm38
REF_GENOME=GRCm38 	
#Number of cores to be used for alignment
NUM_CORES=12

## INITIALIZATION ##

PERL5LIB=/usr/local/lib/czplib.v1.0.6/

mkdir $WORKING_PATH//fastq
mkdir $WORKING_PATH//umi
mkdir $WORKING_PATH//umitrimmed
mkdir $WORKING_PATH//bam
mkdir $WORKING_PATH//bed
mkdir $WORKING_PATH//ctk

## MAIN PROCESS ##

for f in $SRA_FILES
do
	echo "Trimming $f"
	$SRATOOLKIT_PATH//bin/fastq-dump -Z $f > $f.fastq
  	umitools trim --end 5 $f.fastq NNNNN > $f.umi 			#Change the number of Ns to specify the length of UMI barcodes
  	$HOMER_PATH//bin/homerTools trim -5 GGG -mis 0 $f.umi 	#Remark this line if there is no template switching
	
	echo "Aligning $f"
	$ALIGNER_PATH//bowtie -p $NUM_CORES --best --sam $REF_GENOME $f.umi.trimmed > $f.sam

 	echo "Converting $f to CTK-BED"
 	samtools view -Sb $f.sam > $f.bam
	samtools view -bq 1 $f.bam > $f.bam2
	samtools sort $f.bam2 $f.bam.sorted
	samtools index $f.bam.sorted.bam
	bamToBed -i $f.bam.sorted.bam > $f.bed 
	sed -i 's/:UMI_/#1#/g' $f.bed 							#Converting the format of read tag for CTK
	perl $CTK_PATH//tag2collapse.pl $f.bed $f.bed.ctk
	echo "Done with $f"

done

mv $WORKING_PATH//*.fastq $WORKING_PATH//fastq
mv $WORKING_PATH//*.umi $WORKING_PATH//umi
mv $WORKING_PATH//*.sam $WORKING_PATH//bam
mv $WORKING_PATH//*.bam $WORKING_PATH//bam
mv $WORKING_PATH//*.bai $WORKING_PATH//bam
mv $WORKING_PATH//*.bam2 $WORKING_PATH//bam
mv $WORKING_PATH//*.umi.trimmed $WORKING_PATH//umitrimmed
mv $WORKING_PATH//*.lengths $WORKING_PATH//umitrimmed
mv $WORKING_PATH//*.bed $WORKING_PATH//bed
mv $WORKING_PATH//*.bed.ctk $WORKING_PATH//ctk

echo "==DONE=="
