#!/bin/sh
#The Vanilla (Exact UMI Collapsing) SCUMI Pipeline
#Taihsien Ouyang, Yixuan Guo, Darwin Shen, Ranran Hu

## Prerequisites ##
#SRAToolkit, SAMTools, BEDTools, umitools, R
#Aligner can be Bowtie or BWA

## PATHS ##

ALIGNER_PATH=/bowtie/path
SRATOOLKIT_PATH=/sratoolkit/path
HOMER_PATH=/homer/path

WORKING_PATH=/working/path
SRA_FILES=/working/path/*.sra

## CONFIGURATIONS ##
#Options of reference genome: ERCC92, GRCh38, GRCm38
REF_GENOME=GRCm38
#Reference genome for BWA
REF_GENOME_PATH=/ref_genome/path/GRCm38.fa	
#Number of cores to be used for alignment
NUM_CORES=12

## INITIALIZATION ##

mkdir $WORKING_PATH//fastq
mkdir $WORKING_PATH//umi
mkdir $WORKING_PATH//umitrimmed
mkdir $WORKING_PATH//bam
mkdir $WORKING_PATH//bed

## MAIN PROCESS ##

for f in $SRA_FILES
do
	echo "Trimming $f"
	$SRATOOLKIT_PATH//bin/fastq-dump -Z $f > $f.fastq
  	umitools trim --end 5 $f.fastq NNNNN > $f.umi #Change the number of Ns to specify the length of UMI barcodes
  	$HOMER_PATH//bin/homerTools trim -5 GGG -mis 0 $f.umi #Remark this line if there is no template switching
	
	echo "Aligning $f"
	$ALIGNER_PATH//bowtie -p $NUM_CORES --best --sam $REF_GENOME $f.umi.trimmed > $f.sam
    #$BWA_PATH/bwa mem -t $NUM_CORES $REF_GENOME_FILE $f.umi.trimmed > $f.sam

 	echo "Converting $f to BED"
 	samtools view -Sb $f.sam > $f.bam
	samtools view -bq 1 $f.bam > $f.bam2
	samtools sort $f.bam2 $f.bam.sorted
	samtools index $f.bam.sorted.bam
	#bamToBed -i $f.bam.sorted.bam > $f.bed 				#Non-UMI collapsing	
	umitools rmdup $f.bam.sorted.bam $f.bamrm > $f.bed		#Naive UMI Collapsing
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

echo "==DONE=="
