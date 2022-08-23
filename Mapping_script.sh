#!/bin/sh

FOCUS_PLASMID="pCAS-PPO"

if [ -f "shared/"${FOCUS_PLASMID}".fa" ]
then
	echo "Found "${FOCUS_PLASMID}".fa"
else
	echo "Did not find shared/"${FOCUS_PLASMID}", exiting";
	exit 1
fi

if [ -f "shared/Illumina_adapters.fa" ]
then
	echo "Found Illumina_adapters.fa"
else
	echo "Did not find Illumina_adapters.fa, exiting";
	exit 1
fi

echo "Creating the index"
bwa-mem2 index shared/${FOCUS_PLASMID}.fa

LIST_SAMPLES=($(ls shared/*_R1.fastq.gz | tr '\n' ' ' | sed 's/_R1.fastq.gz//g' | sed 's/shared\///g' | sort | uniq))
echo "Processing the following samples:" ${LIST_SAMPLES[*]}




for i in "${LIST_SAMPLES[@]}"  
do
if [ -f "shared/${i}_R1.fastq.gz" ]
then
	echo "Found shared/${i}_R1.fastq.gz"
else
	echo "Did not find shared/${i}_R1.fastq.gz, moving to the next sample";
	continue
fi
if [ -f "shared/${i}_R2.fastq.gz" ]
then
	echo "Found shared/${i}_R2.fastq.gz"
else
	echo "Did not find shared/${i}_R2.fastq.gz, moving to the next sample";
	continue
fi

echo "creating directory shared/${i}"
mkdir shared/${i}

echo "Trimming" ${i}
trimmomatic PE shared/${i}_R1.fastq.gz shared/${i}_R2.fastq.gz shared/${i}/${i}_forward_paired.fastq.gz shared/${i}/${i}_forward_unpaired.fastq.gz shared/${i}/${i}_reverse_paired.fastq.gz shared/${i}/${i}_reverse_unpaired.fastq.gz ILLUMINACLIP:shared/Illumina_adapters.fa:2:30:10:1:TRUE CROP:150 -phred33

echo "Unzipping" ${i}
gunzip shared/${i}/${i}_forward_paired.fastq.gz 
gunzip shared/${i}/${i}_reverse_paired.fastq.gz 

echo "Mapping" ${i}
bwa-mem2 mem shared/${FOCUS_PLASMID}.fa shared/${i}/${i}_forward_paired.fastq shared/${i}/${i}_reverse_paired.fastq >shared/${i}/${i}.sam

samtools view -1 shared/${i}/${i}.sam > shared/${i}/${i}.bam
samtools sort -o shared/${i}/${i}_sorted.bam shared/${i}/${i}.bam
echo "Dedupping" ${i}
picard MarkDuplicates --INPUT shared/${i}/${i}_sorted.bam --OUTPUT shared/${i}/${i}.sorted.bam --METRICS_FILE shared/${i}/${i}_sorted_dedup_metrics.txt --REMOVE_DUPLICATES TRUE
samtools index shared/${i}/${i}.sorted.bam shared/${i}/${i}.sorted.bam.bai


done

