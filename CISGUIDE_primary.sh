#!/bin/sh
#first check whether dependent programs have been installed
if ! command -v picard &> /dev/null
then
    echo "Please first install picardTools before running the CISGUIDE program"
    exit 1
fi
if ! command -v bwa-mem2 &> /dev/null
then
    echo "Please first install bwa-mem2 before running the CISGUIDE program"
    exit 1
fi
if ! command -v samtools &> /dev/null
then
    echo "Please first install samtools before running the CISGUIDE program"
    exit 1
fi
if ! command -v trimmomatic &> /dev/null
then
    echo "Please first install trimmomatic before running the CISGUIDE program"
    exit 1
fi
if ! command -v dos2unix &> /dev/null
then
    echo "Please first install dos2unix before running the CISGUIDE program"
    exit 1
fi

#set variables
WORKPATH=~
FASTASWITCH="FALSE"

#Process Options
Help()
{
   # Display Help
   echo "Options:"
   echo "h     Print this Help."
   echo "p     Set work path. default: home directory"
   echo "f     Switches to fasta mode if TRUE. default: FALSE."
   echo "d     Sets deduplicate option. OFF= no dup filtering, OPT=optical dup filtering, UMI=UMI consolidation. default:OFF"
   echo "t     Sets trimming length. Value indicate maximum number of nt to keep."
   echo
}
while getopts "hp:f:d:k:t:" option; do
   case $option in
      h) # display Help
         Help
         exit 1
		 ;;
      p) # Enter a workdir
         WORKPATH=$OPTARG
		 ;;
      f) # Enter TRUE or FALSE
         FASTASWITCH=$OPTARG
		 ;;
      d) #duplicate filtering options
         DEDUPOPT=$OPTARG
		 ;;
      t) #duplicate filtering options
         CURRENTTRIMLEN=$OPTARG
		 ;;
     \?) # Invalid option
         echo "Error: Invalid option. Exiting."
         exit 1
		 ;;
   esac
done

#check for correct work path
if [[ -d "$WORKPATH" ]]; 
then
	echo "looking in ${WORKPATH}"
else 
	echo "invalid directory ${WORKPATH}. Exiting."
	exit 1
fi

#check whether trimlength is a number
if [[ $CURRENTTRIMLEN =~ ^-?[0-9]+$ ]]
then
	echo "Trim length set at ${CURRENTTRIMLEN}."
else
	echo "Trim length ${CURRENTTRIMLEN} invalid, must be an integer. Exiting."
	exit 1
fi

#check whether fasta switch is acivated
if [[ ${FASTASWITCH} = TRUE ]]
then
	echo "fasta mode set"
else 
	echo "fastq mode set"
fi

#read the sample_information file
if [[ -f "${WORKPATH}/Sample_information.txt" ]]
then
	echo "Found Sample_information.txt"
else
	echo "Did not find Sample_information.txt, exiting";
	exit 1
fi

#get a list of refs and create indices
#REF
readarray -t LIST_REFS  < <(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" 'FNR>1{print $4}' | sort | uniq)
echo "Processing the following fasta reference files:" ${LIST_REFS[*]}
for j in "${LIST_REFS[@]}"
do
if [[ -f "${WORKPATH}/${j}" ]]
then
	echo "Found ${j}"
	if [[ -f "${WORKPATH}/${j}.bwt.2bit.64" ]] && [[ -f "${WORKPATH}/${j}.0123" ]] && [[ -f "${WORKPATH}/${j}.amb" ]] && [[ -f "${WORKPATH}/${j}.ann" ]] && [[ -f "${WORKPATH}/${j}.pac" ]]
	then
		echo "BWA-mem2 index of ${WORKPATH}/${j} found"
	else
		echo "Creating the index for ${j}"
		bwa-mem2 index ${WORKPATH}/${j}
	fi
else
	echo "Did not find ${j}, exiting";
	exit 1
fi
done

#get the list of samples
#fastq_name
readarray -t LIST_SAMPLES  < <(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" 'FNR>1{print $3}' | sort | uniq) 
echo "Processing the following samples:" ${LIST_SAMPLES[*]}

> ${WORKPATH}/file0.temp | awk -v OFS="\t" -v FS="\t" 'BEGIN {print "Sample", "Raw read count", "mapped count", "dedupped count", "preprocessed count"}' > ${WORKPATH}/read_numbers.txt


#process all samples

if [[ ${FASTASWITCH} = TRUE ]]
then
	echo "analysing samples in fasta mode";
	for i in "${LIST_SAMPLES[@]}" 
	do
	CURRENTSAMPLE=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $1}}')
	echo "Analyzing sample: ${CURRENTSAMPLE}"
	echo "looking for directory ${WORKPATH}/${CURRENTSAMPLE}"
	if [[ -d "${WORKPATH}/${CURRENTSAMPLE}" ]]
	then
		echo "replacing existing directory"
		rm -rf ${WORKPATH}/${CURRENTSAMPLE}
		mkdir ${WORKPATH}/${CURRENTSAMPLE}
		mkdir ${WORKPATH}/${CURRENTSAMPLE}/temp
	else
		echo "directory does not exist, creating ${WORKPATH}/${CURRENTSAMPLE}"
		mkdir ${WORKPATH}/${CURRENTSAMPLE}
		mkdir ${WORKPATH}/${CURRENTSAMPLE}/temp
	fi
		if [[ -f "${WORKPATH}/${i}.fasta" ]]
		then
			echo "Found ${WORKPATH}/${i}.fasta"
			if ! dos2unix < "${WORKPATH}/${i}.fasta" | cmp - "${WORKPATH}/${i}.fasta" &> /dev/null
			then
			dos2unix -q ${WORKPATH}/${i}.fasta
			echo "Warning: did you supply ${WORKPATH}/${i}.fasta in Windows format? Attempting to convert to Unix format"
			fi
		else
			echo "Did not find ${WORKPATH}/${i}.fasta, moving to the next sample";
			continue
		fi
		LineTotal=$(wc -l < shared/pBas03416_inserts.fasta)
		declare -i LineCounter=1
		SeqVar=""
		cat shared/pBas03416_inserts.fasta |
		while read -r line
		do
			if [[ $line == '>'* ]]
			then 
				echo "$SeqVar" 
				SeqVar="" 
				echo "$line"  
				LineCounter=$(expr $LineCounter + 1)
			else
				
				if [[ $LineCounter -lt $LineTotal ]]
				then
					SeqVar+=$line 
					LineCounter=$(expr $LineCounter + 1)
				else
					echo "$SeqVar" | sed 's/\n//g'
					SeqVar="";
				fi
			fi 
		done |
			sed '/./,$!d' |
			sed 'N;s/\n/\t/g' |
			awk -v OFS="\t" -v FS="\t" '{sub(/[[:space:]].*$/,"",$1); print $1, $2, "+", $2}' | 
			awk -v OFS="\t" -v FS="\t" '{gsub(/[[:alpha:]]/,"G",$4); print}' |
			sed 's/\t/\n/g' | 
			sed 's/>/@/g' > shared/pBas03416_inserts/pBas03416_inserts.fa
		#get the ref and perform mapping
		CURRENTREF=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $4}}')
		echo "Using ref: $CURRENTREF"
		echo "Mapping" ${i}
		bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}/$i.fa > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam
		samtools view -1 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.bam
		samtools sort -o ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.bam
		cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam
		samtools index ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam.bai

	done
else
	echo "analysing samples in regular fastq mode"
for i in "${LIST_SAMPLES[@]}" 
do
#check for presence of R1 and R2
if [[ -f "${WORKPATH}/${i}_R1.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}_R1.fastq.gz"
else
	echo "Did not find ${WORKPATH}/${i}_R1.fastq.gz, moving to the next sample";
	continue
fi
if [[ -f "${WORKPATH}/${i}_R2.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}_R2.fastq.gz"
else
	echo "Did not find ${WORKPATH}/${i}_R2.fastq.gz, moving to the next sample";
	continue
fi

#fastq_name, sample
CURRENTSAMPLE=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $1}}')
echo "Using sample: ${CURRENTSAMPLE}"
echo "looking for directory ${WORKPATH}/${CURRENTSAMPLE}"
if [[ -d "${WORKPATH}/${CURRENTSAMPLE}" ]]
then
	echo "replacing existing directory"
	rm -rf ${WORKPATH}/${CURRENTSAMPLE}
	mkdir ${WORKPATH}/${CURRENTSAMPLE}
	mkdir ${WORKPATH}/${CURRENTSAMPLE}/temp
else
	echo "directory does not exist, creating ${WORKPATH}/${CURRENTSAMPLE}"
	mkdir ${WORKPATH}/${CURRENTSAMPLE}
	mkdir ${WORKPATH}/${CURRENTSAMPLE}/temp
fi
echo ">p5" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
echo ">p5_RC" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
echo ">p7" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
echo ">p7_RC" >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa

echo "Trimming" ${i}
trimmomatic PE ${WORKPATH}/${i}_R1.fastq.gz ${WORKPATH}/${i}_R2.fastq.gz ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_paired.fastq.gz ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_unpaired.fastq.gz ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_paired.fastq.gz ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_unpaired.fastq.gz ILLUMINACLIP:${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa:2:30:10:1:TRUE CROP:${CURRENTTRIMLEN} -phred33

echo "Unzipping" ${i}
gunzip ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_paired.fastq.gz 
gunzip ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_paired.fastq.gz 

#fastqname, ref
CURRENTREF=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $4}}')
echo "Using ref: $CURRENTREF"
echo "Mapping" ${i}
bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_paired.fastq >${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam

samtools view -1 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.bam
samtools sort -o ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.bam

#dedupping
case $DEDUPOPT in
	OPT) #optical duplicate filtering
		echo "Optical duplicate filtering of ${i}"
		picard MarkDuplicates --INPUT ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam --OUTPUT ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam --METRICS_FILE ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted_dedup_metrics.txt --REMOVE_DUPLICATES TRUE
		samtools index ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam.bai
		;;
	UMI) #UMI consolidation
		#echo "UMI consolidation of ${i}"
		echo "UMI consolidation does not work yet, skipping filtering ${i}"
		cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam
		samtools index ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam.bai
		;;
	OFF) #no duplicate filtering
		echo "Skipping duplicate filtering of ${i}"
		cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam
		samtools index ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam.bai
		;;
	*) #default, no duplicate filtering
		echo "Skipping duplicate filtering of ${i}"
		cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam
		samtools index ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam.bai
		;;
esac

#getting primer sequence
PRIMERSEQFULL=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $2}}')
PRIMERSEQ="$( echo "$PRIMERSEQFULL" | tr "[:lower:]" "[:upper:]" | sed -e 's#^TCAGACGTGTGCTCTTCCGATCT##' )"

if [ -z "$PRIMERSEQ" ]
then
	echo "Primer sequence not found, moving to the next sample"
	continue
else
	echo "Using this primer sequence: ${PRIMERSEQ}"
fi


echo "Creating empty output files"
> ${WORKPATH}/file1.temp | awk -v OFS="\t" -v FS="\t" ' BEGIN{print "QNAME", "RNAME_1", "POS_1", "CIGAR_1", "SEQ_1", "QUAL_1", "SATAG_1", "SEQ_RCed_1", "RNAME_2", "POS_2", "CIGAR_2", "SEQ_2", "QUAL_2", "SATAG_2", "SEQ_RCed_2", "FILE_NAME", "PRIMER_SEQ"}' > ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_A.txt

echo "Processing fw reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800 | samtools view -uf 0x80 |samtools view -uF 0x8 | samtools view -F 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_fw.txt

echo "Processing rev reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x80 |samtools view -uF 0x8 | samtools view -f 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv.txt

echo "Rev complementing rev reads of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv.txt ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCed.txt

echo "Combining fw and rev reads and keeping only those starting with primer seq ${PRIMERSEQ}"
cat ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_fw.txt ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCed.txt | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" -v test="$PRIMERSEQ" '$10 ~ "^"test {print}'  > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all.txt
awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all.txt > ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all_names.txt

echo "Processing fw mates of accepted reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}/Mates_fw.txt

echo "Processing rev mates of accepted reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv.txt

echo "Rev complementing rev mates of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv.txt ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCed.txt

echo "Combining fw and rev mates of ${i}"
cat ${WORKPATH}/${CURRENTSAMPLE}/Mates_fw.txt ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCed.txt | sort -T ${WORKPATH}/${CURRENTSAMPLE}/temp/ -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}/Mates_all.txt

echo "Combining selected reads and mates of ${i}"
join -j 1 -o 1.1,1.3,1.4,1.6,1.10,1.11,1.12,1.13,2.3,2.4,2.6,2.10,2.11,2.12,2.13 -t $'\t' ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all.txt ${WORKPATH}/${CURRENTSAMPLE}/Mates_all.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" -v PRIMERSEQ="$PRIMERSEQ" ' {print $0, i, PRIMERSEQ}'  >> ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_A.txt

echo "counting reads of ${CURRENTSAMPLE}"
RAWNO=$(gunzip -c ${WORKPATH}/${i}_R1.fastq.gz | wc -l)
echo "RAWNO: ${RAWNO}"
MAPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam | grep -Ev '^(\@)' | awk '$3 != "*" {print $0}' | sort -u -t$'\t' -k1,1 | wc -l)
echo "MAPNO: ${MAPNO}"
samtools view -h -o ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.dedup.sam ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sorted.bam
DEDUPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.dedup.sam | grep -Ev '^(\@)' | awk '$3 != "*" {print $0}' | sort -u -t$'\t' -k1,1 | wc -l)
echo "DEDUPNO: ${DEDUPNO}"
PREPRONO=$(cat ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_A.txt | wc -l)
echo "PREPRONO: ${PREPRONO}"
cat ${WORKPATH}/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v RAWNO="${RAWNO}" -v MAPNO="${MAPNO}" -v DEDUPNO="${DEDUPNO}" -v PREPRONO="${PREPRONO}" ' END{print CURRENTSAMPLE, RAWNO / 4, MAPNO, DEDUPNO, PREPRONO - 1}' >> ${WORKPATH}/read_numbers.txt

echo "removing temporary files"
rm ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.bam 
rm ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted.bam 
rm ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.sam 

if [[ $DEDUPOPT = OPT ]]
then
	rm ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}.dedup.sam 
fi
rm ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCseqs.txt
rm ${WORKPATH}/${CURRENTSAMPLE}/Mates_all.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Mates_fw.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Mates_rv_RCed.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_fw.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCed.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_all_names.txt 
rm ${WORKPATH}/${CURRENTSAMPLE}/Primer_reads_rv_RCseqs.txt
rm ${WORKPATH}/${CURRENTSAMPLE}/Illumina_adapters.fa
rm ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_paired.fastq 
rm ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_paired.fastq 
rm ${WORKPATH}/${CURRENTSAMPLE}/${i}_forward_unpaired.fastq.gz 
rm ${WORKPATH}/${CURRENTSAMPLE}/${i}_reverse_unpaired.fastq.gz 
if [[ $DEDUPOPT = OPT ]]
then
	rm ${WORKPATH}/${CURRENTSAMPLE}/${CURRENTSAMPLE}_sorted_dedup_metrics.txt 
fi
rm -r ${WORKPATH}/${CURRENTSAMPLE}/temp
rm ${WORKPATH}/file1.temp 

done

fi

echo "removing temporary files"
rm ${WORKPATH}/file0.temp 

now=$(date)
echo "Finished at $now"

exit 0
