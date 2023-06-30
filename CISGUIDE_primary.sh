#!/bin/sh
################################################################################################################
#Check for dependencies
################################################################################################################
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
################################################################################################################
#Process options
################################################################################################################
WORKPATH=~
FASTASWITCH="FALSE"
CURRENTTRIMLEN=999999
DEDUPOPT=OFF
StartTime=$(date +%s)
Help()
{
   # Display Help
   echo "Options:"
   echo "h     Print this Help."
   echo "p     Set work path. default: home directory"
   echo "f     Switches to fasta mode if TRUE. default: FALSE."
   echo "d     Sets deduplicate option. OFF= no dup filtering, OPT=dup filtering, UMI=UMI consolidation. default:OFF"
   echo "t     Sets trimming length. Value indicate maximum number of nt to keep. default:999999 (meaning no trimming is performed)."
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
	echo "Fasta mode set"
else 
	echo "Fastq mode set"
fi
################################################################################################################
#Directory creation
################################################################################################################
echo "Looking for directory ${WORKPATH}/input"
if [[ -d "${WORKPATH}/input" ]]
then
	echo "Directory already exists. Emptying..."
	rm -r ${WORKPATH}/input
	mkdir ${WORKPATH}/input
	if [ -z "$(ls -A ${WORKPATH}/input)" ]; then
		echo "Old files removed succesfully"
	else
		echo "Cleanup not succesful. Please restart your device and run again."
		exit 1
	fi
else
	echo "Directory does not exist, creating ${WORKPATH}/input..."
	mkdir ${WORKPATH}/input
fi

echo "Looking for directory ${WORKPATH}/bams"
if [[ -d "${WORKPATH}/bams" ]]
then
	echo "Directory already exists. Emptying..."
	rm -r ${WORKPATH}/bams
	mkdir ${WORKPATH}/bams
	if [ -z "$(ls -A ${WORKPATH}/bams)" ]; then
		echo "Old files removed succesfully"
	else
		echo "Cleanup not succesful. Please restart your device and run again."
		exit 1
	fi
else
	echo "Directory does not exist, creating ${WORKPATH}/bams..."
	mkdir ${WORKPATH}/bams
fi

echo "Looking for directory ${WORKPATH}/misc"
if [[ -d "${WORKPATH}/misc" ]]
then
	echo "Directory already exists. Emptying..."
	rm -r ${WORKPATH}/misc
	mkdir ${WORKPATH}/misc
	if [ -z "$(ls -A ${WORKPATH}/misc)" ]; then
		echo "Old files removed succesfully"
	else
		echo "Cleanup not succesful. Please restart your device and run again."
		exit 1
	fi
else
	echo "Directory does not exist, creating ${WORKPATH}/misc..."
	mkdir ${WORKPATH}/misc
fi

################################################################################################################
#Read Sample information
################################################################################################################

if [[ -f "${WORKPATH}/Sample_information.txt" ]]
then
	echo "Found Sample_information.txt"
	if ! dos2unix < "${WORKPATH}/Sample_information.txt" | cmp - "${WORKPATH}/Sample_information.txt" &> /dev/null
			then
			dos2unix -q ${WORKPATH}/Sample_information.txt
			echo "Warning: did you supply ${WORKPATH}/Sample_information.txt in Windows format? Attempting to convert to Unix format"
			fi
else
	echo "Did not find Sample_information.txt, exiting";
	exit 1
fi

#get a list of refs and create indices
readarray -t LIST_REFS  < <(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" 'FNR>1{print $4}' | sort | uniq)
echo "Processing the following fasta reference files:" ${LIST_REFS[*]}
for j in "${LIST_REFS[@]}"
do
if [[ -f "${WORKPATH}/${j}" ]]
then
	echo "Found ${j}"
	NUM_CHROM=$(grep -c '^[\>][0123456789]' ${WORKPATH}/${j})
	if (( NUM_CHROM > 0 ))
	then
		echo "${NUM_CHROM} chromosome names start with a number, exiting";
		exit 1
	else
		if [[ -f "${WORKPATH}/${j}.bwt.2bit.64" ]] && [[ -f "${WORKPATH}/${j}.0123" ]] && [[ -f "${WORKPATH}/${j}.amb" ]] && [[ -f "${WORKPATH}/${j}.ann" ]] && [[ -f "${WORKPATH}/${j}.pac" ]]
		then
			echo "BWA-mem2 index of ${WORKPATH}/${j} found"
		else
			echo "Creating the index for ${j}"
			bwa-mem2 index ${WORKPATH}/${j}
		fi
	fi
else
	echo "Did not find ${j}, exiting";
	exit 1
fi
done

#get the list of files to process
readarray -t LIST_FILES  < <(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" 'FNR>1{print $3}' | sort | uniq) 
echo "Processing the following files:" ${LIST_FILES[*]}

> ${WORKPATH}/file0.temp | awk -v OFS="\t" -v FS="\t" 'BEGIN {print "Sample", "RunID", "File", "Subject", "Type", "Reads"}' > ${WORKPATH}/input/read_numbers.txt




################################################################################################################
#start analysing files
################################################################################################################
echo "Analysing files in regular fastq mode"
for i in "${LIST_FILES[@]}" 
do
################################################################################################################
#Reading variables
################################################################################################################

echo "Acquiring information from Sample_information.txt ..."
CURRENTR1SUFFIX=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $20}}')
echo "Current R1 suffix: ${CURRENTR1SUFFIX}"
CURRENTR2SUFFIX=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $21}}')
echo "Current R2 suffix: ${CURRENTR2SUFFIX}"
CURRENTSAMPLE=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $1}}')
echo "Current sample: ${CURRENTSAMPLE}"
CURRENTREF=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $4}}')
echo "Current reference: $CURRENTREF"
PRIMERSEQFULL=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $2}}')
PRIMERSEQ="$( echo "$PRIMERSEQFULL" | tr "[:lower:]" "[:upper:]" | sed -e 's#^TCAGACGTGTGCTCTTCCGATCT##' )"
if [ -z "$PRIMERSEQ" ]
then
	echo "Primer sequence not found, moving to the next sample ..."
	continue
else
	echo "Current primer sequence: ${PRIMERSEQ}"
fi
CURRENTRUNID=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $16}}')
echo "Current Run ID: ${CURRENTRUNID}"
FOCUSLOCUS=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $10}}')
echo "Current focus locus: ${FOCUSLOCUS}"
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"



################################################################################################################
#Creating output folder
################################################################################################################
echo "Looking for directory ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}"
if [[ -d "${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}" ]]
then
	echo "Directory already exists. Emptying..."
	#add a check to see whether there is anything inside
	rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/*
	if [ -z "$(ls -A ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID})" ]; then
		echo "Old files removed succesfully"
	else
		echo "Cleanup not succesful. Please restart your device and run again."
		exit 1
	fi
else
	echo "Directory does not exist, creating ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}..."
	mkdir ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
fi
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#Looking for data and unzipping
################################################################################################################
echo "checking for the presence of the sequencing data files..."
if [[ -f "${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz, unzipping ..."
	gunzip -c ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.fastq
else
	echo "Did not find ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz, moving to the next sample";
	rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
	continue
fi
if [[ -f "${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz, unzipping ..."
	gunzip -c ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2.fastq
else
	echo "Did not find ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz, moving to the next sample";
	continue
fi

echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#Trimming
################################################################################################################

if [[ ${FASTASWITCH} = FALSE ]]
then
	echo ">p5" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | sed 's/[*]//g' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	echo ">p5_RC" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | sed 's/[*]//g' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	echo ">p7" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | sed 's/[*]//g' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	echo ">p7_RC" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
	cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | sed 's/[*]//g' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa

	cp ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa ${WORKPATH}/misc/${CURRENTSAMPLE}_${CURRENTRUNID}_Illumina_adapters.fa

	echo "Removing adapters, cropping until ${CURRENTTRIMLEN} bp, and discarding reads shorter than 90 bp"
	echo "###########################################################################"
	trimmomatic PE ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_unpaired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_unpaired.fastq ILLUMINACLIP:${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa:2:30:10:1:TRUE CROP:${CURRENTTRIMLEN} MINLEN:90 -phred33
	echo "###########################################################################"
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"
else
	echo "Cropping until ${CURRENTTRIMLEN} bp, and discarding reads shorter than 90 bp"
	echo "###########################################################################"
	trimmomatic PE ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_unpaired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_unpaired.fastq CROP:${CURRENTTRIMLEN} MINLEN:90 -phred33
	echo "###########################################################################"
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"
fi

################################################################################################################
#Mapping
################################################################################################################

echo "Mapping" ${i}
echo "###########################################################################"
bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_paired.fastq >${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam
echo "###########################################################################"

samtools view -1 ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.bam
samtools sort -o ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.bam
samtools index ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam.bai

echo "Creating empty output files"
> ${WORKPATH}/file1.temp | awk -v OFS="\t" -v FS="\t" ' BEGIN{print "QNAME", "RNAME_1", "POS_1", "CIGAR_1", "SEQ_1", "QUAL_1", "SATAG_1", "SEQ_RCed_1", "RNAME_2", "POS_2", "CIGAR_2", "SEQ_2", "QUAL_2", "SATAG_2", "SEQ_RCed_2", "FILE_NAME", "PRIMER_SEQ", "DEDUP_METHOD", "TRIM_LEN"}' > ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#Filtering reads
################################################################################################################

#filtering reads
#note meaning of the flags:
#0x100		SECONDARY		secondary alignment
#0x4		UNMAP			segment unmapped
#0x800		SUPPLEMENTARY	supplementary alignment
#0x80		READ2			the last segment in the template
#0x8		MUNMAP			next segment in the template unmapped
#0x10		REVERSE			SEQ is reverse complemented
#0x40		READ1			the first segment in the template

#filter steps:
#remove supplementary and secondary alignment reads (this information is already present in the SAtag, and in fact the sup/sec reads do not contain the full info)
#remove unmapped reads
#only keep reads with MAPQ above 50

echo "Processing forward reads of ${i}: discarding unmapped reads or reads with unmapped mates, discarding supplementary or secondary alignment reads, discarding reads with MAPQ <51"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800 | samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $0}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_0.txt
#get the SAtag info
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_0.txt | sed -En 's/.*(SA:Z:.*;)|.*/\1/p' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_1.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_0.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_2.txt
paste -d $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_2.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw_1.txt | awk -v OFS="\t" -v FS="\t" '{print $0, "FALSE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Processing reverse reads of ${i}: discarding unmapped reads or reads with unmapped mates, discarding supplementary or secondary alignment reads, discarding reads with MAPQ <51"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $0}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_0.txt
#get the SA tag info
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_0.txt | sed -En 's/.*(SA:Z:.*;)|.*/\1/p' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_1.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_0.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_2.txt
paste -d $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_2.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_1.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Reverse complementing reverse reads of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

if [[ ${FASTASWITCH} = FALSE ]]
then
	echo "Combining forward and reverse reads and keeping only those starting with primer seq ${PRIMERSEQ}"
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" -v test="$PRIMERSEQ" '$10 ~ "^"test {print}'  > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt
	awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"
else
	echo "Combining forward and reverse reads, ignoring primer seq"
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt | sort -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt
	awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"
fi

echo "Processing forward mates of accepted reads of ${i}: discarding unmapped reads or reads with unmapped mates, discarding supplementary or secondary alignment reads, discarding reads with MAPQ <51"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $0}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_0.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_0.txt | sed -En 's/.*(SA:Z:.*;)|.*/\1/p' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_1.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_0.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_2.txt
paste -d $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_2.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw_1.txt | awk -v OFS="\t" -v FS="\t" '{print $0, "FALSE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Processing reverse mates of accepted reads of ${i}: discarding unmapped reads or reads with unmapped mates, discarding supplementary or secondary alignment reads, discarding reads with MAPQ <51"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $0}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_0.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_0.txt | sed -En 's/.*(SA:Z:.*;)|.*/\1/p' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_1.txt
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_0.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_2.txt
paste -d $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_2.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_1.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Reverse complementing reverse mates of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Combining forward and reverse mates of ${i}"
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt | sort -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Combining selected reads and mates of ${i}, and write mapping info plus options to output"
join -j 1 -o 1.1,1.3,1.4,1.6,1.10,1.11,1.12,1.13,2.3,2.4,2.6,2.10,2.11,2.12,2.13 -t $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt | 
awk -v OFS="\t" -v FS="\t" -v i="$i" -v PRIMERSEQ="$PRIMERSEQ" -v CURRENTTRIMLEN="$CURRENTTRIMLEN" ' {print $0, i, PRIMERSEQ, CURRENTTRIMLEN}'  >> ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"



################################################################################################################
#Statistics
################################################################################################################

echo "Counting reads of ${CURRENTSAMPLE} ${CURRENTRUNID}"
RAWNO=$(( $(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq | wc -l)/4 )) 
echo "RAWNO: ${RAWNO}"
cat ${WORKPATH}/input/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v RAWNO="${RAWNO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Raw", RAWNO}' >> ${WORKPATH}/input/read_numbers.txt

MAPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam | grep -Ev '^(\@)' | awk '$3 != "*" {print $0}' | sort -u -t$'\t' -k1,1 | wc -l)
echo "MAPNO: ${MAPNO}"
cat ${WORKPATH}/input/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v MAPNO="${MAPNO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Mapped", MAPNO}' >> ${WORKPATH}/input/read_numbers.txt

PREPRONO=$(( $(cat ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt | wc -l)-1 ))
echo "PREPRONO: ${PREPRONO}"
cat ${WORKPATH}/input/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v PREPRONO="${PREPRONO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Filtered", PREPRONO}' >> ${WORKPATH}/input/read_numbers.txt


echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#Cleanup
################################################################################################################

echo "Removing temporary files"

rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
rm ${WORKPATH}/file1.temp 
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

done

echo "Removing temporary files"
rm ${WORKPATH}/file0.temp 

echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

exit 0

