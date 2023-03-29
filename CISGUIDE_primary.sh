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
StartTime=$(date +%s)
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

#check the duplicate filtering option
case $DEDUPOPT in
	OPT) 
	echo "Program set to filter optical duplicates"
	;;
	UMI) 
	echo "Program set to remove duplicates by UMI consolidation"
	;;
	OFF) 
	echo "Duplicate removal switched off"
	;;
	*) 
	echo "Duplicate removal switched off"
	;;
esac

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

echo "Looking for directory ${WORKPATH}/stats"
if [[ -d "${WORKPATH}/stats" ]]
then
	echo "Directory already exists. Emptying..."
	rm -r ${WORKPATH}/stats
	mkdir ${WORKPATH}/stats
	if [ -z "$(ls -A ${WORKPATH}/stats)" ]; then
		echo "Old files removed succesfully"
	else
		echo "Cleanup not succesful. Please restart your device and run again."
		exit 1
	fi
else
	echo "Directory does not exist, creating ${WORKPATH}/stats..."
	mkdir ${WORKPATH}/stats
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

#get the list of files to process
readarray -t LIST_FILES  < <(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" 'FNR>1{print $3}' | sort | uniq) 
echo "Processing the following files:" ${LIST_FILES[*]}

> ${WORKPATH}/file0.temp | awk -v OFS="\t" -v FS="\t" 'BEGIN {print "Sample", "RunID", "File", "Subject", "Type", "Reads"}' > ${WORKPATH}/stats/read_numbers.txt

################################################################################################################
#FASTA mode: read samples, check refs and directories
################################################################################################################

if [[ ${FASTASWITCH} = TRUE ]]
then
	echo "Analysing samples in fasta mode";
	for i in "${LIST_FILES[@]}" 
	do
	
	CURRENTSAMPLE=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $1}}')
	echo "Analyzing sample: ${CURRENTSAMPLE}"
	echo "Looking for directory ${WORKPATH}/${CURRENTSAMPLE}"
	if [[ -d "${WORKPATH}/${CURRENTSAMPLE}" ]]
	then
		echo "Directory already exists. Emptying."
		rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/*
	else
		echo "Directory does not exist, creating ${WORKPATH}/${CURRENTSAMPLE}"
		mkdir ${WORKPATH}/${CURRENTSAMPLE}
		mkdir ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp
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
	CURRENTREF=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $4}}')
	echo "Using ref: ${CURRENTREF} file"
	CURRENTPLASMID=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $7}}')
	echo "Plasmid name: ${CURRENTPLASMID}"
	CURRENTLBSEQ=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $17}}')
	echo "Current LB sequence: ${CURRENTLBSEQ}"
	CURRENTRBSEQ=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $18}}')
	echo "Current RB sequence: ${CURRENTRBSEQ}"
	
	CURRENTRBPOS=1
	CURRENTLBPOS=749
	
	
################################################################################################################
#FASTA mode: processing the fasta data
################################################################################################################
	
	LineTotal=$(wc -l < ${WORKPATH}/${i}.fasta)
	declare -i LineCounter=1
	SeqVar=""
	cat ${WORKPATH}/${i}.fasta |
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
	sed 's/>/@/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_raw.fa
			
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_raw.fa |
	sort | 
	uniq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa
			
################################################################################################################
#FASTA mode: Premapping first and last 30bp
################################################################################################################
		
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, 1, 30), $3, substr($4, 1, 30)}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_first30.fa
		
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, length($2)-30, length($2)), $3, substr($4, length($4)-30, length($4))}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_last30.fa
		
	#transform the data so that it looks like a fastq file
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_first30.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1" 1:N:0:NNNNNNNN+NNNNNNNN", $2, $3, $4}' |
	sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_first30.fastq
		
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_last30.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1" 1:N:0:NNNNNNNN+NNNNNNNN", $2, $3, $4}' |
	sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_last30.fastq
		
	echo "Premapping 1" ${i}
	echo "###########################################################################"
	bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_first30.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_first30.sam
	bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_last30.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_last30.sam	
	
	#of the "first 30bp" sequences, split the unambiguous sequences into LB and RB files
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_first30.sam | 
	sed -n -E '/^[[:alpha:]]/p' |
	awk -v OFS="\t" -v FS="\t" -v CP="$CURRENTPLASMID" '{if ($3 == CP && $5 > 0) {print $0}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0.txt |
	awk -v OFS="\t" -v FS="\t" -v CLP="$CURRENTLBPOS" -v CRP="$CURRENTRBPOS" '{if ( sqrt(($4-CLP)^2) < sqrt(($4-CRP)^2)) { print $1 }}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_LB_names.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0.txt |
	awk -v OFS="\t" -v FS="\t" -v CLP="$CURRENTLBPOS" -v CRP="$CURRENTRBPOS" '{if ( sqrt(($4-CLP)^2) > sqrt(($4-CRP)^2)) { print $1 }}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_RB_names.txt
	
	#of the "last 30bp" sequences, split the unambiguous sequences into LB and RB files
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_last30.sam | 
	sed -n -E '/^[[:alpha:]]/p' |
	awk -v OFS="\t" -v FS="\t" -v CP="$CURRENTPLASMID" '{if ($3 == CP && $5 > 0) {print $0}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0.txt |
	awk -v OFS="\t" -v FS="\t" -v CLP="$CURRENTLBPOS" -v CRP="$CURRENTRBPOS" '{if ( sqrt(($4-CLP)^2) < sqrt(($4-CRP)^2)) { print $1 }}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_LB_names.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0.txt |
	awk -v OFS="\t" -v FS="\t" -v CLP="$CURRENTLBPOS" -v CRP="$CURRENTRBPOS" '{if ( sqrt(($4-CLP)^2) > sqrt(($4-CRP)^2)) { print $1 }}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_RB_names.txt
	
	
	#of the "first 30bp" sequences, also list the ambiguous sequences
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_first30.sam | 
	sed -n -E '/^[[:alpha:]]/p' |
	awk -v OFS="\t" -v FS="\t" -v CP="$CURRENTPLASMID" '{if ($3 == CP) {print $1, $5}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names.txt |
	awk -v OFS="\t" -v FS="\t" '{if ( $2 == 0 ) {print $1}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQ0.txt
	
	#of the "last 30bp" sequences, also list the ambiguous sequences
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_last30.sam | 
	sed -n -E '/^[[:alpha:]]/p' |
	awk -v OFS="\t" -v FS="\t" -v CP="$CURRENTPLASMID" '{if ($3 == CP) {print $1, $5}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names.txt
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names.txt |
	awk -v OFS="\t" -v FS="\t" '{if ( $2 == 0 ) {print $1}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQ0.txt
	
	
################################################################################################################
#FASTA mode: Premapping first and last 60bp for ambiguous 30bp mappings
################################################################################################################
	
	#for the first 30bp sequences, if there were ambiguous reads, map again but this time 60bp
	if [ -s ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQ0.txt ]
	then
		#first get the list of those that were ambiguous
		readarray -t lines < ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQ0.txt
		for r in "${lines[@]}"
		do
			cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | egrep "${r}"  >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_startT_MQ0.fa
		done	
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_startT_MQ0.fa |
		awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, 1, 60), $3, substr($4, 1, 60)}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_first60.fa
		
		#transform the data so that it looks like a fastq file
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_first60.fa | 
		awk -v OFS="\t" -v FS="\t" '{print $1" 1:N:0:NNNNNNNN+NNNNNNNN", $2, $3, $4}' |
		sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_first60.fastq
	
		#the mapping
		echo "Premapping 2.a" ${i}
		echo "###########################################################################"
		bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_first60.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_first60.sam
	
		#!!! note: change the code below. I need to split that for LB and RB sequences. 
		#acquiring the names of those that end up getting a larger MQ
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_first60.sam | 
		sed -n -E '/^[[:alpha:]]/p' |
		awk -v OFS="\t" -v FS="\t" '{if ($5 > 0) {print $1}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0_60.txt
		
		#then combining the good ones from this set with the good ones from the 30bp set
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0_60.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0_comb.txt
	else
		#also make a comb file even if all were alright in the 30bp set
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0_comb.txt
	fi
	
	#for the last 30bp sequences, if there were ambiguous reads, map again 60bp
	if [ -s ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQ0.txt ]
	then
		#first get the list of those that were ambiguous
		readarray -t lines < ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQ0.txt
		for r in "${lines[@]}"
		do
			cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | egrep ${r}  >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_endT_MQ0.fa
		done	
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_endT_MQ0.fa | 
		awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, length($2)-60, length($2)), $3, substr($4, length($4)-60, length($4))}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_last60.fa
		
		#transform the data so that it looks like a fastq file
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_last60.fa | 
		awk -v OFS="\t" -v FS="\t" '{print $1" 1:N:0:NNNNNNNN+NNNNNNNN", $2, $3, $4}' |
		sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_last60.fastq
	
		#The mapping
		echo "Premapping 2.b" ${i}
		echo "###########################################################################"
		bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_last60.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_last60.sam	

		#acquiring the names of those that end up getting a larger MQ
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_last60.sam | 
		sed -n -E '/^[[:alpha:]]/p' |
		awk -v OFS="\t" -v FS="\t" '{if ($5 > 0) {print $1}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0_60.txt
		

		#then combining the good ones from this set with the good ones from the 30bp set
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0_60.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0_comb.txt
	else
		#also make a comb file even if all were alright in the 30bp set
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0_comb.txt
	fi
	
	#then get back the rest of the data for the combined list of good starting reads
	readarray -t lines < ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_start_w_TDNA_names_MQLT0_comb.txt
	for r in "${lines[@]}"
	do
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | egrep ${r}  >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_start_w_TDNA_names_MQLT0_comb.fa
	done	
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_start_w_TDNA_names_MQLT0_comb.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, length($2)-60, length($2)), $3, substr($4, length($4)-60, length($4))}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_fw_comb.fa
	
	#then do the same for the good ending reads
	readarray -t lines < ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_end_w_TDNA_names_MQLT0_comb.txt
	for r in "${lines[@]}"
	do
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | egrep ${r}  >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_end_w_TDNA_names_MQLT0_comb.fa
	done	
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_end_w_TDNA_names_MQLT0_comb.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1, substr($2, length($2)-60, length($2)), $3, substr($4, length($4)-60, length($4))}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb.fa
	
	#reverse complement the sequences, and combine them with names again
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb.fa |
	awk -v OFS="\t" -v FS="\t" '{print $2}' |
	tr ACGTacgt TGCAtgca | 
	rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb_seq_RC.fa
	paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb.fa ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb_seq_RC.fa |
	awk -v OFS="\t" -v FS="\t" '{print $1, $5, $3, $4}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb_Rced.fa
	
	#put the good starting and good ending reads back together into one file
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_fw_comb.fa ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_rev_comb_Rced.fa |
	sort |
	uniq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_total_combined.fa
	exit 1

	
################################################################################################################
#FASTA mode: Proper mapping
################################################################################################################

	#transform the data so that it looks like a fastq file
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $1" 1:N:0:NNNNNNNN+NNNNNNNN", $2, $3, $4}' |
	sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_R1.fastq
			
	#create the reverse complement data
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa | 
	awk -v OFS="\t" -v FS="\t" '{print $2}' |
	tr ACGTacgt TGCAtgca | 
	rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns_RC.fa
			
	#then perform a transformation to produce fastq format data of the R2 (which is just R1 but RCed)
	paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_columns.fa ${WORKPATH}/${i}/${i}_columns_RC.fa |
	awk -v OFS="\t" -v FS="\t" '{print $1" 2:N:0:NNNNNNNN+NNNNNNNN", $5, $3, $4}' |
	sed 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_R2.fastq

	echo "Mapping" ${i}
	echo "###########################################################################"
	bwa-mem2 mem ${WORKPATH}/${CURRENTREF} ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_R1.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_R2.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam
		
	samtools view -1 ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.bam
	samtools sort -o ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted.bam ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.bam
	
################################################################################################################
#FASTA mode: dedupping
################################################################################################################
	
	case $DEDUPOPT in
		OPT) #optical duplicate filtering
			echo "Optical duplicate filtering of ${i}"
			echo "###########################################################################"
			picard MarkDuplicates --INPUT ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted.bam --OUTPUT ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam --METRICS_FILE ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted_dedup_metrics.txt --REMOVE_DUPLICATES TRUE
			echo "###########################################################################"
			samtools index ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam.bai
			;;
		UMI) #UMI consolidation
			echo "No UMI consolidation during fasta mode, skipping filtering of ${i}"
			cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam
			samtools index ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam.bai
			;;
		OFF) #no duplicate filtering
			echo "Skipping duplicate filtering of ${i}"
			cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam
			samtools index ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam.bai
			;;
		*) #default, no duplicate filtering
			echo "Skipping duplicate filtering of ${i}"
			cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_sorted.bam > ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam
			samtools index ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam.bai
			;;
	esac
	
################################################################################################################	
#FASTA mode: Filtering reads
################################################################################################################

	echo "Creating empty output files"
	> ${WORKPATH}/file1.temp | awk -v OFS="\t" -v FS="\t" ' BEGIN{print "QNAME", "RNAME_1", "POS_1", "CIGAR_1", "SEQ_1", "QUAL_1", "SATAG_1", "SEQ_RCed_1", "RNAME_2", "POS_2", "CIGAR_2", "SEQ_2", "QUAL_2", "SATAG_2", "SEQ_RCed_2", "UMI", "FILE_NAME", "PRIMER_SEQ"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_A.txt


	echo "Processing fw reads of ${i}"
	samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800 | samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt
		
	echo "Processing rev reads of ${i}"
	samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt
		
	echo "Rev complementing rev reads of ${i}"
	awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt
	paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt

	echo "Combining fw and rev reads"
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt
	awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt

	echo "Processing fw mates of accepted reads of ${i}"
	samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt

	echo "Processing rev mates of accepted reads of ${i}"
	samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt

	echo "Rev complementing rev mates of ${i}"
	awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt
	paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt

	echo "Combining fw and rev mates of ${i}"
	cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt | sort -T ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/temp/ -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt

	echo "Combining selected reads and mates of ${i}"
	join -j 1 -o 1.1,1.3,1.4,1.6,1.10,1.11,1.12,1.13,2.3,2.4,2.6,2.10,2.11,2.12,2.13 -t $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" ' {print $0, "NA", i, "NA"}'  >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}_A.txt


################################################################################################################
#FASTA mode: counting reads
################################################################################################################
	
	echo "Counting reads of ${CURRENTSAMPLE} ${CURRENTRUNID}"
	RAWNO=$(( $(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq | wc -l)/4 )) 
	echo "RAWNO: ${RAWNO}"
	DEDUPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/dedup_names.txt  | wc -l)
	echo "DEDUPNO: ${DEDUPNO}"
	MAPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam | grep -Ev '^(\@)' | awk '$3 != "*" {print $0}' | sort -u -t$'\t' -k1,1 | wc -l)
	echo "MAPNO: ${MAPNO}"
	PREPRONO=$(( $(cat ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt | wc -l)-1 ))
	echo "PREPRONO: ${PREPRONO}"
	cat ${WORKPATH}/stats/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v RAWNO="${RAWNO}" -v MAPNO="${MAPNO}" -v DEDUPNO="${DEDUPNO}" -v PREPRONO="${PREPRONO}" ' END{print CURRENTSAMPLE, CURRENTRUNID, RAWNO, MAPNO, DEDUPNO, PREPRONO}' >> ${WORKPATH}/stats/read_numbers.txt
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTA mode: removing temporary files
################################################################################################################
	
	echo "Removing temporary files"
	rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
	rm ${WORKPATH}/file1.temp 
	echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

done
else
################################################################################################################
#FASTQ_mode
################################################################################################################
echo "Analysing files in regular fastq mode"
for i in "${LIST_FILES[@]}" 
do
################################################################################################################
#FASTQ_mode: Reading variables
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
#FASTQ_mode: Creating output folder
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
#FASTQ_mode: looking for data and unzipping
################################################################################################################
echo "checking for the presence of the sequencing data files..."
if [[ -f "${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz, unzipping ..."
	gunzip -c ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq
else
	echo "Did not find ${WORKPATH}/${i}${CURRENTR1SUFFIX}.fastq.gz, moving to the next sample";
	rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
	continue
fi
if [[ -f "${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz" ]]
then
	echo "Found ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz, unzipping ..."
	gunzip -c ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR2SUFFIX}.fastq
else
	echo "Did not find ${WORKPATH}/${i}${CURRENTR2SUFFIX}.fastq.gz, moving to the next sample";
	continue
fi

case $DEDUPOPT in
	UMI)#checking for the presence of the UMI file 
		CURRENTUMISUFFIX=$(cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $22}}')
		echo "Current UMI suffix: ${CURRENTUMISUFFIX}"
		if [[ -f "${WORKPATH}/${i}${CURRENTUMISUFFIX}.fastq.gz" ]]
		then
			echo "Found ${WORKPATH}/${i}${CURRENTUMISUFFIX}.fastq.gz, unzipping..."
			gunzip -c ${WORKPATH}/${i}${CURRENTUMISUFFIX}.fastq.gz > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTUMISUFFIX}.fastq
		else
			echo "Did not find ${WORKPATH}/${i}${CURRENTUMISUFFIX}.fastq.gz, moving to the next sample";
			continue
		fi
		;;
	*) 	#no need to check for the presence of the UMI file 
		echo "Skipping UMI"
		;;
esac
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: fastq to tab separated
################################################################################################################

echo "Reformatting R1 data..."
awk -v FS="\t" 'ORS=NR%4?FS:RS' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq |
sed -E 's/[[:space:]]/\t/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.txt

echo "Reformatting R2 data..."
awk -v FS="\t" 'ORS=NR%4?FS:RS' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR2SUFFIX}.fastq |
sed -E 's/[[:space:]]/\t/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2.txt


case $DEDUPOPT in
	UMI) #UMI consolidation
		echo "Creating list of read names and UMIs..."
		awk -v FS="\t" 'ORS=NR%4?FS:RS' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTUMISUFFIX}.fastq | awk -v FS="\t|[ ]" -v OFS="\t" '{print $1, $3}'> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/names_identifiers.txt
		;;
	OPT)
		#making a list of fake "UMIs", which are all the same. This way proper optical duplicate filtering can take place
		echo "Creating list of names without UMIs..."
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.txt | awk -v OFS="\t" -v FS="\t" '{print $1, "AAAAAAAA"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/names_identifiers.txt
		;;
	*) 	#skip filtering
		;;
esac
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: dedupping
################################################################################################################

case $DEDUPOPT in
	UMI|OPT) #from here UMI and OPT filtering use the same code
		echo "Performing duplicate filtering..."
		#first creating a file that contains the R1, R2 and UMI
		paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/names_identifiers.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined.txt
		#then sort and only leaving the unique ones, based only on the sequences of course
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined.txt | sort -u -k3,3 -k8,8 -k12 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined_dedupped.txt
		awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined_dedupped.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/dedup_names.txt
		#creating a list of duplicate reads, for troubleshooting
		grep -v -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/dedup_names.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined.txt > ${WORKPATH}/stats/${CURRENTSAMPLE}_${CURRENTRUNID}_Duplicate_pairs.txt
		#creating deduplicated R1 file in fastq format
		awk -v OFS="\t" -v FS="\t" '{print $1, $3, $4, $5}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined_dedupped.txt | sed -E 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_dedupped.fastq
		#creating deduplicated R2 file in fastq format
		awk -v OFS="\t" -v FS="\t" '{print $6, $8, $9, $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_R2_UMI_combined_dedupped.txt | sed -E 's/\t/\n/g' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2_dedupped.fastq
		;;
	*) 	#skip filtering
		echo "Skipping duplicate filtering..."
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_dedupped.fastq
		cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR2SUFFIX}.fastq > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2_dedupped.fastq
		;;
esac
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: Trimming
################################################################################################################

echo ">p5" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
echo ">p5_RC" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $5}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
echo ">p7" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
echo ">p7_RC" >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa
cat ${WORKPATH}/Sample_information.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($3==i) {print $6}}' | tr "[:lower:]" "[:upper:]" | tr "RYMKSWBDHV" "NNNNNNNNNN" | tr "ACGT" "TGCA" | rev >> ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa

echo "Trimming" ${i}
echo "###########################################################################"
trimmomatic PE ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R1_dedupped.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/R2_dedupped.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_forward_unpaired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_paired.fastq ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}_reverse_unpaired.fastq ILLUMINACLIP:${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Illumina_adapters.fa:2:30:10:1:TRUE CROP:${CURRENTTRIMLEN} -phred33
echo "###########################################################################"
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: Mapping
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
#FASTQ_mode: Filtering reads
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


echo "Processing fw reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800 | samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Processing rev reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x80 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Rev complementing rev reads of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Combining fw and rev reads and keeping only those starting with primer seq ${PRIMERSEQ}"
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_rv_RCed.txt | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" -v test="$PRIMERSEQ" '$10 ~ "^"test {print}'  > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt
awk -v OFS="\t" -v FS="\t" '{print $1}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Processing fw mates of accepted reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Processing rev mates of accepted reads of ${i}"
samtools view -uF 0x100 ${WORKPATH}/bams/${CURRENTSAMPLE}_${CURRENTRUNID}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' | grep -f ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all_names.txt > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Rev complementing rev mates of ${i}"
awk -v OFS="\t" -v FS="\t" '{print $10}' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt | tr ACGTacgt TGCAtgca | rev > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt
paste ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Combining fw and rev mates of ${i}"
cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_fw.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_rv_RCed.txt | sort -k 1,1 > ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

echo "Combining selected reads and mates of ${i}"
join -j 1 -o 1.1,1.3,1.4,1.6,1.10,1.11,1.12,1.13,2.3,2.4,2.6,2.10,2.11,2.12,2.13 -t $'\t' ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Primer_reads_all.txt ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/Mates_all.txt | 
awk -v OFS="\t" -v FS="\t" -v i="$i" -v PRIMERSEQ="$PRIMERSEQ" -v DEDUPOPT="$DEDUPOPT" -v CURRENTTRIMLEN="$CURRENTTRIMLEN" ' {print $0, i, PRIMERSEQ, DEDUPOPT, CURRENTTRIMLEN}'  >> ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: Statistics
################################################################################################################

echo "Counting reads of ${CURRENTSAMPLE} ${CURRENTRUNID}"
RAWNO=$(( $(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${i}${CURRENTR1SUFFIX}.fastq | wc -l)/4 )) 
echo "RAWNO: ${RAWNO}"
cat ${WORKPATH}/stats/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v RAWNO="${RAWNO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Raw", RAWNO}' >> ${WORKPATH}/stats/read_numbers.txt

case $DEDUPOPT in
	UMI|OPT) 
		DEDUPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/dedup_names.txt  | wc -l)
		;;
	*) 	
		DEDUPNO=RAWNO
		;;
esac
echo "DEDUPNO: ${DEDUPNO}"
cat ${WORKPATH}/stats/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v DEDUPNO="${DEDUPNO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Dedupped", DEDUPNO}' >> ${WORKPATH}/stats/read_numbers.txt

MAPNO=$(cat ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}/${CURRENTSAMPLE}.sam | grep -Ev '^(\@)' | awk '$3 != "*" {print $0}' | sort -u -t$'\t' -k1,1 | wc -l)
echo "MAPNO: ${MAPNO}"
cat ${WORKPATH}/stats/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v MAPNO="${MAPNO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Mapped", MAPNO}' >> ${WORKPATH}/stats/read_numbers.txt

PREPRONO=$(( $(cat ${WORKPATH}/input/${CURRENTSAMPLE}_${CURRENTRUNID}_A.txt | wc -l)-1 ))
echo "PREPRONO: ${PREPRONO}"
cat ${WORKPATH}/stats/read_numbers.txt | awk -v OFS="\t" -v FS="\t" -v CURRENTFILE="${i}" -v CURRENTSAMPLE="${CURRENTSAMPLE}" -v CURRENTRUNID="${CURRENTRUNID}" -v PREPRONO="${PREPRONO}" -v FOCUSLOCUS="${FOCUSLOCUS}" ' END{print CURRENTSAMPLE, CURRENTRUNID, CURRENTFILE, FOCUSLOCUS, "Filtered", PREPRONO}' >> ${WORKPATH}/stats/read_numbers.txt


echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

################################################################################################################
#FASTQ_mode: Cleanup
################################################################################################################

echo "Removing temporary files"

rm -r ${WORKPATH}/${CURRENTSAMPLE}_${CURRENTRUNID}
rm ${WORKPATH}/file1.temp 
echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

done

fi

echo "Removing temporary files"
rm ${WORKPATH}/file0.temp 

echo "## $(( $(date +%s) - ${StartTime} )) seconds elapsed ##"

exit 0

