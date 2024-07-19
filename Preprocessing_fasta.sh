#set these variables
WORKPATH=workdir
echo Workpath =  ${WORKPATH}
GENBANKFILE=Sequences_detailed_all.txt #contains all sequences in genbank format. It should also contain the plasmid name somewhere.
echo Genbank file =  ${GENBANKFILE}
PLASMIDNAMES=(pADIS1 pROK2 pGABI1 pAC106 pAC161) #all possible plasmid names
echo Plasmids = ${PLASMIDNAMES[@]}

#first replace all line separators with spaces
#then turn double forward slash followed by space into line separator
#then turn a variable number of spaces into a single tab
#then delete single spaces
echo parsing genbank file
cat ${WORKPATH}/${GENBANKFILE} | tr '\n' ' ' | sed 's/\/\/ /\n/g' | sed 's/ \+ /\t/g' | tr -d ' ' > ${WORKPATH}/sequences_parsed_messy.tmp
cat ${WORKPATH}/sequences_parsed_messy.tmp | sed -n -e 's/^.*ORIGIN//p' | sed 's/\t//g' | sed 's/[0-9]//g' > ${WORKPATH}/sequences_seqs_only.tmp
cat ${WORKPATH}/sequences_parsed_messy.tmp | awk '{print $2}' > ${WORKPATH}/sequences_accnos_only.tmp
paste ${WORKPATH}/sequences_accnos_only.tmp ${WORKPATH}/sequences_seqs_only.tmp ${WORKPATH}/sequences_parsed_messy.tmp > ${WORKPATH}/sequences_parsed.tmp

#then make separate lists of accessions based on what plasmid is mentioned in the comments
#and prepare the files so that they resemble fastq files from illumina sequencing
#and zip them
for i in "${PLASMIDNAMES[@]}"
do
echo Processing ${i}
cat ${WORKPATH}/sequences_parsed.tmp | grep -F ${i} | awk -v FS='\t' -v OFS='\t' '{print "@"$1" fasta", $2, "+", $2}' | awk -v FS='\t' -v OFS='\t' '{gsub("[ctgan]", "G", $4); print}' | tr "[:lower:]" "[:upper:]" > ${WORKPATH}/sequences_clean_${i}.tmp
LINECOUNT=$(wc -l < ${WORKPATH}/sequences_clean_${i}.tmp)
if [[ ${LINECOUNT} = 0 ]] 
then
	continue
fi
cat ${WORKPATH}/sequences_clean_${i}.tmp | tr "\t" "\n" > ${WORKPATH}/${i}_R2.fastq
cat ${WORKPATH}/sequences_clean_${i}.tmp | awk -v FS='\t' -v OFS='\t' '{print $2}' | tr "ACGTN" "TGCAN" | rev >  ${WORKPATH}/sequences_${i}_RCed.tmp
paste ${WORKPATH}/sequences_${i}_RCed.tmp ${WORKPATH}/sequences_clean_${i}.tmp | awk -v FS='\t' -v OFS='\t' '{print $2, $1, $4, $5}'| tr "\t" "\n" > ${WORKPATH}/${i}_R1.fastq
echo Finished ${i}
done
echo gzipping
gzip ${WORKPATH}/*.fastq
echo cleanup
rm ${WORKPATH}/*.tmp


