#set these variables
WORKPATH=workdir
echo Workpath =  ${WORKPATH}
FASTAFILES=$(ls ${WORKPATH}/*.fasta | sed 's/.fasta//g')
echo Processing ${FASTAFILES}

for i in ${FASTAFILES}
do

#first make the list-like fasta into a file with 2 columns: sequence name and sequence
#then make the file resemble a fastq file from illumina sequencing 

echo parsing fasta file
cat ${i}.fasta | grep -o -P '(?<=\>).*?(?=\s)' > ${WORKPATH}/sequences_names_only.tmp
cat ${i}.fasta | sed 's/>.*$/newline/g' | tr '\n' ' '| sed 's/\s//g' | sed 's/newline/\n/g'| awk 'FNR>1{print}' > ${WORKPATH}/sequences_seqs_only.tmp
paste ${WORKPATH}/sequences_names_only.tmp ${WORKPATH}/sequences_seqs_only.tmp | awk -v FS='\t' -v OFS='\t' '{print "@"$1" fasta", $2, "+", $2}' | awk -v FS='\t' -v OFS='\t' '{gsub("[CTGANctgan]", "G", $4); print}' | tr "[:lower:]" "[:upper:]" > ${WORKPATH}/sequences_clean.tmp
cat ${WORKPATH}/sequences_clean.tmp | tr "\t" "\n" > ${i}_R2.fastq
cat ${WORKPATH}/sequences_clean.tmp | awk -v FS='\t' -v OFS='\t' '{print $2}' | tr "ACGTN" "TGCAN" | rev >  ${WORKPATH}/sequences_RCed.tmp
paste ${WORKPATH}/sequences_RCed.tmp ${WORKPATH}/sequences_clean.tmp | awk -v FS='\t' -v OFS='\t' '{print $2, $1, $4, $5}'| tr "\t" "\n" > ${i}_R1.fastq
echo Finished ${i}

done

#and gzip them
echo gzipping
gzip ${WORKPATH}/*.fastq
echo cleanup
rm ${WORKPATH}/*.tmp


