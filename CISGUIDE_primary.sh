#!/bin/sh

#this code does the following:
#remove secondary and supplementary alignment reads
#removes reads with unmapped mates
#removes reads shorter than 90bp
#the software will process all the bamfiles in the work directory
#sorting happens in a custom temp directory, otherwise it can run out of memory
#it requires likely at least 8gb of RAM
#note that reverse complementing only affects the sequence, not the quals
#note that duplicates are not removed here. They are removed with picardtools (and possibly combining reads with same UMI.

#used commands in samtools:
#-u : uncompressed, use this when piping to another samtools command
#-f : only output reads with the indicated flag
#-F : do not output reads with indicated flag

#flags used:

#0x4 	: 	UNMAP 			segment unmapped
#0x8 	: 	MUNMAP 			next segment in the template unmapped
#0x10 	: 	REVERSE 		seq is reverse complemented
#0x40	:	READ1			the first segment in the template
#0x80 	: 	READ2 			the last segment in the template
#0x100 	:	SECONDARY		secondary alignment
#0x400	:	DUP				PCR or optical duplicate
#0x800 	: 	SUPPLEMENTARY 	supplementary alignment





LIST_SAMPLES=($(ls shared/*.bam | tr '\n' ' ' | sed 's/.sorted.bam//g' | sed 's/shared\///g'))

echo "Processing the following samples:" ${LIST_SAMPLES[*]}
now=$(date)
echo "Starting at $now"
echo "########################################################################"

for i in "${LIST_SAMPLES[@]}" 

do

if [ ! -d "shared/${i}/" ] 
then
	echo "creating dir  ${i}"
	mkdir shared/${i}/
	mkdir shared/${i}/temp/
fi

echo "Looking for the primer sequence"

if [ -f "shared/Sample_information.txt" ]
then
	echo "Found Sample_information.txt"
else
	echo "Did not find Sample_information.txt, exiting";
	exit 1
fi

awk -v OFS="\t" -v FS="\t" -v i="$i" 'FNR>1{if($1==i) {print $2}}' shared/Sample_information.txt > shared/${i}/PS.txt

PRIMERSEQFULL=$(cat shared/${i}/PS.txt)
PRIMERSEQ="$( echo "$PRIMERSEQFULL" | sed -e 's#^TCAGACGTGTGCTCTTCCGATCT##' )"

if [ -z "$PRIMERSEQ" ]
then
	echo "Primer sequence not found, moving to the next sample"
	rm -r shared/${i}
	echo "########################################################################"
	continue
else
	echo "Using this primer sequence: ${PRIMERSEQ}"
fi

echo "Creating empty output files"
> shared/file1.temp | awk -v OFS="\t" -v FS="\t" ' BEGIN{print "QNAME", "RNAME_1", "POS_1", "CIGAR_1", "SEQ_1", "QUAL_1", "SATAG_1", "SEQ_RCed_1", "RNAME_2", "POS_2", "CIGAR_2", "SEQ_2", "QUAL_2", "SATAG_2", "SEQ_RCed_2", "FILE_NAME", "PRIMER_SEQ"}' > shared/${i}_A.txt

echo "Processing fw reads of ${i}"

samtools view -uF 0x100 shared/${i}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800 | samtools view -uf 0x80 |samtools view -uF 0x8 | samtools view -F 0x10 | sort -T shared/${i}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' > shared/${i}/Primer_reads_fw.txt

echo "Processing rev reads of ${i}"

samtools view -uF 0x100 shared/${i}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x80 |samtools view -uF 0x8 | samtools view -f 0x10 | sort -T shared/${i}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50 && length($10)>89){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' > shared/${i}/Primer_reads_rv.txt

echo "Rev complementing rev reads of ${i}"

awk -v OFS="\t" -v FS="\t" '{print $10}' shared/${i}/Primer_reads_rv.txt | tr ACGTacgt TGCAtgca | rev > shared/${i}/Primer_reads_rv_RCseqs.txt

paste shared/${i}/Primer_reads_rv.txt shared/${i}/Primer_reads_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > shared/${i}/Primer_reads_rv_RCed.txt

echo "Combining fw and rev reads and keeping only those starting with primer seq ${PRIMERSEQ}"

cat shared/${i}/Primer_reads_fw.txt shared/${i}/Primer_reads_rv_RCed.txt | sort -T shared/${i}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" -v test="$PRIMERSEQ" '$10 ~ "^"test {print}'  > shared/${i}/Primer_reads_all.txt

awk -v OFS="\t" -v FS="\t" '{print $1}' shared/${i}/Primer_reads_all.txt > shared/${i}/Primer_reads_all_names.txt

echo "Processing fw mates of accepted reads of ${i}"

samtools view -uF 0x100 shared/${i}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -f 0x10 | sort -T shared/${i}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "FALSE"}}' | grep -f shared/${i}/Primer_reads_all_names.txt > shared/${i}/Mates_fw.txt

echo "Processing rev mates of accepted reads of ${i}"

samtools view -uF 0x100 shared/${i}.sorted.bam | samtools view -uF 0x4 | samtools view -uF 0x800| samtools view -uf 0x40 | samtools view -uF 0x8 | samtools view -F 0x10 | sort -T shared/${i}/temp/ -k 1,1 | awk -v OFS="\t" -v FS="\t" '{ if ($5>50){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' | grep -f shared/${i}/Primer_reads_all_names.txt > shared/${i}/Mates_rv.txt

echo "Rev complementing rev mates of ${i}"

awk -v OFS="\t" -v FS="\t" '{print $10}' shared/${i}/Mates_rv.txt | tr ACGTacgt TGCAtgca | rev > shared/${i}/Mates_rv_RCseqs.txt

paste shared/${i}/Mates_rv.txt shared/${i}/Mates_rv_RCseqs.txt | awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $11, $12, "TRUE"}' > shared/${i}/Mates_rv_RCed.txt

echo "Combining fw and rev mates of ${i}"

cat shared/${i}/Mates_fw.txt shared/${i}/Mates_rv_RCed.txt | sort -T shared/${i}/temp/ -k 1,1 > shared/${i}/Mates_all.txt

echo "Combining selected reads and mates of ${i}"

join -j 1 -o 1.1,1.3,1.4,1.6,1.10,1.11,1.12,1.13,2.3,2.4,2.6,2.10,2.11,2.12,2.13 -t $'\t' shared/${i}/Primer_reads_all.txt shared/${i}/Mates_all.txt | awk -v OFS="\t" -v FS="\t" -v i="$i" -v PRIMERSEQ="$PRIMERSEQ" ' {print $0, i, PRIMERSEQ}'  >> shared/${i}_A.txt

echo "Removing dir ${i}"

rm -r shared/${i}

now=$(date)
echo "Finished at $now"
echo "########################################################################"

done

rm shared/file1.temp

exit 0

