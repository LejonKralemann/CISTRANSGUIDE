First you need to create a Sample_information.txt file. It needs to have the following fields, in the following order. Include a header row.

Sample:	The name of the library sample
Primer:	The secondary GUIDEseq primer sequence (GSP2)
Fastq_name:	the part before "_R1.fastq.gz" or "_R2.fastq.gz"
Ref:	the genome reference file name, containing all the chromosomes, as well as the plasmid sequence(s). Note that it is desirable to have 
P5_adapter:	Enter the full P5 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
P7_adapter:	Enter the full p7 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
Plasmid:	the name of the plasmid as written in the reference file.
DSB_chrom:	the chromosome or plasmid on which the primer sequence is located, should agree with the name in the reference fasta.
FlankAUltEnd:	the position on the DSB_chrom where flank A ends (if there are multiple possibilities, take the furthest downstream). Flank A is the flank starting with the primer. Flank A and Flank B may overlap.
Locus_name: a name denoting the locus in which a DSB is created, or T-DNA border that is focused on
FLANK_A_ORIENT: The orientation of the primer and thus of flank A.
FLANKBUltStart:	the position on the DSB chrom on where flank B starts (if there are multiple possibilities, take the furthest upstream). Flank B is the flank after the nick or DSB. Flank B and flank A may overlap.
Genotype: The name of the genotype
Plasmid_alt: the name of the second plasmid, if a transformation with multiple plasmids has been performed. Enter "NA" if you only used 1 plasmid.
DNA: name of the DNA sample from which libraries were made. Often one DNA sample is used to do an LB, RB, FW, and RV reaction.

Place the CISGUIDE_primary.sh script, Sample_information.txt, reference fastas ("PLASMID_NAME.fa"), raw sequencing files ("NAME_R1.fastq.gz" and "NAME_R2.fastq.gz") in your work directory.
Then run the CISGUIDE_primary.sh script, by doing the following:
bash shared/CISGUIDE_primary.sh |& tee -a shared/primary_log.txt

The output is a folder per library, containing a .bam and .bam.bai file so that sequences may be visualized in a genome browser, and a text file ending in "_A.txt" containing all preprocessed reads.

Next, place the files that end in "_A.txt", the reference fastas, and Sample_information.txt in the CISGUIDE input folder. Then run the R script.

The output is one excel file per library showing on every row a unique event. Another excel file contains all data together.