----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Install and run
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- Install bioconda. Download the .sh file from the website www.anaconda.org and run it (bash miniconda.sh). 
- Then install samtools, trimmomatic, and bwa-mem2 via the command "conda install -c bioconda <samtools/picard/trimmomatic/bwa-mem2>"
- Then install dos2unix with "apt install dos2unix"
- Place CISGUIDE_primary.sh, Sample_information.txt, fastq.gz files, and reference fastas ("PLASMID_NAME.fa") in your work directory.
- Run the software with something like: bash shared/CISGUIDE_primary.sh -p shared -f FALSE -d OPT -t 150 |& tee -a shared/primary_log.txt
- The output is a folder per library, containing a .bam and .bam.bai file so that sequences may be visualized in a genome browser, and a text file ending in "_A.txt" containing all preprocessed reads.
- Next, place the files that end in "_A.txt", the reference fastas (change minus sign to underscore), and Sample_information.txt in the CISGUIDE input folder. Then run the R script.
- The output is one excel file per library showing on every row a unique event. Another excel file contains all data together.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Sample_information
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
You need to create a file named "Sample_information.txt". This file needs to have one row per sample, each with the fields indicated below. Include a header row.

Sample:			The name of the library sample. Example: "7C1_LB". This name should in principle be unique. If the samples have been sequenced multiple times, they can have the same name, as long as the RunID is different.
Primer:			The secondary GUIDEseq primer sequence (GSP2).
File_name:		The part before "_R1.fastq.gz" or "_R2.fastq.gz".
Ref:			The genome reference file name, containing all the chromosomes, as well as the plasmid sequence(s). Example: "pCAS-PPO.fa".
P5_adapter:		Enter the full P5 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
P7_adapter:		Enter the full p7 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
Plasmid:		The name of the plasmid as written in the reference file. Example: "pCAS-PPO".
DSB_chrom:		The chromosome or plasmid on which the primer sequence is located, should agree with the name in the reference fasta. Example:"Chr4".
FlankAUltEnd:	The position on the DSB_chrom where flank A ends (if there are multiple possibilities, take the furthest downstream). Flank A is the flank starting with the primer. Flank A and Flank B may overlap.
Locus_name: 	A name denoting the locus in which a DSB is created, or T-DNA border that is focused on.
FLANK_A_ORIENT: The orientation of the primer and thus of flank A, relative to the reference genome sequence. FW or RV.
FlankBUltStart:	The position on the DSB chrom on where flank B starts (if there are multiple possibilities, take the furthest upstream). Flank B is the flank after the nick or DSB. Flank B and flank A may overlap.
Genotype: 		The name of the genotype of the organism that was sampled.
Plasmid_alt: 	The name of the second plasmid, if a transformation with multiple plasmids has been performed. Enter "NA" if you only used 1 plasmid.
DNA: 			Name of the DNA sample from which libraries were made. Often one DNA sample is used to do an LB, RB, FW, and RV reaction. Example: "7C1".
RunID: 			Name of the sequencing run. This allows one to differentiate between samples that were sequenced multiple times, or to limit comparisons within runs to avoid effects of different sequencing depths between runs.
LBSeq: 			NA, only used for fasta mode
RBSeq: 			NA, only used for fasta mode
Ecotype: 		The name of the ecotype of the organism that was samples. E.g. "Col-0".
R1Suffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for read 1.
R2Suffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for read 2.
UMISuffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for the UMI.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
General notes
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- The reference fasta needs to have chromosome names starting with "Chr" (e.g. "Chr4"), not just a number.
- Make sure that the reference fastas don't contain names with minus signs.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Special notes for fasta mode
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Sample:	The name of the list of sequences in fasta format.
Primer:	Use "NA"
File_name:	the part before ".fasta" 
Ref:	Same as in fastq mode
P5_adapter:	Use "NA"
P7_adapter:	Use "NA"
Plasmid:	Same as in fastq mode
DSB_chrom:	Same as in fastq mode
FlankAUltEnd:	Use "NA"
Locus_name: Use "NA"
FLANK_A_ORIENT: Use "NA"
FlankBUltStart:	Use "NA"
Genotype: Same as in fastq mode
Plasmid_alt: Same as in fastq mode
DNA: 
RunID: Same as in fastq mode
LBseq: Sequence of the left border repeat in the binary plasmid
Rbseq: Sequence of the right border repeat in the binary plasmid