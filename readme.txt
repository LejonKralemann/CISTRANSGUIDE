----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Install and run
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- Install conda and packages: 
	Download the .sh file from the website www.anaconda.org and run it (bash miniconda.sh). 
	Add channels for bioconda:
		conda create -n bio
		conda activate bio
		conda config --add channels defaults
		conda config --add channels bioconda
		conda config --add channels conda-forge
		conda config --set channel_priority strict
	Note that conda may start per default in base, in which case you need to do "conda deactivate" and then "conda activate bio"
	Then install samtools, picard, trimmomatic, and bowtie2 via the command "conda install -c bioconda <package>". Replace <package> with the package name.
- Then install dos2unix with "apt install dos2unix"
- Place CISGUIDE_primary2.sh, Sample_information.txt, fastq.gz files, and reference fastas ("PLASMID_NAME.fa") in your work directory.
- Run the software with something like: (replace "workdir" with the name of the folder where the input files are) 
	conda activate bio
	bash workdir/CISGUIDE_primary2.sh -p workdir -f FALSE -t 150 |& tee -a workdir/primary_log.txt
- The output is a folder per library, containing a .bam and .bam.bai file so that sequences may be visualized in a genome browser, and a text file ending in "_A.txt" containing all preprocessed reads.
- Next, place the files that end in "_A.txt", the read_numbers.txt file, the reference fastas (change minus sign to underscore), and Sample_information.txt in the CISGUIDE input folder. Then run the R script (CISGUIDE2.R).
- The output is one excel file per library showing on every row a unique event. Another excel file contains all data together.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Sample_information
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
You need to create a file named "Sample_information.txt". This file needs to have one row per sample, each with the fields indicated below. Include a header row. The columns need to occur in the following order.

Sample:						The name of the library sample. Example: "7C1_LB". This name should in principle be unique. If the samples have been sequenced multiple times, they can have the same name, as long as the RunID is different.
Primer:						The secondary GUIDEseq primer sequence (GSP2). Use 'NA' in fasta mode. For competitive TRANSGUIDE, extend the primer seq with the barcode directly downstream to only acquire the events from 1 sample. Repeat with the second barcode.
File_name:					The part before "_R1.fastq.gz" or "_R2.fastq.gz". This name needs to be unique in this list.
Ref:						The genome reference file name, containing all the chromosomes, as well as the plasmid sequence(s). Example: "pCAS-PPO.fa". Note that the Plasmid sequence in the ref should "break" in the backbone part, not in the T-DNA part.
P5_adapter:					Enter the full P5 adapter primer, without special characters, in upper case. Degenerate bases allowed but will be automatically changed to N. Use 'NA' in fasta mode.
P7_adapter:					Enter the full p7 adapter primer, without special characters, in upper case. Degenerate bases allowed but will be automatically changed to N. Use 'NA' in fasta mode.
Plasmid:					The name of the plasmid as written in the reference file. Example: "pCAS-PPO".
DSB_CONTIG:					The name of the chromosome on which a DSB is induced (e.g. "Chr4"). If no DSB induction use "NA".
DSB_FW_END:					The final position before (on 5' end of) the induced DSB, including overhangs.
Locus_name:					The name of the locus in which a DSB is made. Use NA when no DSB is made.
Focus_contig_name: 			The contig of the primers. Can be the plasmid name (TRANSGUIDE), or the name of a chromosome (CISGUIDE).
DSB_OVERHANG:				If the induced DSB has an overhang, indicate its length here.
Genotype: 					The name of the genotype of the organism that was sampled.
Plasmid_alt: 				The name of the second plasmid, if a transformation with multiple plasmids has been performed. Enter "NA" if you only used 1 plasmid.
DNA: 						Name of the DNA sample from which libraries were made. Often one DNA sample is used to do an LB, RB, FW, and RV reaction. Example: "7C1".
RunID: 						Name of the sequencing run. This allows one to differentiate between samples that were sequenced multiple times, or to limit comparisons within runs to avoid effects of different sequencing depths between runs.
TDNA_LB_END:				The most extreme position on the T-DNA at the LB end. If 'NA', then the program will look for the position.
TDNA_RB_END:				The most extreme position on the T-DNA at the RB end. If 'NA', then the program will look for the position.
Ecotype: 					The name of the ecotype of the organism that was sampled. E.g. "Col-0".
R1Suffix: 					The part of the filename after File_name, but before the extention .fastq.gz, for read 1 (the read sequenced from the p5 adapter).
R2Suffix: 					The part of the filename after File_name, but before the extention .fastq.gz, for read 2 (the read sequenced from the p7 adapter).
UMISuffix: 					The part of the filename after File_name, but before the extention .fastq.gz, for the UMI.
Family:						Indicate here if some samples are related and therefore junctions can occur in multiple samples. Use integers starting from "1". Use "0" in all rows not to include family information. This information is used to keep some junctions which otherwise would be filtered out by the duplicate position filter. 
AgroGeno:					Genotype/strain of Agrobacterium that was used in the transformation.
TDNA_ALT_LB_END:			LB position on the alternative T-DNA plasmid. "NA" if no alternative plasmid, or if you want the program to look for the position.
TDNA_ALT_RB_END:			Same as previous but for RB.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
General notes
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- The reference fasta needs to have chromosome names starting with "Chr" (e.g. "Chr4"), not just a number.
- Make sure that the reference fastas don't contain names with minus signs. Not inside the file, nor in the name of the file.
- You can run the preprocessing software with options:

-p	Set work path. default: home directory
-f	Switches to fasta mode if TRUE. default: FALSE.
-t 	Trimming length. Value indicate maximum number of nt to keep. default: 999999 (meaning no trimming is performed).

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
About the output
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In the case of microhomology, the output gives you the end of flank A and flank B including the microhomology.

countEvents:	the number of read pairs associated with this event
AnchorCount:	the number of distinct non-primer reads. indicates the original number of DNA fragments in your sample.
ANCHOR_DIST:	in case of a primer read with translocation, this is the distance between the start of flank B on the primer read, and end of flank b on the non-primer read. in the case of a primer read that stays on the same contig, this is the distance between the start of the primer read, and end of the non-primer read. Note that if there are deletions or insertions in the contig, then this value does not match the original fragment length. 


The output can be analysed with SIQPlotteR: https://siq.researchlumc.nl/SIQPlotteR/