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
	Then install samtools, picard, trimmomatic, and bwa-mem2 via the command "conda install -c bioconda <package>". Replace <package> with the package name.
- Then install dos2unix with "apt install dos2unix"
- Place CISGUIDE_primary.sh, Sample_information.txt, fastq.gz files, and reference fastas ("PLASMID_NAME.fa") in your work directory.
- Run the software with something like: 
	bash workdir/CISGUIDE_primary2.sh -p workdir -f FALSE -t 150 |& tee -a workdir/primary_log.txt
	Replace "workdir" with the name of the folder in which the files are
- The output is a folder per library, containing a .bam and .bam.bai file so that sequences may be visualized in a genome browser, and a text file ending in "_A.txt" containing all preprocessed reads.
- Next, place the files that end in "_A.txt", the reference fastas (change minus sign to underscore), and Sample_information.txt in the CISGUIDE input folder. Then run the R script (CISGUIDE.R).
- The output is one excel file per library showing on every row a unique event. Another excel file contains all data together.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Sample_information
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
You need to create a file named "Sample_information.txt". This file needs to have one row per sample, each with the fields indicated below. Include a header row.

Sample:			The name of the library sample. Example: "7C1_LB". This name should in principle be unique. If the samples have been sequenced multiple times, they can have the same name, as long as the RunID is different.
Primer:			The secondary GUIDEseq primer sequence (GSP2). Use 'NA' in fasta mode.
File_name:		The part before "_R1.fastq.gz" or "_R2.fastq.gz". This name needs to be unique in this list.
Ref:			The genome reference file name, containing all the chromosomes, as well as the plasmid sequence(s). Example: "pCAS-PPO.fa".
P5_adapter:		Enter the full P5 adapter primer, without special characters, in upper case. Degenerate bases allowed but will be automatically changed to N. Use 'NA' in fasta mode.
P7_adapter:		Enter the full p7 adapter primer, without special characters, in upper case. Degenerate bases allowed but will be automatically changed to N. Use 'NA' in fasta mode.
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
LBSeq: 			Currently not in use
RBSeq: 			Currently not in use
Ecotype: 		The name of the ecotype of the organism that was samples. E.g. "Col-0".
R1Suffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for read 1 (the read sequenced from the p5 adapter).
R2Suffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for read 2 (the read sequenced from the p7 adapter).
UMISuffix: 		The part of the filename after File_name, but before the extention .fastq.gz, for the UMI.
Family:			Indicate here if some samples are related and therefore junctions can occur in multiple samples. Use integers starting from "1". Use "0" in all rows not to include family information. This information is used to keep some junctions which otherwise would be filtered out by the duplicate position filter. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
General notes
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- The reference fasta needs to have chromosome names starting with "Chr" (e.g. "Chr4"), not just a number.
- Make sure that the reference fastas don't contain names with minus signs. Not inside the file, nor in the name of the file.
- You can run the software with options:

-p	Set work path. default: home directory
-f	Switches to fasta mode if TRUE. default: FALSE.
-t 	Trimming length. Value indicate maximum number of nt to keep. default: 999999 (meaning no trimming is performed).

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
About the output
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In the case of microhomology, the output gives you the end of flank A and flank B including the microhomology.