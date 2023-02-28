Install:
- Install bioconda. Download the .sh file from the website www.anaconda.org and run it (bash miniconda.sh). 
- Then install samtools, picardtools, trimmomatic, and bwa-mem2 via the command "conda install -c bioconda <samtools/picard/trimmomatic/bwa-mem2>"
- Then install dos2unix with "apt install dos2unix"

Running:
First you need to create a file named "Sample_information.txt". This file needs to have one row per sample, each with the fields indicated below. Include a header row.

Sample:	The name of the library sample. Example: "7C1_LB".
Primer:	The secondary GUIDEseq primer sequence (GSP2).
File_name:	the part before "_R1.fastq.gz" or "_R2.fastq.gz".
Ref:	the genome reference file name, containing all the chromosomes, as well as the plasmid sequence(s). Example: "pCAS-PPO.fa".
P5_adapter:	Enter the full P5 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
P7_adapter:	Enter the full p7 adapter primer, without special characters. Degenerate bases allowed but will be automatically changed to N.
Plasmid:	the name of the plasmid as written in the reference file. Example: "pCAS-PPO".
DSB_chrom:	the chromosome or plasmid on which the primer sequence is located, should agree with the name in the reference fasta. Example:"Chr4".
FlankAUltEnd:	the position on the DSB_chrom where flank A ends (if there are multiple possibilities, take the furthest downstream). Flank A is the flank starting with the primer. Flank A and Flank B may overlap.
Locus_name: a name denoting the locus in which a DSB is created, or T-DNA border that is focused on.
FLANK_A_ORIENT: The orientation of the primer and thus of flank A, relative to the reference genome sequence. FW or RV.
FlankBUltStart:	the position on the DSB chrom on where flank B starts (if there are multiple possibilities, take the furthest upstream). Flank B is the flank after the nick or DSB. Flank B and flank A may overlap.
Genotype: The name of the genotype of the organism that was sampled.
Plasmid_alt: the name of the second plasmid, if a transformation with multiple plasmids has been performed. Enter "NA" if you only used 1 plasmid.
DNA: name of the DNA sample from which libraries were made. Often one DNA sample is used to do an LB, RB, FW, and RV reaction. Example: "7C1".
RunID: name of the sequencing run. This allows one to differentiate between samples that were sequenced multiple times, or to limit comparisons within runs to avoid effects of different sequencing depths between runs.

Then place the CISGUIDE_primary.sh script, Sample_information.txt, reference fastas ("PLASMID_NAME.fa"), raw sequencing files ("NAME_R1.fastq.gz" and "NAME_R2.fastq.gz") in your work directory.
Then run the CISGUIDE_primary.sh script, by doing the following:
bash CISGUIDE_primary.sh |& tee -a primary_log.txt
Example with options:
bash shared/CISGUIDE_primary.sh -p shared -f FALSE -d OPT -t 150 |& tee -a shared/primary_log.txt


The output is a folder per library, containing a .bam and .bam.bai file so that sequences may be visualized in a genome browser, and a text file ending in "_A.txt" containing all preprocessed reads.

Next, place the files that end in "_A.txt", the reference fastas, and Sample_information.txt in the CISGUIDE input folder. Then run the R script.
Change chromosome names in the fasta so that all minus signs are changed to underscores.

The output is one excel file per library showing on every row a unique event. Another excel file contains all data together.

General notes:
- The reference fasta needs to have chromosome names starting with "Chr" (e.g. "Chr4"), not just a number.
- Make sure that the reference fastas don't contain names with minus signs.


Special notes for fasta mode:

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