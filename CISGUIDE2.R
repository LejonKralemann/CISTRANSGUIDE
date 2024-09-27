###############################################################################
#install and load packages
###############################################################################
if (require(BiocManager)==FALSE) {install.packages("BiocManager")}
if (require(Biostrings)==FALSE){BiocManager::install("Biostrings")}
if (require(stringi)==FALSE){install.packages("stringi", repos = "http://cran.us.r-project.org")}
if (require(stringdist)==FALSE){install.packages("stringdist", repos = "http://cran.us.r-project.org")}
if (require(tidyverse)==FALSE){install.packages("tidyverse", repos = "http://cran.us.r-project.org")}
if (require(openxlsx)==FALSE){install.packages("openxlsx", repos = "http://cran.us.r-project.org")}

###############################################################################
#set mode
###############################################################################
GLOBAL.FASTA_MODE = FALSE #Typically false, if TRANSGUIDE/CISGUIDE library prep and illumina sequencing has been done. TRUE if sequences from another source are being analyzed with this program.
#adjustable parameters are automatically optimized for fasta mode if FASTA_MODE==TRUE

###############################################################################
#set parameters - adjustable
###############################################################################
GLOBAL.input_dir= "./input/"
GLOBAL.output_dir= "./output/"
GLOBAL.GROUPSAMEPOS=FALSE #if true, it combines reads with the same genomic pos, which helps in removing artefacts. Typically used for TRANSGUIDE, but disabled for CISGUIDE.
GLOBAL.REMOVENONTRANS=TRUE #if true, it only considers translocations. Typically used for TRANSGUIDE, but disabled for CISGUIDE. Note that some translocations on the same chromosome will also be removed thusly.
GLOBAL.REMOVEPROBLEMS=FALSE #if true it removes all problematic reads from the combined datafile. Note if this is false, no duplicate filtering will be performed, because first reads due to barcode hopping need to be removed by removing events with few anchors. Cannot be used for CISGUIDE, because duplicate positions between samples are expected.
GLOBAL.ANCHORCUTOFF=3 #each event needs to have at least this number of anchors, otherwise it is marked as problematic (and potentially removed) 
GLOBAL.MINANCHORDIST=150 #should be matching a situation where the mate is 100% flank B (no overlap with flank A).
GLOBAL.MAXANCHORDIST=2000 #the furthest position that the mate anchor can be, except on T-DNA.
GLOBAL.FLANKBEYONDDSB=5000 #how much flank A and flank B are allowed to continue beyond the DSB (not applicable when the focus contig is the T-DNA)
GLOBAL.MINLEN=150 #this is the minimal read length. if you write NA here, then the software will calculate the minimal read length based on the distance to nick/dsb and FLANK_B_LEN_MIN. Should be at the very least 60bp, but 90bp is more common to have as minimum.
GLOBAL.LB_SEQUENCES = c("TGGCAGGATATATTGTGGTGTAAAC", "CGGCAGGATATATTCAATTGTAAAT") #the nick is made after the 3rd nt
GLOBAL.RB_SEQUENCES = c("TGACAGGATATATTGGCGGGTAAAC", "TGGCAGGATATATGCGGTTGTAATT", "TGGCAGGATATATACCGTTGTAATT") #the nick is made after the 3rd nt
GLOBAL.TD_SIZE_CUTOFF = 6 #the smallest TD that is considered as TD (*with regards to the Type variable). Any smaller TD is considered merely an insertion.
GLOBAL.TESTNAME = "GTGM0094-0027-1-003_1" #name of a read, used for testing
GLOBAL.DEBUG = FALSE #If true, only the read with GLOBAL.TESTNAME is processed
GLOBAL.UNGROUPMATES = TRUE #if true, it separates events based on the mate position. FALSE is usually recommended, except in cases where the mate position is really important.

###############################################################################
#set parameters - non-adjustable
###############################################################################
GLOBAL.NF_NUMBER = as.integer(-99999999) #don't change
GLOBAL.ERROR_NUMBER = as.integer(99999999) #don't change
GLOBAL.hash=system("git rev-parse HEAD", intern=TRUE)
GLOBAL.hash_little=substr(GLOBAL.hash, 1, 8)
GLOBAL.sample_info = read.csv(paste0(GLOBAL.input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
GLOBAL.TIME_START=round(as.numeric(Sys.time())*1000, digits=0)
GLOBAL.FLANK_B_LEN_MIN = 30 #minimum length of flank B. Also determines the size of DSB_AREA_SEQ. do not change because the preprocessing program will still be set at 30.
GLOBAL.MINMAPQUALA = 42 #minimum mapping quality (phred). 42 means a perfect, unambiguous match (well, should be)
GLOBAL.MINMAPQUALB = 42 #same, but for flank B
GLOBAL.PercentageDone = 0 #var for indicating progress
GLOBAL.TotalFileSize = 0
#automatically adjust parameters if fasta mode is on
if (GLOBAL.FASTA_MODE==TRUE){
  GLOBAL.ANCHORCUTOFF=1 #each event needs to have at least this number of anchors, otherwise it is marked as problematic (and potentially removed) 
  GLOBAL.MINANCHORDIST=0 #should be matching a situation where the mate is 100% flank B (no overlap with flank A).
  GLOBAL.MAXANCHORDIST=10000 #the furthest position that the mate anchor can be, except on T-DNA.
}

###############################################################################
#Initial checks
###############################################################################
runlog=c()
if (file.exists(paste0(GLOBAL.input_dir, "Sample_information.txt"))==FALSE){
  runlog=rbind(runlog, "Sample information sheet not present, aborting")
  message(runlog[length(runlog),])
  quit()}
if (file.exists(paste0(GLOBAL.input_dir, "read_numbers.txt"))==FALSE){
  runlog=rbind(runlog, "Read numbers file not found, aborting")
  message(runlog[length(runlog),])
  quit()}

###############################################################################
#Functions
###############################################################################
#function: matcher_skipper
#this function matches from the end until it encounters a mismatch. If the continuing sequence is longer than 9, then it jumps over this mismatch, and continues matching until a second mismatch is encountered. Then it outputs the matching sequence.  
matcher_skipper <-function(ref, seq1){
  if (ref != "ERROR"){
    seq1_len = nchar(seq1)                                                          #get the length of seq1
    ref_len = nchar(ref)                                                            #get the length of ref
    match_1_len = lcprefix(ref, seq1)                                               #find a match between ref and read
    match_1_non_ref = if ((match_1_len + 2)< ref_len){                              #if match length + 2 bp is smaller than ref length
      substr(ref, start = match_1_len+2, stop = ref_len)                            #get the non-matching sequence from the ref, minus the first mismatch
    }else{
      ""}                                                                           #otherwise output empty string
    match_1_non_seq1 = if (match_1_non_ref != ""){                                  #if the match sequence is not empty
      substr(seq1, start = match_1_len +2, stop = seq1_len)                         #get the non-matching sequence from the read, minus the first mismatch
    }else{
      ""}                                                                           #otherwise output empty string
    match_2_len = if (match_1_non_ref != ""){                                       #if the match sequence is not empty
      lcprefix(match_1_non_ref, match_1_non_seq1)                                   #find a match between ref_match and read_non_match
    }else{          
      0}                                                                            #otherwise output 0 as length
    subtotal_match = if (match_2_len > 9){                                          #if the match 2 length is larger than 9 bp
      substr(ref, start = 1, stop = match_1_len + 1 + match_2_len)                  #get the subtotal matching sequence
    }else{
      substr(ref, start = 1, stop = match_1_len)}                                   #else get the primary matching sequence only
  
    subtotal_match_non_ref = if(match_2_len >9 & nchar(subtotal_match)+9 < ref_len){#if there is a secondary match, and we didn't go up to the end of the ref yet
      substr(ref, start = match_1_len + match_2_len + 3, stop = seq1_len)           #get the ref sequence minus the subtotal matching sequence, minus the first mismatch
    }else{
      ""}                                                                           #otherwise get an empty string
    subtotal_match_non_seq1 = if(subtotal_match_non_ref != ""){                     #if there is a ref sequence to match with
      substr(seq1, start = match_1_len + match_2_len + 3, stop = seq1_len)          #get the read sequence minus the subtotal matching sequence, minus the first mismatch
    }else{
    ""}                                                                             #otherwise get an empty string
    match_3_len = if (subtotal_match_non_ref != ""){                                #if there is a ref sequence to match with
      lcprefix(subtotal_match_non_ref, subtotal_match_non_seq1)                     #match the pieces of ref and seq1
    }else{
      0}                                                                            #otherwise output 0 as length
    flank_match = if (match_3_len > 9){                                             #if the match 3 is longer than 9bp
      substr(ref, start = 1, stop = match_1_len + 1 + match_2_len + 1 + match_3_len)#get the total matching sequence with match 1, 2 and 3
    }else{
      subtotal_match}                                                               #else only get the matching sequence with match 1 and 2
    return(flank_match)                                                             #return the total match
  }
}
funlog <-function(warningtext){
  runlog<<-rbind(runlog, warningtext)
  message(runlog[length(runlog),])
}
function_time <-function(text){
  TIME_CURRENT=round(as.numeric(Sys.time())*1000, digits=0)
  funlog(paste0(text, (TIME_CURRENT - GLOBAL.TIME_START), " milliseconds"))
  GLOBAL.TIME_START<<-round(as.numeric(Sys.time())*1000, digits=0)
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###############################################################################
#Process data: step 0
#calculating total work
###############################################################################

funlog("calculating total work")

for (i in row.names(GLOBAL.sample_info)){
  FILE.Sample = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Sample))
  FILE.RunID = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(RunID))
  
  if (file.exists(paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))==FALSE){
    next
  }else{
    GLOBAL.TotalFileSize = GLOBAL.TotalFileSize + (file.info((paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))))$size
  }
}

funlog(paste0("Total file size to process: ", GLOBAL.TotalFileSize, " bytes"))

###############################################################################
#Process data: step 1
#checking the file and reading metadata
###############################################################################

funlog("Checking file and reading metadata")

#process all files
for (i in row.names(GLOBAL.sample_info)){
  
  ####################  pre-cleanup  #####################
  objects_to_clean=c("FILE.AgroGeno", "FILE.contig_seq", "FILE.CurrentFileSize", "FILE.DNASample", "FILE.DSB_AREA_SEQ", "FILE.DSB_AREA_SEQ_RC", "FILE.DSB_CONTIG", "FILE.DSB_FW_END", "FILE.Ecotype", "FILE.FLANK_A_ISFORWARD", "FILE.FLANK_A_REF", "FILE.FlankAUltEnd", "FILE.FOCUS_CONTIG", "FILE.FOCUS_LOCUS", "FILE.Genotype", "FILE.LOCUS_NAME", "FILE.PLASMID", "FILE.PLASMID_ALT", "FILE.plasmid_seq", "FILE.Primer_match_perfect", "FILE.Primer_on_TDNA", "FILE.Primer_pos", "FILE.Primer_seq", "FILE.PRIMER_TO_DSB", "FILE.REF", "FILE.RunID", "FILE.Sample", "FILE.TDNA_ALT_IS_LBRB", "FILE.TDNA_ALT_LB_END", "FILE.TDNA_ALT_LB_FW", "FILE.TDNA_ALT_RB_END", "FILE.TDNA_ALT_RB_FW", "FILE.TDNA_IS_LBRB", "FILE.TDNA_LB_END", "FILE.TDNA_LB_FW", "FILE.TDNA_RB_END", "FILE.TDNA_RB_FW", "FILE.TOTAL_REF", "FILE.TOTAL_REF_START", "FILE.TOTAL_REF_STOP", "FILE.data", "FILE.data2", "FILE.data3", "FILE.data4", "FILE.data5", "FILE.data6", "FILE.data7", "FILE.data8", "FILE.data9", "FILE.data10", "FILE.data11", "FILE.data12", "FILE.data13", "FILE.data14", "FILE.data15", "FILE.data16", "FILE.data17", "FILE.data18", "FILE.genomeseq", "FILE.LB_match", "FILE.LB_match_RV", "FILE.LB2_match", "FILE.LB2_match_RV", "FILE.Primer_match", "FILE.Primer_RC_match", "FILE.RB_match", "FILE.RB_match_RV", "FILE.RB2_match", "FILE.RB2_match_RV")
  for (z in objects_to_clean){
    assign(z, NULL)
  }
 
  ####################  general variables acquired from the information sheet  #####################
  
  FILE.Sample = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Sample))
  FILE.RunID = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(RunID))
  
  ####################  check for existence of files  #####################
  
  if (file.exists(paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))==FALSE){
    funlog(paste0("Primary processed file ",FILE.Sample, "_", FILE.RunID, "_A.txt not found, moving to the next sample"))
    next
  }else if (file.exists(paste0(GLOBAL.output_dir, FILE.Sample, "_", FILE.RunID, "_CISTRANSGUIDE_V2.xlsx"))==TRUE){   #check whether file has already been processed
    funlog(paste0("File ", GLOBAL.output_dir, FILE.Sample, "_", FILE.RunID, "_A.txt has already been processed, moving to the next sample"))
    #show progress
    FILE.CurrentFileSize = (file.info((paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))))$size
    GLOBAL.PercentageDone = GLOBAL.PercentageDone + ((FILE.CurrentFileSize/GLOBAL.TotalFileSize)*99)
    funlog(paste0("CISTRANSGUIDE analysis ", round(GLOBAL.PercentageDone, digits=3), "% complete"))
    next
    }else{
      funlog(paste0("Processing ",GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))
      FILE.CurrentFileSize = (file.info((paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))))$size
    }
  
  FILE.data = read.csv(paste0(GLOBAL.input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
  if (nrow(FILE.data)==0){
    funlog("Primary processed file empty, moving to the next sample")
    next
  }
  
  ####################  continue general variables acquired from the information sheet  #####################
  
  FILE.DSB_CONTIG = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(DSB_CONTIG))#chromosome name or NA
  if (is.na(FILE.DSB_CONTIG) | FILE.DSB_CONTIG=="NA"){FILE.DSB_CONTIG=""}
  FILE.Genotype = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Genotype))
  FILE.PLASMID = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Plasmid))
  FILE.PLASMID_ALT = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Plasmid_alt))
  if (is.na(FILE.PLASMID_ALT) | FILE.PLASMID_ALT=="NA"){FILE.PLASMID_ALT=""} #to avoid doing the following check: FLANK_B_CHROM==NA
  
  FILE.REF = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Ref))
  FILE.DNASample = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(DNA))
  FILE.Ecotype = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Ecotype))
  FILE.AgroGeno = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(AgroGeno))
  FILE.DSB_FW_END = as.integer(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(DSB_FW_END)) #end of left flank before DSB
  if (is.na(FILE.DSB_FW_END) | FILE.DSB_FW_END=="NA"){FILE.DSB_FW_END=0}
  FILE.FOCUS_CONTIG = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Focus_contig_name))#Same as DSB_contig, or same as plasmid
  FILE.LOCUS_NAME = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Locus_name))#name of the genomic locus in which a DSB is made
  if (is.na(FILE.LOCUS_NAME) | FILE.LOCUS_NAME=="NA"){FILE.LOCUS_NAME=""}
  FILE.Primer_seq = str_replace_all(toupper(as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Primer))), "TCAGACGTGTGCTCTTCCGATCT", "")
  FILE.Species = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Species))
  FILE.Experiment = as.character(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(Experiment))
  
  #the following variables can be supplied via the sample information sheet, but NA is also allowed, then the software will look.
  FILE.TDNA_LB_END = as.integer(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(TDNA_LB_END))
  FILE.TDNA_RB_END = as.integer(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(TDNA_RB_END))
  FILE.TDNA_ALT_LB_END = as.integer(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(TDNA_ALT_LB_END))
  FILE.TDNA_ALT_RB_END = as.integer(GLOBAL.sample_info %>% filter(row.names(GLOBAL.sample_info) %in% i) %>% select(TDNA_ALT_RB_END))
  
  
  
  ####################  calculated general variables  #####################
  
  ############# REF CHECK ##############
  
  if (file.exists(paste0(GLOBAL.input_dir, FILE.REF))==FALSE){
    funlog("Reference fasta not found. Moving to next sample.")
    next
  }else{
    funlog(paste0("Using ref ", GLOBAL.input_dir, FILE.REF))
  }
  
  FILE.genomeseq = readDNAStringSet(paste0(GLOBAL.input_dir, FILE.REF) , format="fasta")
  
  ############ FINDING LB, RB ############################
  
  FILE.plasmid_seq = as.character(eval(parse(text = paste0("FILE.genomeseq$`", FILE.PLASMID, "`"))))
  
  #emptying some variables
  FILE.LB_match = NULL
  FILE.LB_match_RV = NULL
  FILE.LB2_match = NULL
  FILE.LB2_match_RV = NULL
  FILE.RB_match = NULL
  FILE.RB_match_RV = NULL
  FILE.RB2_match = NULL
  FILE.RB2_match_RV = NULL
  
  #check whether plasmid can be found in ref
  if (isEmpty(FILE.plasmid_seq)==TRUE){
   funlog(paste0("Plasmid name ", FILE.PLASMID, " not found in ", FILE.REF, " . Moving to next sample."))
    next
  }
  
  
  #find the LB
  if (is.na(FILE.TDNA_LB_END)){
    FILE.LB_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[1]], subject = DNAString(FILE.plasmid_seq), max.mismatch = 0, fixed=TRUE))
    FILE.LB_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[1]]))), subject = DNAString(FILE.plasmid_seq), max.mismatch = 0, fixed=TRUE))
    FILE.LB2_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[2]],subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.LB2_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[2]]))),subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.TDNA_LB_END = NA
    FILE.TDNA_LB_FW = NA
    
    for (i in c("FILE.LB_match", "FILE.LB_match_RV", "FILE.LB2_match", "FILE.LB2_match_RV")){
      if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
        
        
        if (i=="FILE.LB_match" | i=="FILE.LB2_match"){
        FILE.TDNA_LB_FW = TRUE
        FILE.TDNA_LB_END = as.numeric(get(i)$start)+3 
        }else{
        FILE.TDNA_LB_FW = FALSE
        FILE.TDNA_LB_END = as.numeric(get(i)$start)+21 
        }
      } else {
        next
      } 
    }
    if (is.na(FILE.TDNA_LB_END)){
      funlog("No single LB sequence found, moving to next sample")
      next
    }
  }else{
    funlog("Using user supplied LB position")
  }
    
    #find the RB
  if (is.na(FILE.TDNA_RB_END)){
    FILE.RB_match = as.data.frame(matchPattern(pattern = GLOBAL.RB_SEQUENCES[[1]], subject = DNAString(FILE.plasmid_seq), max.mismatch = 0, fixed=TRUE))
    FILE.RB_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.RB_SEQUENCES[[1]]))), subject = DNAString(FILE.plasmid_seq), max.mismatch = 0, fixed=TRUE))
    FILE.RB2_match = as.data.frame(matchPattern(pattern = GLOBAL.RB_SEQUENCES[[2]],subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.RB2_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.RB_SEQUENCES[[2]]))),subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.RB3_match = as.data.frame(matchPattern(pattern = GLOBAL.RB_SEQUENCES[[3]],subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.RB3_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.RB_SEQUENCES[[3]]))),subject = DNAString(FILE.plasmid_seq),max.mismatch = 0,fixed = TRUE))
    FILE.TDNA_RB_END = NA
    FILE.TDNA_RB_FW = NA
    
    for (i in c("FILE.RB_match", "FILE.RB_match_RV", "FILE.RB2_match", "FILE.RB2_match_RV", "FILE.RB3_match", "FILE.RB3_match_RV")){
      if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
        
        
        if (i=="FILE.RB_match" | i=="FILE.RB2_match" | i=="FILE.RB3_match"){
          FILE.TDNA_RB_FW = TRUE
          FILE.TDNA_RB_END = as.numeric(get(i)$start)+2 
        }else{
          FILE.TDNA_RB_FW = FALSE
          FILE.TDNA_RB_END = as.numeric(get(i)$start)+22 
        }
      } else {
        next
      } 
    }
    if (is.na(FILE.TDNA_RB_END)){
      funlog("No single RB sequence found, moving to next sample")
      next
    }
  }else{
    funlog("Using user supplied RB position")
  }
    
    #get the LB and RB positions of the alternative plasmid
  FILE.TDNA_ALT_LB_END = NA
  FILE.TDNA_ALT_LB_FW = NA
  FILE.TDNA_ALT_RB_END = NA
  FILE.TDNA_ALT_RB_FW = NA
    
    if (FILE.PLASMID_ALT != ""){
      FILE.plasmid_alt_seq = as.character(eval(parse(text = paste0("FILE.genomeseq$`", FILE.PLASMID_ALT, "`"))))
      
      #check whether alternative plasmid can be found in ref
      if (isEmpty(FILE.plasmid_alt_seq)==TRUE){
        funlog(paste0("Plasmid name ", FILE.PLASMID_ALT, " not found in ", FILE.REF, " . Moving to next sample."))
        next
      }
      
      #find the LB
      FILE.LB_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[1]], subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.LB_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[1]]))), subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.LB2_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[2]],subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      FILE.LB2_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[2]]))),subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      
      
      for (i in c("FILE.LB_match", "FILE.LB_match_RV", "FILE.LB2_match", "FILE.LB2_match_RV")){
        if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
          
          
          if (i=="FILE.LB_match" | i=="FILE.LB2_match"){
            FILE.TDNA_ALT_LB_FW = TRUE
            FILE.TDNA_ALT_LB_END = as.numeric(get(i)$start)+3 
          }else{
            FILE.TDNA_ALT_LB_FW = FALSE
            FILE.TDNA_ALT_LB_END = as.numeric(get(i)$start)+21 
          }
        } else {
          next
        } 
      }
      if (is.na(FILE.TDNA_ALT_LB_END)){
        funlog("No single LB sequence found, moving to next sample")
        next
      }
      
      #find the RB
      FILE.RB_match = as.data.frame(matchPattern(pattern = GLOBAL.RB_SEQUENCES[[1]], subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.RB_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.RB_SEQUENCES[[1]]))), subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.RB2_match = as.data.frame(matchPattern(pattern = GLOBAL.RB_SEQUENCES[[2]],subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      FILE.RB2_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.RB_SEQUENCES[[2]]))),subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      
      
      for (i in c("FILE.RB_match", "FILE.RB_match_RV", "FILE.RB2_match", "FILE.RB2_match_RV")){
        if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
          
          
          if (i=="FILE.RB_match" | i=="FILE.RB2_match"){
            FILE.TDNA_ALT_RB_FW = TRUE
            FILE.TDNA_ALT_RB_END = as.numeric(get(i)$start)+2 
          }else{
            FILE.TDNA_ALT_RB_FW = FALSE
            FILE.TDNA_ALT_RB_END = as.numeric(get(i)$start)+22 
          }
        } else {
          next
        } 
      }
      if (is.na(FILE.TDNA_ALT_RB_END)){
        funlog("No single RB sequence found, moving to next sample")
        next
      }
    }
    
  #determine the orientation of the T-DNA on the main and alternative plasmids
  if (FILE.TDNA_LB_END < FILE.TDNA_RB_END){
    FILE.TDNA_IS_LBRB = TRUE
  }else{
    FILE.TDNA_IS_LBRB = FALSE
  }
    
  if (FILE.PLASMID_ALT != ""){
    if (FILE.TDNA_ALT_LB_END < FILE.TDNA_ALT_RB_END){
      FILE.TDNA_ALT_IS_LBRB = TRUE
    }else{
      FILE.TDNA_ALT_IS_LBRB = FALSE}
  }else{
    FILE.TDNA_ALT_IS_LBRB = NA
  }
  
  ####################  REF check  #####################
  
  
  FILE.contig_seq = as.character(eval(parse(text = paste0("FILE.genomeseq$`", FILE.FOCUS_CONTIG, "`"))))
  if (length(FILE.contig_seq)==0){
    funlog("Focus contig not found in reference fasta. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
    next 
  }
  
  ####################  Primer check  #####################
  
  if (GLOBAL.FASTA_MODE==FALSE){  
    
    #this checks whether primer can be found, checks what the orientation of flank A is, whether the primer matches the plasmid or not
    FILE.Primer_match = as.data.frame(matchPattern(pattern = FILE.Primer_seq,subject = DNAString(FILE.contig_seq),max.mismatch = 0,fixed = TRUE))
    FILE.Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(FILE.Primer_seq))),subject = DNAString(FILE.contig_seq),max.mismatch = 0,fixed = TRUE))
    
    #in case of a T-DNA primer
    if (FILE.FOCUS_CONTIG == FILE.PLASMID) {
      if (nrow(FILE.Primer_match) > 0 & nrow(FILE.Primer_match) < 2) {
        FILE.FLANK_A_ISFORWARD = TRUE
        FILE.Primer_pos = as.numeric(FILE.Primer_match$start)
        if (FILE.TDNA_IS_LBRB == TRUE & FILE.Primer_pos > FILE.TDNA_LB_END & FILE.Primer_pos < FILE.TDNA_RB_END) {
          FILE.Primer_on_TDNA = TRUE
          FILE.FOCUS_LOCUS = "RB"
          FILE.FlankAUltEnd = FILE.TDNA_RB_END
        }else if (FILE.TDNA_IS_LBRB == FALSE & FILE.Primer_pos < FILE.TDNA_LB_END & FILE.Primer_pos > FILE.TDNA_RB_END){
          FILE.Primer_on_TDNA = TRUE
          FILE.FOCUS_LOCUS = "LB"
          FILE.FlankAUltEnd = FILE.TDNA_LB_END
        }else{
          funlog("Primer not within the T-DNA, moving to the next sample")
          next
        }
      } else  if (nrow(FILE.Primer_RC_match) > 0 & nrow(FILE.Primer_RC_match) < 2) {
        FILE.FLANK_A_ISFORWARD = FALSE
        FILE.Primer_pos = as.numeric(FILE.Primer_RC_match$end)
        if (FILE.TDNA_IS_LBRB == TRUE & FILE.Primer_pos > FILE.TDNA_LB_END & FILE.Primer_pos < FILE.TDNA_RB_END) {
          FILE.Primer_on_TDNA = TRUE
          FILE.FOCUS_LOCUS = "LB"
          FILE.FlankAUltEnd = FILE.TDNA_LB_END
        }else if (FILE.TDNA_IS_LBRB == FALSE & FILE.Primer_pos < FILE.TDNA_LB_END & FILE.Primer_pos > FILE.TDNA_RB_END) {
          FILE.Primer_on_TDNA = TRUE
          FILE.FOCUS_LOCUS = "RB"
          FILE.FlankAUltEnd = FILE.TDNA_RB_END
        }else{
          funlog("Primer not within the T-DNA, moving to the next sample")
          next
        }
      } else if (nrow(FILE.Primer_match) > 1 | nrow(FILE.Primer_RC_match) > 1) {
        funlog("Primer found several times on this contig, moving to the next sample")
        next
      } else{
        funlog("Primer not found, moving to the next sample")
        next
      }
      #if the primer is located on a chromosome (next to an induced DSB)
    } else{
      FILE.FOCUS_LOCUS = FILE.LOCUS_NAME
      FILE.Primer_on_TDNA = FALSE
      
      if (nrow(FILE.Primer_match) > 0 & nrow(FILE.Primer_match) < 2) {
        FILE.FLANK_A_ISFORWARD = TRUE
        FILE.FlankAUltEnd = FILE.DSB_FW_END
        FILE.Primer_pos = as.numeric(FILE.Primer_match$start)
      } else if (nrow(FILE.Primer_RC_match) > 0 & nrow(FILE.Primer_RC_match) < 2) {
        FILE.FLANK_A_ISFORWARD = FALSE
        FILE.FlankAUltEnd = FILE.DSB_FW_END + 1
        FILE.Primer_pos = as.numeric(FILE.Primer_RC_match$end)
      } else if (nrow(FILE.Primer_match) > 1 | nrow(FILE.Primer_RC_match) > 1) {
        funlog("Primer found several times on this contig, moving to the next sample")
        next
      } else{
        funlog("Primer not found, moving to the next sample")
        next
      }
    }
    
  ####################  acquire the DSB AREA SEQ  #####################
  
 
  FILE.DSB_AREA_SEQ = if (FILE.DSB_CONTIG == FILE.FOCUS_CONTIG){ (if (FILE.FLANK_A_ISFORWARD == TRUE){
    substr(FILE.contig_seq, start= FILE.FlankAUltEnd - ((GLOBAL.FLANK_B_LEN_MIN/2)-1), stop= FILE.FlankAUltEnd + (GLOBAL.FLANK_B_LEN_MIN/2))
    }else if (FILE.FLANK_A_ISFORWARD == FALSE){
    as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= FILE.FlankAUltEnd - (GLOBAL.FLANK_B_LEN_MIN/2), stop= FILE.FlankAUltEnd + ((GLOBAL.FLANK_B_LEN_MIN/2)-1)))))
      })
  }else{
    ""
  }

  FILE.DSB_AREA_SEQ_RC = as.character(reverseComplement(DNAString(FILE.DSB_AREA_SEQ)))
    
  ####################  continue calculated general variables  #####################  
    
  #calculate the distance from primer to DSB/nick
  FILE.PRIMER_TO_DSB = if (FILE.FLANK_A_ISFORWARD == TRUE){
    FILE.FlankAUltEnd - (FILE.Primer_pos -1)
  }else{
    FILE.Primer_pos - (FILE.FlankAUltEnd -1)
  }

  #usually you don't want the primer to be so far away. But sometimes when you cut away the end of the T-DNA for instance, then you may want to keep the T-DNA end position.
  if (FILE.PRIMER_TO_DSB>300){
    funlog("Warning: primer is more than 300 bp away from the indicated end of FLANK A. Continuing anyway.")
  }

 
  #get the start and stop positions for the ref sequence, in case it goes beyond the end of the chromosome
  FILE.TOTAL_REF_START = if(FILE.FlankAUltEnd-GLOBAL.FLANKBEYONDDSB < 1){
    1
  }else{
    FILE.FlankAUltEnd-GLOBAL.FLANKBEYONDDSB
  }
  FILE.TOTAL_REF_STOP = if(FILE.FlankAUltEnd+GLOBAL.FLANKBEYONDDSB > nchar(FILE.contig_seq)){
    nchar(FILE.contig_seq)
  }else{
    FILE.FlankAUltEnd+GLOBAL.FLANKBEYONDDSB
  }
  
  #get the REF seq for flank A.
  FILE.FLANK_A_REF = if (FILE.FLANK_A_ISFORWARD == TRUE){
    substr(FILE.contig_seq, start= FILE.Primer_pos, stop= FILE.TOTAL_REF_STOP)
  }else{
    as.character(reverseComplement(DNAString(
      substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.Primer_pos))
    ))
  }
  #get the REF seq that includes both flank a and b, later used to get the ref for flank b  
  FILE.TOTAL_REF = if (FILE.FLANK_A_ISFORWARD == TRUE){substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.TOTAL_REF_STOP)
  }else{
    as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.TOTAL_REF_STOP))))
  }
  }
  
  function_time("Step 1 took ")

  ###############################################################################
  #Process data: step 2
  ###############################################################################
  
  if (GLOBAL.FASTA_MODE == TRUE){
    
    FILE.data1 = FILE.data  %>%
    #first determine if the read is T-DNA:genome or genome:T-DNA. If genome:T-DNA, reverse the sequence and change the B_POS
    mutate(fasta_orient = case_when(A_CHROM == FILE.FOCUS_CONTIG & FLANK_B_CHROM != FILE.FOCUS_CONTIG ~ "TDNA-genome",
                                    A_CHROM != FILE.FOCUS_CONTIG & FLANK_B_CHROM == FILE.FOCUS_CONTIG ~ "genome-TDNA",
                                    TRUE ~ "Other"))%>%
    rowwise()%>%
    mutate(SEQ_1_NEW = if_else(fasta_orient == "genome-TDNA",
                               SEQ_2,
                               SEQ_1))%>%
    mutate(SEQ_2_NEW = if_else(fasta_orient == "genome-TDNA",
                                 SEQ_1,
                                 SEQ_2))%>%
    mutate(B_POS_NEW = if_else(fasta_orient == "genome-TDNA",
                               A_POS,
                               B_POS))%>%
    mutate(FLANK_B_ORIENT_NEW = case_when(fasta_orient == "genome-TDNA" & A_ORIENT == "FW" ~ "RV",
                          fasta_orient == "genome-TDNA" & A_ORIENT == "RV" ~ "FW",
                          TRUE ~ FLANK_B_ORIENT))%>%
    mutate(FLANK_B_CHROM_NEW = if_else(fasta_orient == "genome-TDNA",
                                       A_CHROM,
                                       FLANK_B_CHROM))%>%
    mutate(MATE_FLANK_B_CHROM_NEW = if_else(fasta_orient == "genome-TDNA",
                                            A_CHROM,
                                            MATE_FLANK_B_CHROM))%>%
    mutate(MATE_B_ORIENT_NEW = if_else(FLANK_B_ORIENT_NEW == "FW",
                                       "RV",
                                       "FW"))%>%
    mutate(MATE_B_POS_NEW = B_POS_NEW)%>%
    mutate(SEQ_1 = SEQ_1_NEW,
           SEQ_2 = SEQ_2_NEW,
           B_POS = B_POS_NEW,
           FLANK_B_ORIENT = FLANK_B_ORIENT_NEW,
           FLANK_B_CHROM = FLANK_B_CHROM_NEW,
           MATE_FLANK_B_CHROM = MATE_FLANK_B_CHROM_NEW,
           MATE_B_ORIENT = MATE_B_ORIENT_NEW,
           MATE_B_POS = MATE_B_POS_NEW)%>%
    #make the first 30 bp a fake primer, and find where the primer matches  
    rowwise()%>%
    mutate(PRIMER_SEQ = substr(SEQ_1, 1, 30))%>%
    #initialize columns
    mutate(Primer_pos = GLOBAL.ERROR_NUMBER,
           Primer_OK = FALSE,
           FLANK_A_ISFORWARD = FALSE,
           FlankAUltEnd = GLOBAL.ERROR_NUMBER,
           FOCUS_LOCUS = "",
           Primer_on_TDNA = FALSE)
    
    for (j in row.names(FILE.data1)){
      j_int=as.integer(j)
      Primer_match = as.data.frame(matchPattern(pattern = FILE.data1$PRIMER_SEQ[[j_int]], subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE))
      Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(FILE.data1$PRIMER_SEQ[[j_int]]))), subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE))
    if (FILE.FOCUS_CONTIG == FILE.PLASMID){
      if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
        FILE.data1$FLANK_A_ISFORWARD[[j_int]]=TRUE
        FILE.data1$Primer_pos[[j_int]] = as.numeric(Primer_match$start)
        if (FILE.TDNA_IS_LBRB == TRUE & FILE.data1$Primer_pos[[j_int]] > FILE.TDNA_LB_END & FILE.data1$Primer_pos[[j_int]] < FILE.TDNA_RB_END){
          FILE.data1$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data1$FOCUS_LOCUS[[j_int]]="RB"
          FILE.data1$FlankAUltEnd[[j_int]] = FILE.TDNA_RB_END
          FILE.data1$Primer_OK[[j_int]]=TRUE
        }else if (FILE.TDNA_IS_LBRB == FALSE & FILE.data1$Primer_pos[[j_int]] < FILE.TDNA_LB_END & FILE.data1$Primer_pos[[j_int]] > FILE.TDNA_RB_END){
          FILE.data1$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data1$FOCUS_LOCUS[[j_int]]="LB"
          FILE.data1$FlankAUltEnd[[j_int]] = FILE.TDNA_LB_END
          FILE.data1$Primer_OK[[j_int]]=TRUE
        }else{
          #primer not on T-DNA
          if (GLOBAL.DEBUG==TRUE){
            funlog("Warning: primer not on T-DNA. Continuing anyway.")
          }
          FILE.data1$Primer_OK[[j_int]]=FALSE
        }
        
      }else  if (nrow(Primer_RC_match) > 0 & nrow(Primer_RC_match) < 2){
        FILE.data1$FLANK_A_ISFORWARD[[j_int]]=FALSE
        FILE.data1$Primer_pos[[j_int]] = as.numeric(Primer_RC_match$end)
        if (FILE.TDNA_IS_LBRB == TRUE & FILE.data1$Primer_pos[[j_int]] > FILE.TDNA_LB_END & FILE.data1$Primer_pos[[j_int]] < FILE.TDNA_RB_END){
          FILE.data1$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data1$FOCUS_LOCUS[[j_int]]="LB"
          FILE.data1$FlankAUltEnd[[j_int]] = FILE.TDNA_LB_END
          FILE.data1$Primer_OK[[j_int]]=TRUE
        }else if (FILE.TDNA_IS_LBRB == FALSE & FILE.data1$Primer_pos[[j_int]] < FILE.TDNA_LB_END & FILE.data1$Primer_pos[[j_int]] > FILE.TDNA_RB_END){
          FILE.data1$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data1$Primer_OK[[j_int]]=TRUE
          FILE.data1$FOCUS_LOCUS[[j_int]]="RB"
          FILE.data1$FlankAUltEnd[[j_int]] = FILE.TDNA_RB_END
        }else{
          #primer not on T-DNA
          if (GLOBAL.DEBUG==TRUE){
            funlog("Warning: primer not on T-DNA. Continuing anyway.")
          }
          FILE.data1$Primer_OK[[j_int]]=FALSE
        }
        

      }else if (nrow(Primer_match) >1 | nrow(Primer_RC_match) >1){
        #primer found multiple times
        if (GLOBAL.DEBUG==TRUE){
        funlog("Warning: primer found multiple times. Continuing anyway.")
        }
        FILE.data1$Primer_OK[[j_int]]=FALSE
      }else{
        #primer not found
        if (GLOBAL.DEBUG==TRUE){
          funlog("Warning: primer not found. Continuing anyway.")
        }
        FILE.data1$Primer_OK[[j_int]]=FALSE
      }
    }else{
      FILE.data1$FOCUS_LOCUS[[j_int]]=FILE.LOCUS_NAME
      FILE.data1$Primer_on_TDNA[[j_int]]=FALSE
      
      if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
        FILE.data1$Primer_pos[[j_int]] = as.numeric(Primer_match$start)
        FILE.data1$FLANK_A_ISFORWARD[[j_int]]=TRUE
        FILE.data1$FlankAUltEnd[[j_int]] = FILE.DSB_FW_END
        FILE.data1$Primer_OK[[j_int]]=TRUE
      }else if (nrow(Primer_RC_match) > 0 & nrow(Primer_RC_match) < 2){
        FILE.data1$Primer_pos[[j_int]] = as.numeric(Primer_RC_match$end)
        FILE.data1$FLANK_A_ISFORWARD[[j_int]]=FALSE
        FILE.data1$FlankAUltEnd[[j_int]] = FILE.DSB_FW_END+1
        FILE.data1$Primer_OK[[j_int]]=TRUE
      }else if (nrow(Primer_match) >1 | nrow(Primer_RC_match) >1){
        #primer found multiple times
        if (GLOBAL.DEBUG==TRUE){
          funlog("Warning: primer found multiple times. Continuing anyway.")
        }
        FILE.data1$Primer_OK[[j_int]]=FALSE
      }else{
        #primer not found
        if (GLOBAL.DEBUG==TRUE){
          funlog("Warning: primer not found. Continuing anyway.")
        }
        FILE.data1$Primer_OK[[j_int]]=FALSE
      }
    }
    }
    FILE.data1_filter = FILE.data1  %>%
      #remove rows from fasta_mode == TRUE where the first bit matches the plasmid but not T-DNA, or not matches anything, or matches multiple things, or matches genome.
      filter(Primer_OK == TRUE & fasta_orient != "Other")
    
    #check whether any reads survived
    if (nrow(FILE.data1_filter)==0){
      funlog(paste0("No reads surviving for sample ", FILE.Sample))
      #show progress
      GLOBAL.PercentageDone = GLOBAL.PercentageDone + ((FILE.CurrentFileSize/GLOBAL.TotalFileSize)*100)
      funlog(paste0("CISTRANSGUIDE analysis ", round(GLOBAL.PercentageDone, digits=3), "% complete"))
      next
    }
    
    FILE.data2 = FILE.data1_filter  %>%
      rowwise()%>%
      mutate(DSB_AREA_SEQ = case_when(FILE.DSB_CONTIG == FILE.FOCUS_CONTIG & FLANK_A_ISFORWARD == TRUE ~ substr(FILE.contig_seq, start= FlankAUltEnd - ((GLOBAL.FLANK_B_LEN_MIN/2)-1), stop= FlankAUltEnd + (GLOBAL.FLANK_B_LEN_MIN/2)),
                                      FILE.DSB_CONTIG == FILE.FOCUS_CONTIG & FLANK_A_ISFORWARD == FALSE ~ as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= FlankAUltEnd - (GLOBAL.FLANK_B_LEN_MIN/2), stop= FlankAUltEnd + ((GLOBAL.FLANK_B_LEN_MIN/2)-1))))),
                                      TRUE ~ ""
               )) %>%
      mutate(DSB_AREA_SEQ_RC = as.character(reverseComplement(DNAString(DSB_AREA_SEQ)))) %>%
      mutate(PRIMER_TO_DSB = case_when( FLANK_A_ISFORWARD == TRUE ~   FlankAUltEnd - (Primer_pos -1),
                                          TRUE ~ Primer_pos - (FlankAUltEnd -1))) %>%
      mutate(TOTAL_REF_START = case_when(FlankAUltEnd-GLOBAL.FLANKBEYONDDSB < 1 ~ 1,
                                         TRUE ~ FlankAUltEnd-GLOBAL.FLANKBEYONDDSB)) %>%
      mutate(TOTAL_REF_STOP = case_when(FlankAUltEnd+GLOBAL.FLANKBEYONDDSB > nchar(FILE.contig_seq) ~ nchar(FILE.contig_seq),
                                        TRUE ~ FlankAUltEnd+GLOBAL.FLANKBEYONDDSB)) %>%
      
      mutate(FLANK_A_REF = case_when(FLANK_A_ISFORWARD == TRUE ~ substr(FILE.contig_seq, start= Primer_pos, stop= TOTAL_REF_STOP),
                                     TRUE ~ as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= TOTAL_REF_START, stop= Primer_pos))))))%>%
    mutate(TOTAL_REF = case_when(FLANK_A_ISFORWARD == TRUE ~ substr(FILE.contig_seq, start= TOTAL_REF_START, stop= TOTAL_REF_STOP),
                                 TRUE ~ as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= TOTAL_REF_START, stop= TOTAL_REF_STOP)))))) %>%
      ungroup()
    
  }else{
    FILE.data2 = FILE.data %>%
      #put a bunch of file-level variables into each row.
      mutate(FLANK_A_ISFORWARD = FILE.FLANK_A_ISFORWARD) %>%
      mutate(FlankAUltEnd = FILE.FlankAUltEnd) %>% 
      mutate(FOCUS_LOCUS = FILE.FOCUS_LOCUS) %>% 
      mutate(PRIMER_SEQ = FILE.Primer_seq) %>%
      mutate(PRIMER_TO_DSB = FILE.PRIMER_TO_DSB) %>%
      mutate(FLANK_A_REF = FILE.FLANK_A_REF)  %>%
      mutate(TOTAL_REF = FILE.TOTAL_REF) %>%
      mutate(Primer_pos = FILE.Primer_pos) %>%
      mutate(DSB_AREA_SEQ = FILE.DSB_AREA_SEQ)%>%
      mutate(DSB_AREA_SEQ_RC = FILE.DSB_AREA_SEQ_RC)%>%
      mutate(TOTAL_REF_START = FILE.TOTAL_REF_START)%>%
      mutate(TOTAL_REF_STOP = FILE.TOTAL_REF_STOP)%>%
      mutate(Primer_on_TDNA = FILE.Primer_on_TDNA)%>%
      mutate(Primer_OK = TRUE)%>%
      rowwise() %>%
      ungroup()
  }
  
  if (GLOBAL.DEBUG==TRUE){
  FILE.data2_2  = FILE.data2 %>%
    filter(QNAME == GLOBAL.TESTNAME)
  }else{
    FILE.data2_2  = FILE.data2
  }
    
    FILE.data3  = FILE.data2_2 %>%
    #get the length of the ref seq
    mutate(TOTAL_REF_LEN = nchar(TOTAL_REF)) %>%
    #Count number of Ns and remove any reads with Ns
    mutate(NrN = str_count(SEQ_1, pattern = "N"),
           SEQ_1_LEN = nchar(SEQ_1)) %>%
    filter(NrN < 1) %>%
    
    #remove reads that do not have perfectly mapped ends
    filter(A_MAPQ >= GLOBAL.MINMAPQUALA,
           B_MAPQ >= GLOBAL.MINMAPQUALB)%>%
    
    #filter away reads that are too short
    filter(!(SEQ_1_LEN < GLOBAL.MINLEN)) %>%
    
    #calculate the average base quality
    rowwise() %>%
    mutate(AvgBaseQual_1 = mean(utf8ToInt(QUAL_1)-33)) %>%
    mutate(AvgBaseQual_2 = mean(utf8ToInt(QUAL_2)-33)) %>%
    ungroup()

  #here remove all reads where flank A and flank B are on the same chromosome
  #for transguide these are typically removed because we are interested in T-DNA-genome junctions
  if (GLOBAL.REMOVENONTRANS==TRUE){
    FILE.data4 = FILE.data3 %>%
      filter(FLANK_B_CHROM != FILE.FOCUS_CONTIG)
  }else{
    FILE.data4 = FILE.data3
  }
  
  
  FILE.data5 = FILE.data4 %>%
    
    #turn string variables into logical variables
    mutate(FLANK_B_ISFORWARD = if_else(FLANK_B_ORIENT=="FW",
                                       TRUE,
                                       FALSE)) %>%
    mutate(MATE_B_ISFORWARD = if_else(MATE_B_ORIENT=="FW",
                                       TRUE,
                                       FALSE)) %>%
    
    
    #calculate the end position of the B flank in the mate. This is the end the furthest away from the junction.
    mutate(MATE_B_END_POS = case_when(MATE_B_ISFORWARD==TRUE ~ as.integer(MATE_B_POS),
                                      MATE_B_ISFORWARD==FALSE ~ as.integer(MATE_B_POS + 29),
                                      TRUE ~ GLOBAL.ERROR_NUMBER))  %>%
    
    #test whether B flanks from read and mate are mapped to same chromosome, and agree on orientation.
    #note that the two reads are in two different orientations, so the B orientation should not agree
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_ISFORWARD != MATE_B_ISFORWARD,
                                              TRUE,
                                              FALSE)) 
  
  
  
  FILE.data6 = FILE.data5 %>%
    
    #count reads, taking the highest quals
    #acts as basic dupfilter
    #and also summarizes anchors, but puts them in a list
    group_by(
      A_CHROM,
      A_POS,
      FLANK_B_CHROM,
      B_POS,
      FLANK_B_ISFORWARD,
      SEQ_1,
      FILE_NAME,
      PRIMER_SEQ,
      TRIM_LEN,
      SEQ_1_LEN,
      MATE_FLANK_B_CHROM_AGREE,
      MATE_FLANK_B_CHROM,
      FLANK_A_ISFORWARD,
      FlankAUltEnd,
      FOCUS_LOCUS,
      FLANK_A_REF,
      Primer_pos,
      DSB_AREA_SEQ,
      DSB_AREA_SEQ_RC,
      PRIMER_TO_DSB,
      TOTAL_REF,
      TOTAL_REF_START,
      TOTAL_REF_STOP,
      TOTAL_REF_LEN,
      Primer_on_TDNA) %>%
    arrange(desc(AvgBaseQual_1), desc(AvgBaseQual_2)) %>%
    summarize(QNAME_first = dplyr::first(QNAME),
              SEQ_2_first = dplyr::first(SEQ_2),
              SEQ_2_list = toString(SEQ_2),
              MATE_B_END_POS_list = toString(MATE_B_END_POS),
              .groups="drop"
              )%>%

    #then remove reads that don't begin with the primer
    mutate(PRIMER_SEQ_LEN = nchar(PRIMER_SEQ)) %>%
    mutate(SEQ_1_start = substr(SEQ_1, 1, PRIMER_SEQ_LEN)) %>%
    filter(SEQ_1_start == PRIMER_SEQ) 
  
  #check if any reads have survived
  if (nrow(FILE.data6)==0){
    funlog(paste0("No reads surviving for sample ", FILE.Sample))
    #show progress
    GLOBAL.PercentageDone = GLOBAL.PercentageDone + ((FILE.CurrentFileSize/GLOBAL.TotalFileSize)*100)
    funlog(paste0("CISTRANSGUIDE analysis ", round(GLOBAL.PercentageDone, digits=3), "% complete"))
    next
  }
  
  function_time("Step 2 took ")
  
  ###############################################################################
  #Process data: step 3
  ###############################################################################
  
  FILE.data7 = FILE.data6 %>%
    
    #check for expected position and orientation base on primer seqs
    rowwise()%>%
    mutate(FLANK_A_START_POS =  Primer_pos)%>%
    ungroup() %>%
    
    #calculate the end position of the B flank, in order to get the ref seq. This is the end the furthest away from the junction.
    mutate(FLANK_B_END_POS = case_when(FLANK_B_ISFORWARD==FALSE ~ as.integer(B_POS),
                                       FLANK_B_ISFORWARD==TRUE ~ as.integer(B_POS + 29),
                                       TRUE ~ GLOBAL.ERROR_NUMBER))
    

  function_time("Step 3 took ")
  
  ###############################################################################
  #Process data: step 4
  ###############################################################################
  
  FILE.data8 = FILE.data7 %>%

    #Find how much SEQ_1 matches with FLANK_A_REF. Allow 1 bp mismatch somewhere, if the alignment after that continues for at least another 10 bp.
    rowwise() %>%
    mutate(FLANK_A_MATCH = matcher_skipper(FLANK_A_REF, SEQ_1)) %>%
    ungroup() %>%
    mutate(FLANK_A_LEN = nchar(FLANK_A_MATCH)) %>%
    mutate(FLANK_A_END_POS = case_when(FLANK_A_ISFORWARD == TRUE ~ as.integer(FLANK_A_START_POS + (FLANK_A_LEN -1)),
                                       FLANK_A_ISFORWARD == FALSE ~ as.integer(FLANK_A_START_POS - (FLANK_A_LEN -1)),
                                       TRUE ~ GLOBAL.ERROR_NUMBER)) %>%
    
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    mutate(SEQ_1_WO_A_LEN = nchar(SEQ_1_WO_A))%>%
    
    #if the complete read is flank A (or if there are a few bases at the end that dont match) then make the type WT
    mutate(Type=case_when(FLANK_A_MATCH == SEQ_1 ~ "WT", #if the complete read is flank A
                           FLANK_A_MATCH != SEQ_1 & SEQ_1_WO_A_LEN < 2 ~ "SNV", #or if there are a couple of mismatches
                            TRUE~ "OTHER"))
  
  
   if (FILE.DSB_CONTIG == FILE.FOCUS_CONTIG){
    

     #############
    
    #Then with exception of the reads thus far called "WT", determine whether the area surrounding the DSB is intact
    #if this is the case, then the read is probably a WT read with seq errors.
    #check in both seq1 and seq2
    #and check allowing for 1 mismatch, though output this as "SNV" because it may be a real 1bp substitution at the DSB.
  FILE.data9 = FILE.data8 %>%  
  rowwise() %>%
    mutate(DSB_AREA_CHECK = if_else(Type != "WT",
                                    list(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1)),
                                    list("")))%>%
    mutate(DSB_AREA_COUNT = if (Type!="WT"){
      length(DSB_AREA_CHECK@ranges)}else{
        0})%>%
    mutate(DSB_AREA2_CHECK = if_else(Type != "WT" & DSB_AREA_COUNT==0,
                                     list(matchPattern(DNAString(DSB_AREA_SEQ_RC), DNAString(SEQ_2_first), max.mismatch = 1)),
                                     list("")))%>%
    mutate(DSB_AREA2_COUNT = if (Type!="WT"& DSB_AREA_COUNT==0){
      length(DSB_AREA2_CHECK@ranges)}else{
        0})%>%
    ungroup()%>%
    mutate(DSB_AREA_HIT = "",
           DSB_AREA2_HIT = "",
           DSB_AREA_INTACT = FALSE,
           DSB_AREA_1MM = FALSE)
  
  for (j in row.names(FILE.data9)){
    j_int=as.integer(j)
    if (FILE.data9$DSB_AREA_COUNT[[j_int]]==1){
      if ((FILE.data9$SEQ_1_LEN[[j_int]] >= (FILE.data9$DSB_AREA_CHECK[[j_int]]@ranges@start[1]+FILE.data9$DSB_AREA_CHECK[[j_int]]@ranges@width[1]-1)) & (FILE.data9$DSB_AREA_CHECK[[j_int]]@ranges@start[1]>0)){
        FILE.data9$DSB_AREA_HIT[[j_int]] = as.character(FILE.data9$DSB_AREA_CHECK[[j_int]])
      }else{
        FILE.data9$DSB_AREA_HIT[[j_int]] = ""
      }
    }else{
      FILE.data9$DSB_AREA_HIT[[j_int]] = ""
    }
    if (FILE.data9$DSB_AREA2_COUNT[[j_int]]==1){
      if ((nchar(FILE.data9$SEQ_2_first[[j_int]]) >= (FILE.data9$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]+FILE.data9$DSB_AREA2_CHECK[[j_int]]@ranges@width[1]-1)) & (FILE.data9$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]>0)){
        FILE.data9$DSB_AREA2_HIT[[j_int]] = as.character(FILE.data9$DSB_AREA2_CHECK[[j_int]])
      }else{
        FILE.data9$DSB_AREA2_HIT[[j_int]] = ""
      }
    }else{
      FILE.data9$DSB_AREA2_HIT[[j_int]] = ""
    }
  }
  
  FILE.data10 =   FILE.data9 %>%
    mutate(DSB_AREA_INTACT_SEQ1 = if_else(DSB_AREA_HIT == DSB_AREA_SEQ,
                                     "TRUE",
                                     "FALSE"))%>%
    mutate(DSB_AREA_INTACT_SEQ2 = if_else(DSB_AREA_HIT == DSB_AREA_SEQ | DSB_AREA2_HIT == DSB_AREA_SEQ_RC ,
                                           "TRUE",
                                           "FALSE"))%>%
    mutate(DSB_AREA_1MM_SEQ1 = if_else(
      DSB_AREA_INTACT_SEQ1 == FALSE & (DSB_AREA_COUNT>0),
      "TRUE",
      "FALSE")) %>%
    mutate(DSB_AREA_1MM_SEQ2 = if_else(
      DSB_AREA_INTACT_SEQ2 == FALSE & (DSB_AREA2_COUNT>0),
      "TRUE",
      "FALSE")) %>%
    mutate(DSB_HIT_MULTI_SEQ1 = if_else(DSB_AREA_COUNT>1,
                                   "TRUE",
                                   "FALSE")) %>%
    mutate(DSB_HIT_MULTI_SEQ2 = if_else(DSB_AREA2_COUNT>1,
                                        "TRUE",
                                        "FALSE")) %>%
    #set the type again
    mutate(Type = case_when(Type=="WT" | DSB_AREA_INTACT_SEQ1 == TRUE ~ "WT",
                            DSB_AREA_INTACT_SEQ1 == FALSE & DSB_AREA_INTACT_SEQ2 == FALSE & (DSB_AREA_1MM_SEQ1 == TRUE | DSB_AREA_INTACT_SEQ2 == TRUE | DSB_AREA_1MM_SEQ2 == TRUE) ~ "SNV",
                            TRUE ~ Type))
  
   }else{
     FILE.data10 = FILE.data8 %>%
       mutate(
  #the following variables are set to FALSE. This does not mean that there is a DSB area that does not match wt sequence, but rather that it has not been determined.
       DSB_AREA_INTACT_SEQ1 = FALSE,
       DSB_AREA_INTACT_SEQ2 = FALSE,
       DSB_AREA_1MM_SEQ1 = FALSE,
       DSB_AREA_1MM_SEQ2 = FALSE,
       DSB_HIT_MULTI_SEQ1 = FALSE,
       DSB_HIT_MULTI_SEQ2 = FALSE)
     }
  
  
  
    
  FILE.data11 = FILE.data10 %>%
    
    #check whether FLANK_A ends within the T-DNA. (IF LB or RB transguide reaction. if not, limit to the T-DNA)
    mutate(FLANK_A_ENDS_ON_TDNA = case_when(
      
      FILE.FOCUS_CONTIG == FILE.PLASMID & FOCUS_LOCUS=="LB" & FILE.TDNA_IS_LBRB == TRUE & FLANK_A_END_POS >= FILE.TDNA_LB_END ~ TRUE,
      FILE.FOCUS_CONTIG == FILE.PLASMID & FOCUS_LOCUS=="LB" & FILE.TDNA_IS_LBRB == FALSE & FLANK_A_END_POS <= FILE.TDNA_LB_END ~ TRUE,
      FILE.FOCUS_CONTIG == FILE.PLASMID & FOCUS_LOCUS=="RB" & FILE.TDNA_IS_LBRB == TRUE & FLANK_A_END_POS <= FILE.TDNA_RB_END ~ TRUE,
      FILE.FOCUS_CONTIG == FILE.PLASMID & FOCUS_LOCUS=="RB" & FILE.TDNA_IS_LBRB == FALSE & FLANK_A_END_POS >= FILE.TDNA_RB_END ~ TRUE,
      
      TRUE ~ FALSE))%>%
    
    #fix end pos for the WT/SNV case
    mutate(FLANK_A_LEN = if_else(Type=="WT"| Type=="SNV" | (FLANK_A_ENDS_ON_TDNA==FALSE & FILE.FOCUS_CONTIG == FILE.PLASMID),
                                  PRIMER_TO_DSB,
                                  FLANK_A_LEN))%>%
    mutate(FLANK_A_END_POS = if_else(Type=="WT"| Type=="SNV"| (FLANK_A_ENDS_ON_TDNA==FALSE & FILE.FOCUS_CONTIG == FILE.PLASMID),
                                     FlankAUltEnd,
                                     FLANK_A_END_POS))%>%
    #also determine again what the part of the read is without FLANK_A
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    mutate(SEQ_1_WO_A_LEN = nchar(SEQ_1_WO_A))%>%
    
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        FLANK_A_LEN != GLOBAL.ERROR_NUMBER & FLANK_A_LEN != GLOBAL.NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
        TRUE ~ GLOBAL.ERROR_NUMBER))


  function_time("Step 4 took ")
  
  ###############################################################################
  #Process data: step 5
  ###############################################################################
  
  FILE.data12 = FILE.data11 %>%
    mutate(FLANK_B_CLOSE_BY = case_when(FLANK_A_ISFORWARD == TRUE & TOTAL_REF_STOP > FLANK_B_END_POS &  FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_B_CHROM == FILE.FOCUS_CONTIG ~ TRUE,
                                        FLANK_A_ISFORWARD == FALSE & TOTAL_REF_START < FLANK_B_END_POS &  FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_B_CHROM == FILE.FOCUS_CONTIG ~ TRUE,
                                        TRUE ~ FALSE))%>%

    #FLANK_B_REF. This ref includes homology.
    rowwise() %>%
    mutate(FLANK_B_REF =
             if (Type == "WT" | Type=="SNV"){ #if we already know it's a wt read, don't get the ref
               ""
             }else{
             if (FLANK_B_CLOSE_BY==TRUE){
               if (FLANK_B_ISFORWARD == TRUE){
                 substr(TOTAL_REF, start = 1, stop = (TOTAL_REF_LEN - (TOTAL_REF_STOP-FLANK_B_END_POS)))
               }else{
                 substr(TOTAL_REF, start = 1, stop = (TOTAL_REF_LEN - (FLANK_B_END_POS-TOTAL_REF_START)))
               } 
             }else{ #translocations
               if (FLANK_B_ISFORWARD == TRUE){
                 substr(as.character(eval(parse(text = paste0("FILE.genomeseq$`", FLANK_B_CHROM, "`")))), FLANK_B_END_POS-(SEQ_1_LEN-1), FLANK_B_END_POS)
               }else{
                 as.character(reverseComplement(DNAString(substr(as.character(eval(parse(text = paste0("FILE.genomeseq$`", FLANK_B_CHROM, "`")))), FLANK_B_END_POS, FLANK_B_END_POS+(SEQ_1_LEN-1)))))
               }}})%>%
    #find the full flank b match, while skipping over seq errors. But causes a problem on focus contig when del=1 en ins=1
    #skip if the read is wt
    mutate(FLANK_B_MATCH = if_else(Type != "WT" & Type != "SNV",
                                   stri_reverse(matcher_skipper(stri_reverse(FLANK_B_REF), stri_reverse(SEQ_1))),
                                    ""))%>%     
    ungroup() %>%
    mutate(FLANK_B_MATCH_LEN = nchar(FLANK_B_MATCH)) %>%
  


    #flank b start position including MH
    mutate(FLANK_B_START_POS_MH = case_when(
      
      FLANK_B_ISFORWARD == TRUE & Type!= "WT" & Type!= "SNV" ~ as.integer(FLANK_B_END_POS-(FLANK_B_MATCH_LEN-1)),
      FLANK_B_ISFORWARD == FALSE & Type!= "WT" & Type!= "SNV" ~ as.integer(FLANK_B_END_POS+(FLANK_B_MATCH_LEN-1)),
      
      #if WT reads, FLANK B should start the base after the final allowed base of flank A
      (Type=="WT" | Type=="SNV") & FLANK_A_ISFORWARD == TRUE ~ as.integer(FLANK_A_END_POS+1),
      (Type=="WT" | Type=="SNV") & FLANK_A_ISFORWARD == FALSE ~ as.integer(FLANK_A_END_POS-1),

                                            TRUE ~ GLOBAL.ERROR_NUMBER)) %>%

    
    #get the FLANK B len including MH, unless the read is WT
    mutate(FLANK_B_LEN_MH = case_when(Type!="WT" & Type!="SNV" & FLANK_B_ISFORWARD == TRUE ~ FLANK_B_END_POS-(FLANK_B_START_POS_MH-1),
                                      Type!="WT" & Type!="SNV" & FLANK_B_ISFORWARD == FALSE ~ FLANK_B_START_POS_MH-(FLANK_B_END_POS-1),
                                      TRUE ~ SEQ_1_WO_A_LEN))%>%
             
    rowwise() %>%
    #Extract the MH sequence (if there is a TD, this is wrong)
    mutate(MH_TD = if_else(SEQ_1_WO_A_LEN < FLANK_B_LEN_MH,
                        substr(SEQ_1, start = (1 + SEQ_1_LEN - FLANK_B_LEN_MH), stop = FLANK_A_LEN),
                        "")) %>%
    ungroup()%>%
    mutate(MH_TD_LEN = nchar(MH_TD))%>%
    
    #if there is MH, but A+B-MH does not equal seq1, then there is TD
    mutate(hasTandemDuplication = case_when(
      FLANK_B_CHROM==FILE.FOCUS_CONTIG & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD == TRUE & FLANK_A_END_POS >= (FLANK_B_START_POS_MH+MH_TD_LEN) ~ TRUE,
      FLANK_B_CHROM==FILE.FOCUS_CONTIG & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD == FALSE & FLANK_A_END_POS <= (FLANK_B_START_POS_MH-MH_TD_LEN) ~ TRUE,
      
      
      TRUE ~ FALSE))%>%
             

    #adjust the MH
    mutate(MH = if_else(hasTandemDuplication == TRUE,
                        "",
                        MH_TD))%>%
    
    #flank B start pos excluding MH or TD, for del calculation
    mutate(FLANK_B_START_POS = case_when(Type!="WT" & Type!="SNV" & hasTandemDuplication == TRUE & FLANK_B_ISFORWARD==TRUE ~ FLANK_A_END_POS +1 ,
                                         Type!="WT" & Type!="SNV" & hasTandemDuplication == TRUE & FLANK_B_ISFORWARD==FALSE ~ FLANK_A_END_POS - 1,
                                         hasTandemDuplication == FALSE & nchar(MH)>0 & FLANK_B_ISFORWARD == TRUE ~ FLANK_B_START_POS_MH + nchar(MH),
                                         hasTandemDuplication == FALSE & nchar(MH)>0 & FLANK_B_ISFORWARD == FALSE ~ FLANK_B_START_POS_MH - nchar(MH),
                                         TRUE ~ FLANK_B_START_POS_MH))%>%

    #calculate length of flank B minus the MH and TD
    mutate(FLANK_B_LEN = case_when(Type!="WT" & Type!="SNV" & FLANK_B_ISFORWARD== TRUE ~ FLANK_B_END_POS-(FLANK_B_START_POS-1),
                                   Type!="WT" & Type!="SNV" & FLANK_B_ISFORWARD== FALSE ~ FLANK_B_START_POS-(FLANK_B_END_POS-1),
                                   TRUE ~ SEQ_1_WO_A_LEN))
            

  function_time("Step 5 took ")
  
  ###############################################################################
  #Process data: step 6
  ###############################################################################
  
  FILE.data13 = FILE.data12 %>%
    
    #determine the filler sequence 
    mutate(
      FILLER = if_else(
        FLANK_B_LEN < SEQ_1_WO_A_LEN,
        substr(SEQ_1_WO_A, start = 1, stop = SEQ_1_WO_A_LEN - FLANK_B_LEN),	 
        "" #no filler
      )) %>%
    mutate(insSize = nchar(FILLER)) %>%
    

    #get the tandem duplication out of the filler
    #match with total ref, because tandem duplication can extend beyond start of read
    rowwise()%>%
    mutate(tandemDuplicationLength = if(hasTandemDuplication == TRUE){
      if (FLANK_A_ISFORWARD==TRUE){
        lcprefix(reverse(substr(TOTAL_REF, start=1, stop= 1+FLANK_A_END_POS-TOTAL_REF_START)), reverse(FILLER))
      }else{
        lcprefix(reverse(substr(TOTAL_REF, start=1, stop= 1+TOTAL_REF_STOP-FLANK_A_END_POS )), reverse(FILLER))
      }
    }else{
      0
    })%>%
    
    mutate(TANDEM_DUPLICATION = if_else(hasTandemDuplication == TRUE,
                                        substr(FILLER, start=(1+insSize-tandemDuplicationLength), stop=insSize),
                                        ""))%>%
    ungroup()%>%

  
    #determine how much from flank B has been deleted
    #only report deletion when ends are known
    mutate(FLANK_B_ON_TDNA = case_when(FLANK_B_CHROM == FILE.PLASMID & FILE.TDNA_IS_LBRB == TRUE & FLANK_B_START_POS >= FILE.TDNA_LB_END & FLANK_B_START_POS <= FILE.TDNA_RB_END ~ TRUE,
                                       FLANK_B_CHROM == FILE.PLASMID & FILE.TDNA_IS_LBRB == FALSE & FLANK_B_START_POS <= FILE.TDNA_LB_END & FLANK_B_START_POS >= FILE.TDNA_RB_END ~ TRUE,
                                       FLANK_B_CHROM == FILE.PLASMID_ALT & FILE.TDNA_ALT_IS_LBRB == TRUE & FLANK_B_START_POS >= FILE.TDNA_ALT_LB_END & FLANK_B_START_POS <= FILE.TDNA_ALT_RB_END ~ TRUE,
                                       FLANK_B_CHROM == FILE.PLASMID_ALT & FILE.TDNA_ALT_IS_LBRB == FALSE & FLANK_B_START_POS <= FILE.TDNA_ALT_LB_END & FLANK_B_START_POS >= FILE.TDNA_ALT_RB_END ~ TRUE,
                                       TRUE ~ FALSE
                                       ))%>%
    mutate(FLANK_B_DEL = case_when(FLANK_B_CHROM == FILE.DSB_CONTIG & FLANK_B_ISFORWARD == TRUE   ~ as.integer(FLANK_B_START_POS - (FILE.DSB_FW_END+1)),
                                   FLANK_B_CHROM == FILE.DSB_CONTIG & FLANK_B_ISFORWARD == FALSE  ~ as.integer(FILE.DSB_FW_END - FLANK_B_START_POS),
                                   
                                   FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - FILE.TDNA_LB_END,
                                   FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - FILE.TDNA_RB_END,
                                   FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FILE.TDNA_RB_END - FLANK_B_START_POS,
                                   FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FILE.TDNA_LB_END - FLANK_B_START_POS,
                                   
                                   FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - FILE.TDNA_ALT_LB_END,
                                   FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - FILE.TDNA_ALT_RB_END,
                                   FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FILE.TDNA_ALT_RB_END - FLANK_B_START_POS,
                                   FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FILE.TDNA_ALT_LB_END - FLANK_B_START_POS,
                                   
                                   TRUE ~ GLOBAL.ERROR_NUMBER))%>%
    
    #also report what side of the T-DNA flank B is
    mutate(FLANK_B_TDNA_SIDE = case_when(FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == FILE.PLASMID & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         
                                         FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & FILE.TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == FILE.PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & FILE.TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         
                                         TRUE ~ NA))%>%
    
      
    #calculate total deletion length
    #only for deletions around the DSB with no translocations
    mutate(delSize = case_when(Type=="WT" | Type=="SNV" ~ 0,
                               Type!= "WT" & Type!="SNV" & FLANK_B_CHROM == FILE.FOCUS_CONTIG & FLANK_B_ISFORWARD == FLANK_A_ISFORWARD & FLANK_B_ISFORWARD == TRUE & FLANK_A_END_POS < FLANK_B_START_POS & FLANK_B_DEL != GLOBAL.ERROR_NUMBER & FLANK_B_CHROM != FILE.PLASMID ~ as.integer(FLANK_A_DEL + FLANK_B_DEL),
                               Type!= "WT" & Type!="SNV" & FLANK_B_CHROM == FILE.FOCUS_CONTIG & FLANK_B_ISFORWARD == FLANK_A_ISFORWARD & FLANK_B_ISFORWARD == FALSE & FLANK_A_END_POS > FLANK_B_START_POS & FLANK_B_DEL != GLOBAL.ERROR_NUMBER & FLANK_B_CHROM != FILE.PLASMID ~ as.integer(FLANK_A_DEL + FLANK_B_DEL),
                               TRUE ~ GLOBAL.ERROR_NUMBER))%>%

    #correct MH if delSize == 0
    mutate(MH = if_else(delSize <= 0,
                        "",
                        MH)) %>%
    mutate(homologyLength = if_else(insSize == 0,
                                    nchar(MH),
                                    as.integer(-1))) %>% 
    

    #remove reads with FLANK_B 's shorter than the minimum.
    filter(FLANK_B_LEN_MH >= GLOBAL.FLANK_B_LEN_MIN)
  
  function_time("Step 6 took ")
  
  ###############################################################################
  #Process data: step 7
  ###############################################################################
  
  FILE.data14 = FILE.data13 %>%
  
    #for SIQPlotter delRelativeStart is FLANK_A_DEL (but negative) and delRelativeEnd is FLANK_B_DEL.
    #delRelativeEndTD should include the tandem duplication.

    mutate(
      delRelativeStart = as.integer(FLANK_A_DEL*-1),
      delRelativeStartTD = as.integer(FLANK_A_DEL*-1),
      delRelativeEnd = if_else(abs(as.integer(FLANK_B_DEL))<=GLOBAL.FLANKBEYONDDSB | FLANK_B_CHROM == FILE.PLASMID | FLANK_B_CHROM == FILE.PLASMID_ALT,
                               as.integer(FLANK_B_DEL),
                               GLOBAL.ERROR_NUMBER),
      delRelativeEndTD = if_else( (abs(as.integer(FLANK_B_DEL))<=GLOBAL.FLANKBEYONDDSB  | FLANK_B_CHROM == FILE.PLASMID | FLANK_B_CHROM == FILE.PLASMID_ALT) & FLANK_B_DEL!= GLOBAL.ERROR_NUMBER,
                               as.integer(FLANK_B_DEL-tandemDuplicationLength),
                               GLOBAL.ERROR_NUMBER)) %>%

  
    
    #then check for fake fillers that are caused by seq errors.
    #if delsize == inssize this may be the case 
    
    mutate(POT_FAKE_INS = case_when(insSize == delSize & delSize > 9 ~ as.character(substr(TOTAL_REF, start = FLANK_A_LEN, stop = FLANK_A_LEN + insSize -1)),
                                    TRUE ~ "")) %>%
    rowwise() %>%
    mutate(FAKE_DELIN_CHECK = if (insSize == delSize & delSize > 0) {
      if ((countPattern(
        DNAString(FILLER),
        DNAString(POT_FAKE_INS),
        max.mismatch = (abs(insSize * 0.1))
      )) > 0) {
        TRUE
      }
      else{
        FALSE
      }
    }
    else{
      FALSE
    })%>%
    
    #calculate the minimum length of the read to capture the entire outcome
    #then trim the reads to that minimum.
    #any that have less than this minimum get removed
    mutate(read_minimum_length = FLANK_A_LEN + insSize + 30) %>%
    filter(SEQ_1_LEN >= read_minimum_length)%>%
    mutate(SEQ_1_trimmed = substr(SEQ_1, 1, read_minimum_length))
  
  #here grouping will occur to determine the number of anchors.
  #for TRANSGUIDE a consensus outcome will be determined
         if (GLOBAL.GROUPSAMEPOS == FALSE){
           if (GLOBAL.UNGROUPMATES == FALSE){
      
           FILE.data15 =FILE.data14 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
      group_by(
        FILE_NAME,
        PRIMER_SEQ,
        delRelativeStart,
        delRelativeStartTD,
        delRelativeEnd,
        delRelativeEndTD,
        FLANK_A_LEN,
        FLANK_A_END_POS,
        FLANK_B_CHROM,
        FLANK_B_START_POS,
        MATE_FLANK_B_CHROM_AGREE,
        MATE_FLANK_B_CHROM,
        FLANK_B_ISFORWARD,
        FILLER,
        MH,
        insSize,
        delSize,
        homologyLength,
        FAKE_DELIN_CHECK,
        DSB_AREA_INTACT_SEQ1,
        DSB_AREA_INTACT_SEQ2,
        DSB_AREA_1MM_SEQ1,
        DSB_AREA_1MM_SEQ2,
        DSB_HIT_MULTI_SEQ1,
        DSB_HIT_MULTI_SEQ2,
        Type,
        TRIM_LEN,
        TANDEM_DUPLICATION,
        tandemDuplicationLength,
        FLANK_B_TDNA_SIDE,
        FLANK_A_ISFORWARD,
        FOCUS_LOCUS,
        Primer_on_TDNA,
        FlankAUltEnd,
        FLANK_A_START_POS
      )%>%
      summarize(
        ReadCount = n(),
        SEQ_1_con = names(which.max(table(SEQ_1_trimmed))),
        Name = dplyr::first(QNAME_first),
        Count_consensus = max(table(SEQ_1)),
        AnchorCount = n_distinct(MATE_B_END_POS_list),
        SEQ_2_con = names(which.max(table(SEQ_2_first))),
        #select the extreme mate positions
        MATE_B_END_POS_max = max(as.integer(MATE_B_END_POS_list)),
        MATE_B_END_POS_min = min(as.integer(MATE_B_END_POS_list)),
        SEQ_1_LEN_max = max(SEQ_1_LEN),
        .groups="drop"
      )%>%
          mutate(Consensus_freq = 1)
           #if it is important that the mates are kept separate
           }else{
             FILE.data15 =FILE.data14 %>%
               ungroup()%>%
               separate_longer_delim(cols=c("MATE_B_END_POS_list","SEQ_2_list"), delim = ",") %>%
               group_by(
                 FILE_NAME,
                 PRIMER_SEQ,
                 delRelativeStart,
                 delRelativeStartTD,
                 delRelativeEnd,
                 delRelativeEndTD,
                 FLANK_A_LEN,
                 FLANK_A_END_POS,
                 FLANK_B_CHROM,
                 FLANK_B_START_POS,
                 MATE_FLANK_B_CHROM_AGREE,
                 MATE_FLANK_B_CHROM,
                 FLANK_B_ISFORWARD,
                 FILLER,
                 MH,
                 insSize,
                 delSize,
                 homologyLength,
                 FAKE_DELIN_CHECK,
                 DSB_AREA_INTACT_SEQ1,
                 DSB_AREA_INTACT_SEQ2,
                 DSB_AREA_1MM_SEQ1,
                 DSB_AREA_1MM_SEQ2,
                 DSB_HIT_MULTI_SEQ1,
                 DSB_HIT_MULTI_SEQ2,
                 Type,
                 TRIM_LEN,
                 TANDEM_DUPLICATION,
                 tandemDuplicationLength,
                 FLANK_B_TDNA_SIDE,
                 FLANK_A_ISFORWARD,
                 FOCUS_LOCUS,
                 Primer_on_TDNA,
                 FlankAUltEnd,
                 FLANK_A_START_POS,
                 MATE_B_END_POS_list,
                 SEQ_2_list
               )%>%
               summarize(
                 ReadCount = n(),
                 SEQ_1_con = names(which.max(table(SEQ_1_trimmed))),
                 Name = dplyr::first(QNAME_first),
                 Count_consensus = max(table(SEQ_1)),
                 SEQ_1_LEN_max = max(SEQ_1_LEN),
                 .groups="drop"
               )%>%
               mutate(Consensus_freq = 1,
                      MATE_B_END_POS_max = as.integer(MATE_B_END_POS_list),
                      MATE_B_END_POS_min = as.integer(MATE_B_END_POS_list),
                      AnchorCount = 1,
                      SEQ_2_con = SEQ_2_list)
           }
           
           
           
      }else{
        
        FILE.data15 =FILE.data14 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
          group_by(
            FILE_NAME,
            PRIMER_SEQ,
            FLANK_B_CHROM,
            FLANK_B_START_POS,
            FLANK_B_ISFORWARD,
            TRIM_LEN,
            MATE_FLANK_B_CHROM_AGREE,
            MATE_FLANK_B_CHROM,
            FLANK_B_TDNA_SIDE,
            FLANK_A_ISFORWARD,
            FOCUS_LOCUS,
            Primer_on_TDNA,
            FlankAUltEnd,
            FLANK_A_START_POS
          )%>%
          summarize(
            ReadCount = n(),
            Name = dplyr::first(QNAME_first),
            Count_consensus = max(table(SEQ_1_trimmed)),
            AnchorCount = n_distinct(MATE_B_END_POS_list),
            SEQ_1_con = names(which.max(table(SEQ_1_trimmed))),
            SEQ_2_con = names(which.max(table(SEQ_2_first))),
            #select the extreme mate positions
            MATE_B_END_POS_max = max(as.integer(MATE_B_END_POS_list)),
            MATE_B_END_POS_min = min(as.integer(MATE_B_END_POS_list)),
            delRelativeStart_con = as.integer(names(which.max(table(delRelativeStart)))),
            delRelativeEnd_con = as.integer(names(which.max(table(delRelativeEnd)))),
            delRelativeStartTD_con = as.integer(names(which.max(table(delRelativeStartTD)))),
            delRelativeEndTD_con = as.integer(names(which.max(table(delRelativeEndTD)))),
            FLANK_A_END_POS_con = as.integer(names(which.max(table(FLANK_A_END_POS)))),
            FILLER_con = names(which.max(table(FILLER))),
            MH_con =names(which.max(table(MH))),
            insSize_con = as.integer(names(which.max(table(insSize)))),
            delSize_con = as.integer(names(which.max(table(delSize)))),
            homologyLength_con = as.integer(names(which.max(table(homologyLength)))),
            FAKE_DELIN_CHECK_con = as.logical(names(which.max(table(FAKE_DELIN_CHECK)))),
            TANDEM_DUPLICATION_con = names(which.max(table(TANDEM_DUPLICATION))),
            tandemDuplicationLength_con = as.integer(names(which.max(table(tandemDuplicationLength)))),
            DSB_AREA_INTACT_SEQ1_con = as.logical(names(which.max(table(DSB_AREA_INTACT_SEQ1)))),
            DSB_AREA_INTACT_SEQ2_con = as.logical(names(which.max(table(DSB_AREA_INTACT_SEQ2)))),
            DSB_AREA_1MM_SEQ1_con = as.logical(names(which.max(table(DSB_AREA_1MM_SEQ1)))),
            DSB_AREA_1MM_SEQ2_con = as.logical(names(which.max(table(DSB_AREA_1MM_SEQ2)))),
            DSB_HIT_MULTI_SEQ1_con = as.logical(names(which.max(table(DSB_HIT_MULTI_SEQ1)))),
            DSB_HIT_MULTI_SEQ2_con = as.logical(names(which.max(table(DSB_HIT_MULTI_SEQ2)))),
            Type = names(which.max(table(Type))),
            SEQ_1_LEN_max = max(SEQ_1_LEN),
            .groups="drop"
          )%>%
          mutate(Consensus_freq = Count_consensus/ReadCount)%>%
          #rename columns to that code below is compatible with TRANS and CISGUIDE
          rename(FILLER = FILLER_con,
                 MH = MH_con,
                 FAKE_DELIN_CHECK = FAKE_DELIN_CHECK_con,
                 DSB_AREA_INTACT_SEQ1 = DSB_AREA_INTACT_SEQ1_con,
                 DSB_AREA_INTACT_SEQ2 = DSB_AREA_INTACT_SEQ2_con,
                 DSB_AREA_1MM_SEQ1 = DSB_AREA_1MM_SEQ1_con,
                 DSB_AREA_1MM_SEQ2 = DSB_AREA_1MM_SEQ2_con,
                 DSB_HIT_MULTI_SEQ1 = DSB_HIT_MULTI_SEQ1_con,
                 DSB_HIT_MULTI_SEQ2 = DSB_HIT_MULTI_SEQ2_con,
                 delRelativeStart = delRelativeStart_con,
                 delRelativeEnd = delRelativeEnd_con,
                 delRelativeStartTD = delRelativeStartTD_con,
                 delRelativeEndTD = delRelativeEndTD_con,
                 insSize = insSize_con,
                 delSize = delSize_con,
                 homologyLength = homologyLength_con,
                 FLANK_A_END_POS = FLANK_A_END_POS_con,
                 TANDEM_DUPLICATION = TANDEM_DUPLICATION_con,
                 tandemDuplicationLength = tandemDuplicationLength_con)
      }
        
  FILE.data16 =FILE.data15 %>%
      
    #determine whether a translocation has occurred
    mutate(Translocation = case_when(Type=="WT" | Type=="SNV" ~ FALSE,
                                     FILE.FOCUS_CONTIG ==  FLANK_B_CHROM & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD==TRUE & FLANK_A_END_POS < FLANK_B_START_POS ~ FALSE,
                                     FILE.FOCUS_CONTIG ==  FLANK_B_CHROM & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD==FALSE & FLANK_B_START_POS < FLANK_A_END_POS ~ FALSE,
                                     
                                     TRUE ~ TRUE))%>% 
    
     #then calculate the difference between start of flank B (in the read) and the position of of flank B at the end of the mate.
      mutate(ANCHOR_DIST = case_when(Translocation == TRUE & FLANK_B_ISFORWARD == TRUE & MATE_FLANK_B_CHROM_AGREE == TRUE & MATE_B_END_POS_max > FLANK_B_START_POS ~ as.integer(1+MATE_B_END_POS_max - FLANK_B_START_POS),
                                     Translocation == TRUE & FLANK_B_ISFORWARD == FALSE & MATE_FLANK_B_CHROM_AGREE == TRUE & FLANK_B_START_POS > MATE_B_END_POS_min ~ as.integer(1+FLANK_B_START_POS - MATE_B_END_POS_min),
                                     Translocation == FALSE & FLANK_A_ISFORWARD == TRUE & MATE_FLANK_B_CHROM_AGREE == TRUE ~ as.integer(1+MATE_B_END_POS_max - FLANK_A_START_POS),
                                     Translocation == FALSE & FLANK_A_ISFORWARD == FALSE & MATE_FLANK_B_CHROM_AGREE == TRUE ~ as.integer(1+FLANK_A_START_POS - MATE_B_END_POS_min),
                                     TRUE ~ GLOBAL.NF_NUMBER)) %>%
      #adjust the anchor dist for when the anchor is impossibly far away
      mutate(ANCHOR_DIST = if_else(ANCHOR_DIST > GLOBAL.MAXANCHORDIST & MATE_FLANK_B_CHROM!=FILE.PLASMID & MATE_FLANK_B_CHROM!= FILE.PLASMID_ALT,
                                   GLOBAL.NF_NUMBER,
                                   ANCHOR_DIST)) %>%

     
      
      #note whether the deletion of the flank B (in case of translocation) has been determined
      mutate(Translocation_del_resolved = if_else(Translocation == TRUE & delRelativeEnd != GLOBAL.ERROR_NUMBER,
                                                  TRUE,
                                                  FALSE))%>%
      #if the deletion size cannot be determined it becomes 0 for plotting purposes
      mutate(delRelativeEnd = if_else(delRelativeEnd == GLOBAL.ERROR_NUMBER,
                                      0,
                                      delRelativeEnd),
             delRelativeEndTD = if_else(delRelativeEndTD == GLOBAL.ERROR_NUMBER,
                                      0,
                                      delRelativeEndTD))%>%
    
    
              
      #Add a Subject column
      #adjust the type column again.
    mutate(
      Subject = FOCUS_LOCUS,
      Type = case_when(
        Type=="WT"  ~ "WT",
        Type != "WT" & (Type=="SNV" | FAKE_DELIN_CHECK == TRUE)  ~ "SNV",

        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize != 0 & insSize != 0 & tandemDuplicationLength==0 ~ "DELINS",
        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize != 0 & insSize != 0 & tandemDuplicationLength!=0 ~ "TANDEMDUPLICATION_COMPOUND",
        

        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize == 0 & insSize != 0 & tandemDuplicationLength < GLOBAL.TD_SIZE_CUTOFF ~ "INSERTION",
        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize == 0 & insSize != 0 & tandemDuplicationLength >= GLOBAL.TD_SIZE_CUTOFF ~ "TANDEMDUPLICATION",
        
        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & insSize == 0 ~ "DELETION",
        
        TRUE ~ "OTHER") #this "other" is a catch category for bugs. Nothing should have this value.
    ) %>%
      
    
    
    #select the most important columns
    select(
      Name,
      AnchorCount,
      ANCHOR_DIST,
      ReadCount,
      delRelativeStart,
      delRelativeStartTD,
      delRelativeEnd,
      delRelativeEndTD,
      delSize,
      FILLER,
      MH,
      TANDEM_DUPLICATION,
      insSize,
      homologyLength,
      tandemDuplicationLength,
      FLANK_B_CHROM,
      FLANK_B_START_POS,
      FLANK_B_ISFORWARD,
      Type,
      SEQ_1_con,
      SEQ_2_con,
      DSB_AREA_INTACT_SEQ1,
      DSB_AREA_INTACT_SEQ2,
      DSB_AREA_1MM_SEQ1,
      DSB_AREA_1MM_SEQ2,
      DSB_HIT_MULTI_SEQ1,
      DSB_HIT_MULTI_SEQ2,
      Consensus_freq,
      FAKE_DELIN_CHECK,
      MATE_FLANK_B_CHROM_AGREE,
      MATE_FLANK_B_CHROM,
      FLANK_A_END_POS,
      Translocation,
      Translocation_del_resolved,
      FLANK_B_TDNA_SIDE,
      Subject,
      FILE_NAME,
      PRIMER_SEQ,
      TRIM_LEN,
      FLANK_A_ISFORWARD,
      FOCUS_LOCUS,
      Primer_on_TDNA,
      FlankAUltEnd,
      SEQ_1_LEN_max,
      MATE_B_END_POS_max
    ) 

   

  function_time("Step 7 took ")
  
  ###############################################################################
  #Process data: step 8
  ###############################################################################
  
  #calculate the fraction of reads with a certain outcome within a library
  FILE.data17 = FILE.data16 %>% group_by(FILE_NAME) %>% summarize(ReadCountTotal =
                                                                                  sum(ReadCount),
                                                                        .groups="drop")
  function_time("Step 8 took ")
  
  ###############################################################################
  #Process data: step 9
  ###############################################################################
  
  FILE.data18 = left_join(FILE.data16, FILE.data17, by = "FILE_NAME") %>%
    mutate(fraction = ReadCount / ReadCountTotal) %>%
    #add columns for all the input options and software version
    mutate(countReadsTotal = NULL,
           Focus_contig = FILE.FOCUS_CONTIG,
           Genotype = FILE.Genotype,
           DNASample = FILE.DNASample,
           Ecotype = FILE.Ecotype,
           RunID = FILE.RunID,
           Plasmid = FILE.PLASMID,
           Plasmid_alt = FILE.PLASMID_ALT,
           AgroGeno = FILE.AgroGeno,
           Species = FILE.Species,
           Alias = paste0(FILE.Sample, "_", FILE.RunID),
           DSB_FW_END = FILE.DSB_FW_END,
           DSB_CONTIG = FILE.DSB_CONTIG,
           TDNA_LB_END = FILE.TDNA_LB_END,
           TDNA_RB_END = FILE.TDNA_RB_END,
           TDNA_ALT_LB_END = FILE.TDNA_ALT_LB_END,
           TDNA_ALT_RB_END = FILE.TDNA_ALT_RB_END,
           TDNA_IS_LBRB = FILE.TDNA_IS_LBRB,
           TDNA_ALT_IS_LBRB = FILE.TDNA_ALT_IS_LBRB,
           Experiment = FILE.Experiment,
           MinumumReadLength = GLOBAL.MINLEN,
           program_version = GLOBAL.hash,
           TD_SIZE_CUTOFF = GLOBAL.TD_SIZE_CUTOFF,
           FLANK_B_LEN_MIN = GLOBAL.FLANK_B_LEN_MIN,
           FLANKBEYONDDSB = GLOBAL.FLANKBEYONDDSB,
           RemoveNonTranslocation = GLOBAL.REMOVENONTRANS,
           GroupSamePosition = GLOBAL.GROUPSAMEPOS,
           MAXANCHORDIST = GLOBAL.MAXANCHORDIST)
  

  
  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, FILE.data18)
  saveWorkbook(work_book, file = paste0(GLOBAL.output_dir, FILE.Sample, "_", FILE.RunID, "_CISTRANSGUIDE_V2.xlsx"), overwrite = TRUE)
  
  
  function_time("Step 9 took ")
  
  
  #show progress
  GLOBAL.PercentageDone = GLOBAL.PercentageDone + ((FILE.CurrentFileSize/GLOBAL.TotalFileSize)*100)
  funlog(paste0("CISTRANSGUIDE analysis ", round(GLOBAL.PercentageDone, digits=3), "% complete"))
  
}


###############################################################################
#Combine data: step 10
###############################################################################

GLOBAL.sample_list = list.files(path=GLOBAL.output_dir, pattern = "CISTRANSGUIDE_V2.xlsx")
GLOBAL.wb_pre = tibble()

for (i in GLOBAL.sample_list){
  GLOBAL.wb_pre=bind_rows(GLOBAL.wb_pre, read.xlsx(paste0(GLOBAL.output_dir, i)) %>% select_if(function(x) !(all(is.na(x)) | all(x==""))))

}
#in case the following columns are NA, they will not have been imported from the excel. I therefore need to add them. But also allowing the possibility that they are already there.
GLOBAL.missing_columns = GLOBAL.wb_pre %>%
  select(Name) %>%
  mutate(Plasmid_alt = NA,
         DSB_FW_END = NA,
         DSB_OVERHANG = NA,
         DSB_CONTIG = NA,
         TDNA_ALT_LB_END = NA,
         TDNA_ALT_RB_END = NA,
         TDNA_ALT_IS_LBRB = NA,
         TANDEM_DUPLICATION = NA,
         FLANK_B_TDNA_SIDE = NA)

GLOBAL.wb = left_join(GLOBAL.wb_pre, GLOBAL.missing_columns)

#remove previously marked problematic events as well as duplicate positions
if (GLOBAL.REMOVEPROBLEMS == TRUE) {
  funlog("removing problematic events")
  GLOBAL.total_data_positioncompare_pre = GLOBAL.wb %>%
    #first remove problematic events based on characteristics of the events themselves
    filter(AnchorCount >= GLOBAL.ANCHORCUTOFF,
           ANCHOR_DIST >= GLOBAL.MINANCHORDIST,
           ANCHOR_DIST <= MAXANCHORDIST,
           Consensus_freq >= 0.75,
           MATE_FLANK_B_CHROM_AGREE == TRUE,
           FAKE_DELIN_CHECK == FALSE) %>% 
    #then sort the data based on genomic location
    arrange(Alias, FLANK_B_CHROM, FLANK_B_START_POS)%>%
    #then check whether for each event, the one on the previous row is close by (within 10bp).
    mutate(previous_pos = lag(FLANK_B_START_POS),
           prev_Alias = lag(Alias),
           prev_chrom = lag(FLANK_B_CHROM),
           prev_orient = lag(FLANK_B_ISFORWARD))%>%
    mutate(pos_dif = abs(previous_pos - FLANK_B_START_POS))%>%
    mutate(Alias_compare = if_else(prev_Alias == Alias,
                                    TRUE,
                                    FALSE))%>%
    mutate(chrom_compare = if_else(prev_chrom == FLANK_B_CHROM,
                                   TRUE,
                                   FALSE))%>%
    mutate(orient_compare = if_else(prev_orient == FLANK_B_ISFORWARD,
                                    TRUE,
                                    FALSE))%>%
    mutate(same_as_prev = if_else(pos_dif < 11 & chrom_compare==TRUE & Alias_compare == TRUE & orient_compare == TRUE,
                                  TRUE,
                                  FALSE))%>%
    mutate(ID = 0)%>%
    mutate(Family=NULL)
  

  
  
  #assign IDs, each representing a separate event (meaning that are not too close to be considered the same event)
  GLOBAL.ID_prev = 0
  for (i in 2:length(GLOBAL.total_data_positioncompare_pre$Alias)){
    if (GLOBAL.total_data_positioncompare_pre$same_as_prev[i] == TRUE){
      GLOBAL.total_data_positioncompare_pre$ID[i] = GLOBAL.ID_prev
      GLOBAL.ID_prev = GLOBAL.total_data_positioncompare_pre$ID[i]
    }else{
      GLOBAL.total_data_positioncompare_pre$ID[i] = GLOBAL.ID_prev + 1
      GLOBAL.ID_prev = GLOBAL.total_data_positioncompare_pre$ID[i]
    }
  }
  
  #add family info
  GLOBAL.sample_info2  = GLOBAL.sample_info %>%
    select(Family, RunID, Sample)%>%
    rowwise()%>%
    mutate(Alias = paste(Sample, RunID, sep="_"))%>%
    select(Alias, Family)%>%
    ungroup()
  
  GLOBAL.total_data_positioncompare = left_join(GLOBAL.sample_info2, GLOBAL.total_data_positioncompare_pre, by=c("Alias"))%>%
    filter(!is.na(Family) & !is.na(Name))
  
  funlog("combining junctions with similar positions")
  #combine junctions with similar positions and get the characteristics of the consensus event from the event the most anchors 
  GLOBAL.total_data_near_positioncombined = GLOBAL.total_data_positioncompare %>%
    group_by(Alias, FLANK_B_CHROM, Plasmid, FLANK_B_ISFORWARD, DNASample, Subject, ID, Focus_contig, Genotype, Ecotype, Plasmid_alt, Family, FlankAUltEnd, AgroGeno, Species, RemoveNonTranslocation, GroupSamePosition, Translocation, Translocation_del_resolved, TANDEM_DUPLICATION, DSB_FW_END, DSB_OVERHANG, DSB_CONTIG, FLANK_B_TDNA_SIDE, Experiment)%>%
    summarize(AnchorCountSum = sum(AnchorCount), 
              ReadCountSum = sum(ReadCount),
              FLANK_B_START_POS_CON = as.integer(FLANK_B_START_POS[which.max(AnchorCount)]),
              delRelativeStart_CON = as.integer(delRelativeStart[which.max(AnchorCount)]),
              delRelativeStartTD_CON = as.integer(delRelativeStartTD[which.max(AnchorCount)]),
              delRelativeEnd_CON = as.integer(delRelativeEnd[which.max(AnchorCount)]),
              delRelativeEndTD_CON = as.integer(delRelativeEndTD[which.max(AnchorCount)]),
              ANCHOR_DIST_CON = as.integer(ANCHOR_DIST[which.max(AnchorCount)]),
              FILLER_CON = FILLER[which.max(AnchorCount)],
              MH_CON = MH[which.max(AnchorCount)],
              insSize_CON = as.integer(insSize[which.max(AnchorCount)]),
              homologyLength_CON = as.integer(homologyLength[which.max(AnchorCount)]),
              delSize_CON = as.integer(delSize[which.max(AnchorCount)]),
              Type_CON = Type[which.max(AnchorCount)],
              SEQ_1_con_CON = SEQ_1_con[which.max(AnchorCount)],
              SEQ_2_con_CON = SEQ_2_con[which.max(AnchorCount)],
              TRIM_LEN_CON = as.integer(TRIM_LEN[which.max(AnchorCount)]),
              RunID_CON = RunID[which.max(AnchorCount)],
              SEQ_1_LEN_max_CON = max(SEQ_1_LEN_max),
              .groups="drop"
              )%>%
    rename(FLANK_B_START_POS = FLANK_B_START_POS_CON,
           delRelativeStart = delRelativeStart_CON,
           delRelativeStartTD = delRelativeStartTD_CON,
           delRelativeEnd = delRelativeEnd_CON,
           delRelativeEndTD = delRelativeEndTD_CON,
           ANCHOR_DIST = ANCHOR_DIST_CON,
           FILLER = FILLER_CON,
           MH = MH_CON,
           insSize = insSize_CON,
           homologyLength = homologyLength_CON,
           delSize = delSize_CON,
           Type = Type_CON,
           SEQ_1_con = SEQ_1_con_CON,
           SEQ_2_con = SEQ_2_con_CON,
           RunID = RunID_CON,
           TRIM_LEN = TRIM_LEN_CON,
           ReadCount = ReadCountSum,
           SEQ_1_LEN_max = SEQ_1_LEN_max_CON)

  #get a list of families
  GLOBAL.wb_family = GLOBAL.sample_info2 %>% select(Family) %>% distinct() %>% filter(Family!=0)
  if (nrow(GLOBAL.wb_family) == 0) {
    funlog("no family info detected")
    #if no families are indicated

    GLOBAL.wb_flag = GLOBAL.total_data_near_positioncombined %>%
        #then examine positions across samples and remove those that occur multiple times
        group_by(FLANK_B_START_POS) %>%
        mutate(duplicate_position = if_else(n() > 1 & Focus_contig != FLANK_B_CHROM,
                                            TRUE,
                                            FALSE)) %>%
        ungroup() %>%
        filter(duplicate_position == FALSE)
      
  
  } else{
    funlog("taking family into consideration")
    #take families into account
    GLOBAL.wb_filter_total = GLOBAL.total_data_near_positioncombined %>% filter(Family == 99999999) #make an empty file
    
    for (i in GLOBAL.wb_family$Family) {
      funlog(paste0("Cleanup family ", i))
      #cleanup per family
      GLOBAL.wb_current_family = GLOBAL.total_data_near_positioncombined %>% filter(Family == i) %>% select(Alias) %>% distinct() #make a list of aliases within the current family
      GLOBAL.wb_filter_subtotal = GLOBAL.total_data_near_positioncombined %>% filter(Family == 99999999) #make an empty file
      
      for (j in GLOBAL.wb_current_family$Alias) {
        #per alias in that family
        
        GLOBAL.wb_filter_current = GLOBAL.total_data_near_positioncombined %>%
          filter(Family != i | Alias == j) %>% #events are either not of the current family, or they belong to the current alias
          group_by(FLANK_B_START_POS) %>%
          mutate(duplicate_position = if_else(n() > 1 & Focus_contig != FLANK_B_CHROM,
                                              TRUE,
                                              FALSE)) %>%
          
          ungroup() %>%
          filter(duplicate_position == FALSE &
                   Family == i) #keep only events belonging to the current file, and remove duplicate positions
        
        GLOBAL.wb_filter_subtotal = rbind(GLOBAL.wb_filter_subtotal, GLOBAL.wb_filter_current) #combine surviving events from the current family
        
      }
      GLOBAL.wb_filter_total = rbind(GLOBAL.wb_filter_total, GLOBAL.wb_filter_subtotal) #combining surviving events from all families
      
    }
    funlog("Fetching non-duplicate position events not belonging to a family")
    GLOBAL.wb_nonfamily = GLOBAL.total_data_near_positioncombined %>%
      filter(Family == 0) %>%
      group_by(FLANK_B_START_POS) %>%
      mutate(duplicate_position = if_else(n() > 1 & Focus_contig != FLANK_B_CHROM,
                                          TRUE,
                                          FALSE)) %>%
      
      ungroup() %>%
      filter(duplicate_position == FALSE)
    funlog("combining surviving family and nonfamily events")
    GLOBAL.wb_flag = rbind(GLOBAL.wb_filter_total, GLOBAL.wb_nonfamily)  
    
  }
} else{
  funlog("flagging problems only")
  GLOBAL.wb_flag = GLOBAL.wb %>%
    group_by(FLANK_B_START_POS) %>%
    mutate(duplicate_position = if_else(n() > 1 & Focus_contig != FLANK_B_CHROM,
                                        TRUE,
                                        FALSE)) %>%
    
    ungroup() 
}

GLOBAL.wb_flag1 = GLOBAL.wb_flag %>%
  #add the options
  mutate(RemoveProblematicEvents = GLOBAL.REMOVEPROBLEMS,
         MINANCHORDIST = GLOBAL.MINANCHORDIST,
         ANCHORCUTOFF = GLOBAL.ANCHORCUTOFF)%>%
  #add/ change several things for compatibility with SIQplotteR
  mutate(getHomologyColor = "dummy",
         Barcode = "dummy")%>%
  rename(countEvents = ReadCount,
         insertion = FILLER,
         homology = MH)


funlog("Writing output")
GLOBAL.work_book2 <- createWorkbook()
addWorksheet(GLOBAL.work_book2, "rawData")
writeData(GLOBAL.work_book2, sheet = 1, GLOBAL.wb_flag1)

#Write an additional sheet with read number info
GLOBAL.read_numbers_info = read.csv(paste0(GLOBAL.input_dir, "read_numbers.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
GLOBAL.wb_numbers = GLOBAL.read_numbers_info %>% 
  mutate(Alias = paste0(Sample, "_", RunID))%>%
  mutate(Sample = NULL,
         RunID = NULL)
addWorksheet(GLOBAL.work_book2, "readNumbers")
writeData(GLOBAL.work_book2, sheet = 2, GLOBAL.wb_numbers)

#write another sheet with all the messages
addWorksheet(GLOBAL.work_book2, "runLog")
writeData(GLOBAL.work_book2, sheet = 3, runlog)

saveWorkbook(GLOBAL.work_book2, file = paste0(GLOBAL.output_dir, "Data_combined_CISTRANSGUIDE_V2_", as.integer(Sys.time()), ".xlsx"), overwrite = TRUE)

message("CISTRANSGUIDE analysis has completed")
