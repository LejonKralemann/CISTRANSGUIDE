###############################################################################
#install and load packages
###############################################################################
if (require(Biostrings)==FALSE){install.packages("Biostrings", repos = "http://cran.us.r-project.org")}
if (require(stringi)==FALSE){install.packages("stringi", repos = "http://cran.us.r-project.org")}
if (require(stringdist)==FALSE){install.packages("stringdist", repos = "http://cran.us.r-project.org")}
if (require(tidyverse)==FALSE){install.packages("tidyverse", repos = "http://cran.us.r-project.org")}
if (require(openxlsx)==FALSE){install.packages("openxlsx", repos = "http://cran.us.r-project.org")}


###############################################################################
#set parameters - adjustable
###############################################################################
GLOBAL.input_dir= "./input/"
GLOBAL.output_dir= "./output/"
GLOBAL.GROUPSAMEPOS=FALSE #if true, it combines reads with the same genomic pos, which helps in removing artefacts. Typically used for TRANSGUIDE, but disabled for CISGUIDE.
GLOBAL.REMOVENONTRANS=FALSE #if true, it only considers translocations. Typically used for TRANSGUIDE, but disabled for CISGUIDE. Note that some translocations on the same chromosome will also be removed thusly.
GLOBAL.REMOVEPROBLEMS=FALSE #if true it removes all problematic reads from the combined datafile. Note if this is false, no duplicate filtering will be performed, because first reads due to barcode hopping need to be removed by removing events with few anchors.
GLOBAL.ANCHORCUTOFF=3 #each event needs to have at least this number of anchors, otherwise it is marked as problematic (and potentially removed) 
GLOBAL.MINANCHORDIST=150 #should be matching a situation where the mate is 100% flank B.
GLOBAL.MAXANCHORDIST=2000 #the furthest position that the mate anchor can be, except on T-DNA.
GLOBAL.FLANKBEYONDDSB=5000 #how much flank A and flank B are allowed to continue beyond the DSB (not applicable when the focus contig is the T-DNA)
GLOBAL.MINLEN=90 #this is the minimal read length. if you write NA here, then the software will calculate the minimal read length based on the distance to nick/dsb and FLANK_B_LEN_MIN. Should be at the very least 60bp, but 90bp is more common to have as minimum.
GLOBAL.LB_SEQUENCES = c("TGGCAGGATATATTGTGGTGTAAAC", "CGGCAGGATATATTCAATTGTAAAT") #the nick is made after the 3rd nt
GLOBAL.RB_SEQUENCES = c("TGACAGGATATATTGGCGGGTAAAC", "TGGCAGGATATATGCGGTTGTAATT") #the nick is made after the 3rd nt
GLOBAL.TD_SIZE_CUTOFF = 6 #the smallest TD that is considered as TD (*with regards to the Type variable). Any smaller TD is considered merely an insertion.
GLOBAL.FASTA_MODE = TRUE #Typically false, if TRANSGUIDE/CISGUIDE library prep and illumina sequencing has been done. TRUE if sequences from another source are being analyzed with this program.

###############################################################################
#set parameters - non-adjustable
###############################################################################
GLOBAL.NF_NUMBER = as.integer(-99999999) #don't change
GLOBAL.ERROR_NUMBER = as.integer(99999999) #don't change
GLOBAL.hash=system("git rev-parse HEAD", intern=TRUE)
GLOBAL.hash_little=substr(hash, 1, 8)
GLOBAL.sample_info = read.csv(paste0(input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
GLOBAL.TIME_START=round(as.numeric(Sys.time())*1000, digits=0)
GLOBAL.FLANK_B_LEN_MIN = 30 #minimum length of flank B. Also determines the size of DSB_AREA_SEQ. do not change because the preprocessing program will still be set at 30.
GLOBAL.MINMAPQUALA = 42 #minimum mapping quality (phred). 42 means a perfect, unambiguous match (well, should be)
GLOBAL.MINMAPQUALB = 42 #same, but for flank B
GLOBAL.PercentageDone = 0 #var for indicating progress
GLOBAL.TotalFileSize = 0

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

for (i in row.names(sample_info)){
  FILE.Sample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  FILE.RunID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(RunID))
  
  if (file.exists(paste0(input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))==FALSE){
    next
  }else{
    GLOBAL.TotalFileSize = GLOBAL.TotalFileSize + (file.info((paste0(input_dir, FILE.Sample, "_", FILE.RunID, "_A.txt"))))$size
  }
}

funlog(paste0("Total file size to process: ", GLOBAL.TotalFileSize, " bytes"))

###############################################################################
#Process data: step 1
#checking the file and reading metadata
###############################################################################

funlog("Checking file and reading metadata")

#process all files
for (i in row.names(sample_info)){
  
  ####################  pre-cleanup  #####################
  objects_to_clean=c("data", "data_improved_a", "data_improved_b", "data_improved_c", "data_improved_1", "data_improved_2", "data_improved_3", "data_improved_3b", "data_improved_4", "data_improved_5", "data_improved_5b", "data_improved_6", "data_improved_8pre2", "data_improved_8", "data_improved_9", "data_improved_10", "AgroGeno", "contig_seq", "DNASample", "DSB_AREA_SEQ", "DSB_AREA_SEQ_RC", "DSB_CONTIG", "DSB_FW_END", "DSB_OVERHANG", "Ecotype", "FLANK_A_ISFORWARD", "FLANK_A_REF_GLOBAL", "FlankAUltEnd", "FOCUS_CONTIG", "FOCUS_LOCUS", "Genotype", "GLOBAL_TOTAL_REF", "Library", "PLASMID", "PLASMID_ALT", "plasmid_alt_seq", "plasmid_seq", "Primer_match", "Primer_match_3", "Primer_RC_match", "Primer_match_perfect", "Primer_pos", "Primer_seq", "Primer_seq_len", "PRIMER_TO_DSB_GLOBAL", "REF", "RunID", "Sample", "TDNA_ALT_IS_LBRB", "TDNA_ALT_LB_END", "TDNA_ALT_LB_FW", "TDNA_ALT_RB_END", "TDNA_ALT_RB_FW", "TDNA_IS_LBRB", "TDNA_LB_END,TDNA_LB_FW", "TDNA_RB_END", "TDNA_RB_FW", "genomeseq", "LB_match", "LB_match_RV", "LB2_match", "LB2_match_RV", "RB_match", "RB_match_RV", "RB2_match", "RB2_match_RV")
  for (z in objects_to_clean){
    assign(z, NULL)
  }
 
  ####################  general variables acquired from the information sheet  #####################
  
  FILE.Sample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  FILE.RunID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(RunID))
  
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
  if (nrow(data)==0){
    funlog("Primary processed file empty, moving to the next sample")
    next
  }
  
  ####################  continue general variables acquired from the information sheet  #####################
  
  FILE.DSB_CONTIG = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_CONTIG))#chromosome name or NA
  FILE.Genotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Genotype))
  FILE.PLASMID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid))
  FILE.PLASMID_ALT = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid_alt))
  if (is.na(FILE.PLASMID_ALT) | FILE.PLASMID_ALT=="NA"){FILE.PLASMID_ALT=""} #to avoid doing the following check: FLANK_B_CHROM==NA
  
  FILE.REF = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ref))
  FILE.DNASample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DNA))
  FILE.Ecotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ecotype))
  FILE.Library = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  FILE.AgroGeno = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(AgroGeno))
  FILE.DSB_FW_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_FW_END)) #end of left flank before DSB
  FILE.FOCUS_CONTIG = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Focus_contig_name))#Same as DSB_contig, or same as plasmid
  FILE.LOCUS_NAME = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Locus_name))#LB or RB if TRANSGUIDE, or a name of a locus if CISGUIDE
  FILE.Primer_seq = str_replace_all(toupper(as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Primer))), "TCAGACGTGTGCTCTTCCGATCT", "")
  
  #the following variables can be supplied via the sample information sheet, but NA is also allowed, then the software will look.
  FILE.TDNA_LB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_LB_END))
  FILE.TDNA_RB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_RB_END))
  FILE.TDNA_ALT_LB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_ALT_LB_END))
  FILE.TDNA_ALT_RB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_ALT_RB_END))
  
  
  
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
  
  FILE.plasmid_seq = as.character(eval(parse(text = paste0("genomeseq$`", FILE.PLASMID, "`"))))
  
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
    
    for (i in c("LB_match", "LB_match_RV", "LB2_match", "LB2_match_RV")){
      if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
        
        
        if (i=="LB_match" | i=="LB2_match"){
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
    FILE.TDNA_RB_END = NA
    FILE.TDNA_RB_FW = NA
    
    for (i in c("RB_match", "RB_match_RV", "RB2_match", "RB2_match_RV")){
      if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
        
        
        if (i=="RB_match" | i=="RB2_match"){
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
      FILE.plasmid_alt_seq = as.character(eval(parse(text = paste0("genomeseq$`", FILE.PLASMID_ALT, "`"))))
      
      #find the LB
      FILE.LB_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[1]], subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.LB_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[1]]))), subject = DNAString(FILE.plasmid_alt_seq), max.mismatch = 0, fixed=TRUE))
      FILE.LB2_match = as.data.frame(matchPattern(pattern = GLOBAL.LB_SEQUENCES[[2]],subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      FILE.LB2_match_RV = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(GLOBAL.LB_SEQUENCES[[2]]))),subject = DNAString(FILE.plasmid_alt_seq),max.mismatch = 0,fixed = TRUE))
      
      
      for (i in c("LB_match", "LB_match_RV", "LB2_match", "LB2_match_RV")){
        if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
          
          
          if (i=="LB_match" | i=="LB2_match"){
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
      
      
      for (i in c("RB_match", "RB_match_RV", "RB2_match", "RB2_match_RV")){
        if (nrow(get(i)) > 0 & nrow(get(i)) < 2) {
          
          
          if (i=="RB_match" | i=="RB2_match"){
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
  
  
  FILE.contig_seq = as.character(eval(parse(text = paste0("genomeseq$`", FILE.FOCUS_CONTIG, "`"))))
  if (length(FILE.contig_seq)==0){
    funlog("Focus contig not found in reference fasta. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
    next 
  }
  
  ####################  Primer check  #####################
  
  if (GLOBAL.FASTA_MODE==FALSE){  
    
  #this checks whether primer can be found, checks what the orientation of flank A is, whether the primer matches the plasmid or not
  FILE.Primer_match = as.data.frame(matchPattern(pattern = FILE.Primer_seq, subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE))
  FILE.Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(FILE.Primer_seq))), subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE))
  
  if (FILE.FOCUS_CONTIG == FILE.Plasmid){
  if (nrow(FILE.Primer_match) > 0 & nrow(FILE.Primer_match) < 2){
    FILE.Primer_pos = as.numeric(FILE.Primer_match$start)
    FILE.Primer_match_perfect=TRUE
    if (FILE.Primer_pos > FILE.TDNA_LB_END & FILE.Primer_pos < FILE.TDNA_RB_END){
      FILE.Primer_on_TDNA=TRUE
      FILE.FLANK_A_ISFORWARD=TRUE
      if (FILE.TDNA_IS_LBRB == TRUE){
      FILE.FOCUS_LOCUS="RB"
      FILE.FlankAUltEnd = FILE.TDNA_RB_END  
      }else{
        FILE.FOCUS_LOCUS="LB"
        FILE.FlankAUltEnd = FILE.TDNA_LB_END
      }
    }else{
      FILE.Primer_on_TDNA=FALSE
      FILE.FOCUS_LOCUS="OTHER"
    }
  }else  if (nrow(FILE.Primer_RC_match) > 0 & nrow(FILE.Primer_RC_match) < 2){
    FILE.Primer_pos = as.numeric(FILE.Primer_match$start)
    FILE.Primer_match_perfect=TRUE
    if (FILE.Primer_pos > FILE.TDNA_LB_END & FILE.Primer_pos < FILE.TDNA_RB_END){
      FILE.Primer_on_TDNA=TRUE
      FILE.FLANK_A_ISFORWARD=FALSE
      if (FILE.TDNA_IS_LBRB == TRUE){
        FILE.FOCUS_LOCUS="LB"
        FILE.FlankAUltEnd = FILE.TDNA_LB_END
      }else{
        FILE.FOCUS_LOCUS="RB"
        FILE.FlankAUltEnd = FILE.TDNA_RB_END  
      }
    }else{
      FILE.Primer_on_TDNA=FALSE
      FILE.FOCUS_LOCUS="OTHER"
    }
    }else if (nrow(FILE.Primer_match) >1 | nrow(FILE.Primer_RC_match) >1){
    funlog("Primer found several times on this contig")
    next
  }else{
    funlog("Primer not found")
    next
  }
  }else{
    FILE.FOCUS_LOCUS=FILE.LOCUS_NAME
    FILE.Primer_on_TDNA=FALSE
    
      if (nrow(FILE.Primer_match) > 0 & nrow(FILE.Primer_match) < 2){
        FILE.FLANK_A_ISFORWARD=TRUE
        FILE.FlankAUltEnd = FILE.DSB_FW_END
      }else if (nrow(FILE.Primer_RC_match) > 0 & nrow(FILE.Primer_RC_match) < 2){
        FILE.FLANK_A_ISFORWARD=FALSE
        FILE.FlankAUltEnd = FILE.DSB_FW_END+1
      }else if (nrow(FILE.Primer_match) >1 | nrow(FILE.Primer_RC_match) >1){
        funlog("Primer found several times on this contig")
        next
      }else{
        funlog("Primer not found")
        next
      }
  }
  
  ####################  acquire the DSB AREA SEQ  #####################
  
 
  FILE.DSB_AREA_SEQ = (if (FILE.FLANK_A_ISFORWARD == TRUE){
    substr(FILE.contig_seq, start= FILE.FlankAUltEnd - ((GLOBAL.FLANK_B_LEN_MIN/2)-1), stop= FILE.FlankAUltEnd + (GLOBAL.FLANK_B_LEN_MIN/2))
    }else if (FILE.FLANK_A_ISFORWARD == FALSE){
    as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= FILE.FlankAUltEnd - (GLOBAL.FLANK_B_LEN_MIN/2), stop= FILE.FlankAUltEnd + ((GLOBAL.FLANK_B_LEN_MIN/2)-1)))))
      }else{
    ""
        })
  if (FILE.DSB_AREA_SEQ == ""){
    FILE.EXECUTE_DSB_AREA_CHECK=FALSE
  }else{
    FILE.EXECUTE_DSB_AREA_CHECK=TRUE
  }
  
  FILE.DSB_AREA_SEQ_RC = as.character(reverseComplement(DNAString(FILE.DSB_AREA_SEQ)))
    
  ####################  continue calculated general variables  #####################  
    
  #calculate the distance from primer to DSB/nick
  FILE.PRIMER_TO_DSB_GLOBAL = if (FILE.FLANK_A_ISFORWARD == TRUE){
    FILE.FlankAUltEnd - (FILE.Primer_pos -1)
  }else{
    FILE.Primer_pos - (FILE.FlankAUltEnd -1)
  }

  #usually you don't want the primer to be so far away. But sometimes when you cut away the end of the T-DNA for instance, then you may want to keep the T-DNA end position.
  if (FILE.PRIMER_TO_DSB_GLOBAL>300){
    funlog("Warning: primer is more than 300 bp away from the indicated end of FLANK A. Continuing anyway.")
  }

 
  #get the start and stop positions for the ref sequence, in case it goes beyond the end of the chromosome
  FILE.TOTAL_REF_START = if(FILE.FlankAUltEnd-FILE.FLANKBEYONDDSB < 1){
    1
  }else{
    FILE.FlankAUltEnd-FILE.FLANKBEYONDDSB
  }
  FILE.TOTAL_REF_STOP = if(FILE.FlankAUltEnd+FILE.FLANKBEYONDDSB > nchar(FILE.contig_seq)){
    nchar(FILE.contig_seq)
  }else{
    FILE.FlankAUltEnd+FILE.FLANKBEYONDDSB
  }
  
  #get the REF seq for flank A.
  FILE.FLANK_A_REF_GLOBAL = if (FILE.FLANK_A_ISFORWARD == TRUE){
    substr(FILE.contig_seq, start= FILE.Primer_pos, stop= FILE.TOTAL_REF_STOP)
  }else{
    as.character(reverseComplement(DNAString(
      substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.Primer_pos))
    ))
  }
  #get the REF seq that includes both flank a and b, later used to get the ref for flank b  
  FILE.TOTAL_REF_GLOBAL = if (FILE.FLANK_A_ISFORWARD == TRUE){substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.TOTAL_REF_STOP)
  }else{
    as.character(reverseComplement(DNAString(substr(FILE.contig_seq, start= FILE.TOTAL_REF_START, stop= FILE.TOTAL_REF_STOP))))
  }
  }
  
  function_time("Step 1 took ")

  ###############################################################################
  #Process data: step 2
  ###############################################################################
  
  if (FASTA_MODE == TRUE){
    FILE.data_pre = FILE.data  %>%
    rowwise()%>%
    mutate(PRIMER_SEQ = substr(SEQ_1, 1, 30))%>%
    mutate(Primer_match = as.data.frame(matchPattern(pattern = PRIMER_SEQ, subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE)))%>%
    mutate(Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(PRIMER_SEQ))), subject = DNAString(FILE.contig_seq), max.mismatch = 0, fixed=TRUE))
    )
    
    for (j in row.names(FILE.data_pre)){
      j_int=as.integer(j)
    if (FILE.FOCUS_CONTIG == FILE.Plasmid){
      if (nrow(FILE.data_pre$Primer_match[[j_int]]) > 0 & nrow(FILE.data_pre$Primer_match[[j_int]]) < 2){
        FILE.data_pre$Primer_pos[[j_int]] = as.numeric(FILE.data_pre$Primer_match$start[[j_int]])
        FILE.data_pre$Primer_match_perfect[[j_int]]=TRUE
        if (FILE.data_pre$Primer_pos[[j_int]] > FILE.data_pre$TDNA_LB_END[[j_int]] & FILE.data_pre$Primer_pos[[j_int]] < FILE.data_pre$TDNA_RB_END[[j_int]]){
          FILE.data_pre$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data_pre$FLANK_A_ISFORWARD[[j_int]]=TRUE
          if (FILE.TDNA_IS_LBRB == TRUE){
            FILE.data_pre$FOCUS_LOCUS[[j_int]]="RB"
            FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.TDNA_RB_END  
          }else{
            FILE.data_pre$FOCUS_LOCUS[[j_int]]="LB"
            FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.TDNA_LB_END
          }
        }else{
          FILE.data_pre$Primer_on_TDNA[[j_int]]=FALSE
          FILE.data_pre$FOCUS_LOCUS[[j_int]]="OTHER"
        }
      }else  if (nrow(FILE.data_pre$Primer_RC_match[[j_int]]) > 0 & nrow(FILE.data_pre$Primer_RC_match[[j_int]]) < 2){
        FILE.data_pre$Primer_pos[[j_int]] = as.numeric(FILE.data_pre$Primer_match[[j_int]]$start)
        FILE.data_pre$Primer_match_perfect[[j_int]]=TRUE
        if (FILE.data_pre$Primer_pos[[j_int]] > FILE.TDNA_LB_END & FILE.data_pre$Primer_pos[[j_int]] < FILE.TDNA_RB_END){
          FILE.data_pre$Primer_on_TDNA[[j_int]]=TRUE
          FILE.data_pre$FLANK_A_ISFORWARD[[j_int]]=FALSE
          if (FILE.TDNA_IS_LBRB == TRUE){
            FILE.data_pre$FOCUS_LOCUS[[j_int]]="LB"
            FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.TDNA_LB_END
          }else{
            FILE.data_pre$FOCUS_LOCUS[[j_int]]="RB"
            FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.TDNA_RB_END  
          }
        }else{
          FILE.data_pre$Primer_on_TDNA[[j_int]]=FALSE
          FILE.data_pre$FOCUS_LOCUS[[j_int]]="OTHER"
        }
      }else if (nrow(FILE.data_pre$Primer_match[[j_int]]) >1 | nrow(FILE.data_pre$Primer_RC_match[[j_int]]) >1){
        funlog("Primer found several times on this contig")
        next
      }else{
        funlog("Primer not found")
        next
      }
    }else{
      FILE.data_pre$FOCUS_LOCUS[[j_int]]=FILE.LOCUS_NAME
      FILE.data_pre$Primer_on_TDNA[[j_int]]=FALSE
      
      if (nrow(FILE.data_pre$Primer_match[[j_int]]) > 0 & nrow(FILE.data_pre$Primer_match[[j_int]]) < 2){
        FILE.data_pre$FLANK_A_ISFORWARD[[j_int]]=TRUE
        FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.DSB_FW_END
      }else if (nrow(FILE.data_pre$Primer_RC_match[[j_int]]) > 0 & nrow(FILE.data_pre$Primer_RC_match[[j_int]]) < 2){
        FILE.data_pre$FLANK_A_ISFORWARD[[j_int]]=FALSE
        FILE.data_pre$FlankAUltEnd[[j_int]] = FILE.DSB_FW_END+1
      }else if (nrow(FILE.data_pre$Primer_match[[j_int]]) >1 | nrow(FILE.data_pre$Primer_RC_match[[j_int]]) >1){
        funlog("Primer found several times on this contig")
        next
      }else{
        funlog("Primer not found")
        next
      }
    }
    } #continue here, in fasta mode the things from the previous sections need to be calculated here, line by line
    #other things to change still: there are still many places where I forgot to add "GLOBAL." or "FILE.". 
    #there may also be places where I dont want to add these things, but instead create a copy of this variable within the tibble.
 
    
  }else{
    FILE.data_pre = FILE.data %>%
    mutate(FLANK_A_ISFORWARD = FILE.FLANK_A_ISFORWARD) 
  }
  

  FILE.data_improved_a  = FILE.data_pre %>%
    
    
    #filter(QNAME == "A00379:436:H3CHWDMXY:1:2127:30852:25300")%>%
    
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
    FILE.data_improved_b = FILE.data_improved_a %>%
      filter(FLANK_B_CHROM != FILE.FOCUS_CONTIG)
  }else{
    FILE.data_improved_b = FILE.data_improved_a
  }
  
  
  FILE.data_improved_c = FILE.data_improved_b %>%
    
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
                                      TRUE ~ ERROR_NUMBER))  %>%
    
    #test whether B flanks from read and mate are mapped to same chromosome, and agree on orientation.
    #note that the two reads are in two different orientations, so the B orientation should not agree
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_ISFORWARD != MATE_B_ISFORWARD,
                                              TRUE,
                                              FALSE)) 
  
  
  
  
  
  FILE.data_improved1 = FILE.data_improved_c %>%
    
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
      MATE_FLANK_B_CHROM) %>%
    arrange(desc(AvgBaseQual_1), desc(AvgBaseQual_2)) %>%
    summarize(QNAME_first = dplyr::first(QNAME),
              SEQ_2_first = dplyr::first(SEQ_2),
              MATE_B_END_POS_list = toString(MATE_B_END_POS),
              .groups="drop"
              )%>%

    #then remove reads that don't begin with the primer
    mutate(PRIMER_SEQ_LEN = nchar(PRIMER_SEQ)) %>%
    mutate(SEQ_1_start = substr(SEQ_1, 1, PRIMER_SEQ_LEN)) %>%
    filter(SEQ_1_start == PRIMER_SEQ) %>%
    
    #the lines below until the grouping are under investigation. Maybe doesn't work well yet. First make the code so that the output is identical to the output before the changes of a non fasta sample. Then test the fasta data.
    mutate(PRIMER_SEQ = if (FASTA_MODE == FALSE){
      PRIMER_SEQ
      }else{
        substr(SEQ_1, 1, 30)
      }) %>%
    mutate(PRIMER_POS_FAKE_match = if (FASTA_MODE == FALSE){
      NA
    }else{
      if (FLANK_A_ISFORWARD == TRUE){
        list(matchPattern(DNAString(PRIMER_SEQ), DNAString(contig_seq), max.mismatch = 0))
      }else{
        list(matchPattern(as.character(reverseComplement(DNAString(PRIMER_SEQ))), DNAString(contig_seq), max.mismatch = 0))
      }}) %>%
    mutate(PRIMER_POS_FAKE = if (FASTA_MODE == FALSE){
      as.integer(ERROR_NUMBER)
    }else{
      if (FLANK_A_ISFORWARD == TRUE){
      as.integer(PRIMER_POS_FAKE_match@ranges@start)
      }else{
      as.integer(PRIMER_POS_FAKE_match@ranges@start)+29
      }
    }) %>%
    mutate(PRIMER_TO_DSB = if (FASTA_MODE == FALSE){
      as.integer(PRIMER_TO_DSB_GLOBAL)
    }else{
      if (FLANK_A_ISFORWARD == TRUE){
        FlankAUltEnd - (PRIMER_POS_FAKE -1)
      }else{
        PRIMER_POS_FAKE - (FlankAUltEnd -1)
      }}) %>%
    mutate(FLANK_A_REF = if (FASTA_MODE == FALSE){
      FLANK_A_REF_GLOBAL
    }else{
      if (FLANK_A_ISFORWARD == TRUE){
        substr(contig_seq, start= PRIMER_POS_FAKE, stop= FlankAUltEnd)
      }else{
        as.character(reverseComplement(DNAString(substr(contig_seq, start= FlankAUltEnd, stop= PRIMER_POS_FAKE))))
      }}) %>%
    mutate(TOTAL_REF = if (FASTA_MODE == FALSE){
      TOTAL_REF_GLOBAL
    }else{
      if (FLANK_A_ISFORWARD == TRUE){
        substr(contig_seq, start= PRIMER_POS_FAKE, stop= PRIMER_POS_FAKE+MAX_DIST_FLANK_B_END+PRIMER_TO_DSB)
      }else{
        as.character(reverseComplement(DNAString(
          substr(contig_seq, start= PRIMER_POS_FAKE-(MAX_DIST_FLANK_B_END+PRIMER_TO_DSB), stop= PRIMER_POS_FAKE))))
      }}) %>%
    mutate(TOTAL_REF_LEN = nchar(TOTAL_REF))%>%
    ungroup() 
  
  #check if any reads have survived
  if (nrow(FILE.data_improved1)==0){
    funlog(paste0("No reads surviving for sample ", FILE.Library))
    #show progress
    GLOBAL.PercentageDone = GLOBAL.PercentageDone + ((FILE.CurrentFileSize/GLOBAL.TotalFileSize)*100)
    funlog(paste0("CISTRANSGUIDE analysis ", round(GLOBAL.PercentageDone, digits=3), "% complete"))
    next
  }
  
  function_time("Step 2 took ")
  
  ###############################################################################
  #Process data: step 3
  ###############################################################################
  
  FILE.data_improved2 = FILE.data_improved1 %>%
    
    #check for expected position and orientation base on primer seqs
    rowwise()%>%
    mutate(FLANK_A_START_POS =
      if (FASTA_MODE == FALSE){
        Primer_pos
      }else{
        PRIMER_POS_FAKE
      }
    )%>%
    ungroup() %>%
    
    #calculate the end position of the B flank, in order to get the ref seq. This is the end the furthest away from the junction.
    mutate(FLANK_B_END_POS = case_when(FLANK_B_ISFORWARD==FALSE ~ as.integer(B_POS),
                                       FLANK_B_ISFORWARD==TRUE ~ as.integer(B_POS + 29),
                                       TRUE ~ ERROR_NUMBER))
    

  function_time("Step 3 took ")
  
  ###############################################################################
  #Process data: step 4
  ###############################################################################
  
  FILE.data_improved2b = FILE.data_improved2 %>%

    #Find how much SEQ_1 matches with FLANK_A_REF. Allow 1 bp mismatch somewhere, if the alignment after that continues for at least another 10 bp.
    rowwise() %>%
    mutate(FLANK_A_MATCH = matcher_skipper(FLANK_A_REF, SEQ_1)) %>%
    ungroup() %>%
    mutate(FLANK_A_LEN = nchar(FLANK_A_MATCH)) %>%
    mutate(FLANK_A_END_POS = case_when(FLANK_A_ISFORWARD == TRUE ~ as.integer(FLANK_A_START_POS + (FLANK_A_LEN -1)),
                                       FLANK_A_ISFORWARD == FALSE ~ as.integer(FLANK_A_START_POS - (FLANK_A_LEN -1)),
                                       TRUE ~ ERROR_NUMBER)) %>%
    
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    mutate(SEQ_1_WO_A_LEN = nchar(SEQ_1_WO_A))%>%
    
    #if the complete read is flank A (or if there are a few bases at the end that dont match) then make the type WT
    mutate(Type=case_when(FLANK_A_MATCH == SEQ_1 ~ "WT", #if the complete read is flank A
                           FLANK_A_MATCH != SEQ_1 & SEQ_1_WO_A_LEN < 2 ~ "SNV", #or if there are a couple of mismatches
                            TRUE~ "OTHER"))
  
  
   if (EXECUTE_DSB_AREA_CHECK==TRUE){
    

     #############
    
    #Then with exception of the reads thus far called "WT", determine whether the area surrounding the DSB is intact
    #if this is the case, then the read is probably a WT read with seq errors.
    #check in both seq1 and seq2
    #and check allowing for 1 mismatch, though output this as "SNV" because it may be a real 1bp substitution at the DSB.
  data_improved3 = data_improved2b %>%  
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
  
  for (j in row.names(data_improved3)){
    j_int=as.integer(j)
    if (data_improved3$DSB_AREA_COUNT[[j_int]]==1){
      if ((data_improved3$SEQ_1_LEN[[j_int]] >= (data_improved3$DSB_AREA_CHECK[[j_int]]@ranges@start[1]+data_improved3$DSB_AREA_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved3$DSB_AREA_CHECK[[j_int]]@ranges@start[1]>0)){
        data_improved3$DSB_AREA_HIT[[j_int]] = as.character(data_improved3$DSB_AREA_CHECK[[j_int]])
      }else{
        data_improved3$DSB_AREA_HIT[[j_int]] = ""
      }
    }else{
      data_improved3$DSB_AREA_HIT[[j_int]] = ""
    }
    if (data_improved3$DSB_AREA2_COUNT[[j_int]]==1){
      if ((nchar(data_improved3$SEQ_2_first[[j_int]]) >= (data_improved3$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]+data_improved3$DSB_AREA2_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved3$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]>0)){
        data_improved3$DSB_AREA2_HIT[[j_int]] = as.character(data_improved3$DSB_AREA2_CHECK[[j_int]])
      }else{
        data_improved3$DSB_AREA2_HIT[[j_int]] = ""
      }
    }else{
      data_improved3$DSB_AREA2_HIT[[j_int]] = ""
    }
  }
  
  FILE.data_improved3b =   FILE.data_improved3 %>%
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
     FILE.data_improved3b = FILE.data_improved2b %>%
       mutate(
         #the following variables are set to FALSE. This does not mean that there is a DSB area that does not match wt sequence, but rather that it has not been determined.
       DSB_AREA_INTACT_SEQ1 = FALSE,
     DSB_AREA_INTACT_SEQ2 = FALSE,
     DSB_AREA_1MM_SEQ1 = FALSE,
     DSB_AREA_1MM_SEQ2 = FALSE,
     DSB_HIT_MULTI_SEQ1 = FALSE,
     DSB_HIT_MULTI_SEQ2 = FALSE)
     }
  
  
  
    
  FILE.data_improved3c = FILE.data_improved3b %>%
    
    #check whether FLANK_A ends within the T-DNA. (IF LB or RB transguide reaction. if not, limit to the T-DNA)
    mutate(FLANK_A_ENDS_ON_TDNA = case_when(
      
      FOCUS_CONTIG == PLASMID & FOCUS_LOCUS=="LB" & TDNA_IS_LBRB == TRUE & FLANK_A_END_POS >= TDNA_LB_END ~ TRUE,
      FOCUS_CONTIG == PLASMID & FOCUS_LOCUS=="LB" & TDNA_IS_LBRB == FALSE & FLANK_A_END_POS <= TDNA_LB_END ~ TRUE,
      FOCUS_CONTIG == PLASMID & FOCUS_LOCUS=="RB" & TDNA_IS_LBRB == TRUE & FLANK_A_END_POS <= TDNA_RB_END ~ TRUE,
      FOCUS_CONTIG == PLASMID & FOCUS_LOCUS=="RB" & TDNA_IS_LBRB == FALSE & FLANK_A_END_POS >= TDNA_RB_END ~ TRUE,
      
      TRUE ~ FALSE))%>%
    
    #fix end pos for the WT/SNV case
    mutate(FLANK_A_LEN = if_else(Type=="WT"| Type=="SNV" | (FLANK_A_ENDS_ON_TDNA==FALSE & FOCUS_CONTIG == PLASMID),
                                  PRIMER_TO_DSB,
                                  FLANK_A_LEN))%>%
    mutate(FLANK_A_END_POS = if_else(Type=="WT"| Type=="SNV"| (FLANK_A_ENDS_ON_TDNA==FALSE & FOCUS_CONTIG == PLASMID),
                                     FlankAUltEnd,
                                     FLANK_A_END_POS))%>%
    #also determine again what the part of the read is without FLANK_A
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    mutate(SEQ_1_WO_A_LEN = nchar(SEQ_1_WO_A))%>%
    
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        FLANK_A_LEN != ERROR_NUMBER & FLANK_A_LEN != NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
        TRUE ~ ERROR_NUMBER))


  function_time("Step 4 took ")
  
  ###############################################################################
  #Process data: step 5
  ###############################################################################
  
  FILE.data_improved4 = FILE.data_improved3c %>%
    mutate(FLANK_B_CLOSE_BY = case_when(FLANK_A_ISFORWARD == TRUE & TOTAL_REF_STOP > FLANK_B_END_POS &  FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_B_CHROM == FOCUS_CONTIG ~ TRUE,
                                        FLANK_A_ISFORWARD == FALSE & TOTAL_REF_START < FLANK_B_END_POS &  FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_B_CHROM == FOCUS_CONTIG ~ TRUE,
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
                 substr(as.character(eval(parse(text = paste0("genomeseq$`", FLANK_B_CHROM, "`")))), FLANK_B_END_POS-(SEQ_1_LEN-1), FLANK_B_END_POS)
               }else{
                 as.character(reverseComplement(DNAString(substr(as.character(eval(parse(text = paste0("genomeseq$`", FLANK_B_CHROM, "`")))), FLANK_B_END_POS, FLANK_B_END_POS+(SEQ_1_LEN-1)))))
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

                                            TRUE ~ ERROR_NUMBER)) %>%

    
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
      FLANK_B_CHROM==FOCUS_CONTIG & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD == TRUE & FLANK_A_END_POS >= (FLANK_B_START_POS_MH+MH_TD_LEN) ~ TRUE,
      FLANK_B_CHROM==FOCUS_CONTIG & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & FLANK_A_ISFORWARD == FALSE & FLANK_A_END_POS <= (FLANK_B_START_POS_MH-MH_TD_LEN) ~ TRUE,
      
      
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
  
  FILE.data_improved5 = FILE.data_improved4 %>%
    
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
    mutate(FLANK_B_ON_TDNA = case_when(FLANK_B_CHROM == PLASMID & TDNA_IS_LBRB == TRUE & FLANK_B_START_POS >= TDNA_LB_END & FLANK_B_START_POS <= TDNA_RB_END ~ TRUE,
                                       FLANK_B_CHROM == PLASMID & TDNA_IS_LBRB == FALSE & FLANK_B_START_POS <= TDNA_LB_END & FLANK_B_START_POS >= TDNA_RB_END ~ TRUE,
                                       FLANK_B_CHROM == PLASMID_ALT & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_START_POS >= TDNA_ALT_LB_END & FLANK_B_START_POS <= TDNA_ALT_RB_END ~ TRUE,
                                       FLANK_B_CHROM == PLASMID_ALT & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_START_POS <= TDNA_ALT_LB_END & FLANK_B_START_POS >= TDNA_ALT_RB_END ~ TRUE,
                                       TRUE ~ FALSE
                                       ))%>%
    mutate(FLANK_B_DEL = case_when(FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ISFORWARD == TRUE   ~ as.integer(FLANK_B_START_POS - (DSB_FW_END+1)),
                                   FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ISFORWARD == FALSE  ~ as.integer(DSB_FW_END - FLANK_B_START_POS),
                                   
                                   FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - TDNA_LB_END,
                                   FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - TDNA_RB_END,
                                   FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ TDNA_RB_END - FLANK_B_START_POS,
                                   FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ TDNA_LB_END - FLANK_B_START_POS,
                                   
                                   FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - TDNA_ALT_LB_END,
                                   FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ FLANK_B_START_POS - TDNA_ALT_RB_END,
                                   FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ TDNA_ALT_RB_END - FLANK_B_START_POS,
                                   FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ TDNA_ALT_LB_END - FLANK_B_START_POS,
                                   
                                   TRUE ~ ERROR_NUMBER))%>%
    
    #also report what side of the T-DNA flank B is
    mutate(FLANK_B_TDNA_SIDE = case_when(FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         
                                         FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == TRUE & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_ON_TDNA==TRUE ~ "RB",
                                         FLANK_B_CHROM == PLASMID_ALT & FLANK_B_ISFORWARD == FALSE & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_ON_TDNA==TRUE ~ "LB",
                                         
                                         TRUE ~ NA))%>%
    
      
    #calculate total deletion length
    #only for deletions around the DSB with no translocations
    mutate(delSize = case_when(Type=="WT" | Type=="SNV" ~ 0,
                               Type!= "WT" & Type!="SNV" & FLANK_B_CHROM == FOCUS_CONTIG & FLANK_B_ISFORWARD == FLANK_A_ISFORWARD & FLANK_B_ISFORWARD == TRUE & FLANK_A_END_POS < FLANK_B_START_POS & FLANK_B_DEL != ERROR_NUMBER & FLANK_B_CHROM != PLASMID ~ as.integer(FLANK_A_DEL + FLANK_B_DEL),
                               Type!= "WT" & Type!="SNV" & FLANK_B_CHROM == FOCUS_CONTIG & FLANK_B_ISFORWARD == FLANK_A_ISFORWARD & FLANK_B_ISFORWARD == FALSE & FLANK_A_END_POS > FLANK_B_START_POS & FLANK_B_DEL != ERROR_NUMBER & FLANK_B_CHROM != PLASMID ~ as.integer(FLANK_A_DEL + FLANK_B_DEL),
                               TRUE ~ ERROR_NUMBER))%>%

    #correct MH if delSize == 0
    mutate(MH = if_else(delSize <= 0,
                        "",
                        MH)) %>%
    mutate(homologyLength = if_else(insSize == 0,
                                    nchar(MH),
                                    as.integer(-1))) %>% 
    

    #remove reads with FLANK_B 's shorter than the minimum.
    filter(FLANK_B_LEN_MH >= FLANK_B_LEN_MIN)
  
  function_time("Step 6 took ")
  
  ###############################################################################
  #Process data: step 7
  ###############################################################################
  
  FILE.data_improved6 = FILE.data_improved5 %>%
  
    #for SIQPlotter delRelativeStart is FLANK_A_DEL (but negative) and delRelativeEnd is FLANK_B_DEL.
    #delRelativeEndTD should include the tandem duplication.

    mutate(
      delRelativeStart = as.integer(FLANK_A_DEL*-1),
      delRelativeStartTD = as.integer(FLANK_A_DEL*-1),
      delRelativeEnd = if_else(abs(as.integer(FLANK_B_DEL))<=FLANKBEYONDDSB | FLANK_B_CHROM == PLASMID | FLANK_B_CHROM == PLASMID_ALT,
                               as.integer(FLANK_B_DEL),
                               ERROR_NUMBER),
      delRelativeEndTD = if_else( (abs(as.integer(FLANK_B_DEL))<=FLANKBEYONDDSB  | FLANK_B_CHROM == PLASMID | FLANK_B_CHROM == PLASMID_ALT) & FLANK_B_DEL!=ERROR_NUMBER,
                               as.integer(FLANK_B_DEL-tandemDuplicationLength),
                               ERROR_NUMBER)) %>%

  
    
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
      
           FILE.data_improved8pre2 =FILE.data_improved6 %>%
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
        FLANK_B_TDNA_SIDE
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
        .groups="drop"
      )%>%
          mutate(Consensus_freq = 1)
      }else{
        
        FILE.data_improved8pre2 =FILE.data_improved6 %>%
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
            FLANK_B_TDNA_SIDE
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
            Type = as.logical(names(which.max(table(Type)))),
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
        
  FILE.data_improved8 =FILE.data_improved8pre2 %>%
      
      #then calculate the difference between start of flank B (in the read) and the position of of flank B at the end of the mate.
      mutate(ANCHOR_DIST = case_when(FLANK_B_ISFORWARD=TRUE & FLANK_B_CHROM == MATE_FLANK_B_CHROM & MATE_B_END_POS_max > FLANK_B_START_POS ~ as.integer(1+MATE_B_END_POS_max - FLANK_B_START_POS),
                                     FLANK_B_ISFORWARD==FALSE & FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_START_POS > MATE_B_END_POS_min ~ as.integer(1+FLANK_B_START_POS - MATE_B_END_POS_min),
                                     TRUE ~ NF_NUMBER)) %>%
      #adjust the anchor dist for when the anchor is impossibly far away
      mutate(ANCHOR_DIST = if_else(ANCHOR_DIST > MAXANCHORDIST & MATE_FLANK_B_CHROM!=PLASMID & MATE_FLANK_B_CHROM!=PLASMID_ALT,
                                   NF_NUMBER,
                                   ANCHOR_DIST)) %>%

      #determine whether a translocation has occurred
      mutate(Translocation = case_when(Type=="WT" | Type=="SNV" ~ FALSE,
                                       FOCUS_CONTIG ==  FLANK_B_CHROM & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & delRelativeEnd != ERROR_NUMBER & FLANK_A_ISFORWARD==TRUE & FLANK_A_END_POS < FLANK_B_START_POS & FLANK_B_CHROM != PLASMID   ~ FALSE,
                                       FOCUS_CONTIG ==  FLANK_B_CHROM & FLANK_A_ISFORWARD == FLANK_B_ISFORWARD & delRelativeEnd != ERROR_NUMBER & FLANK_A_ISFORWARD==FALSE & FLANK_B_START_POS < FLANK_A_END_POS & FLANK_B_CHROM != PLASMID ~ FALSE,
                                       TRUE ~ TRUE))%>%
      
      #note whether the deletion of the flank B (in case of translocation) has been determined
      mutate(Translocation_del_resolved = if_else(Translocation == TRUE & delRelativeEnd != ERROR_NUMBER,
                                                  TRUE,
                                                  FALSE))%>%
      #if the deletion size cannot be determined it becomes 0 for plotting purposes
      mutate(delRelativeEnd = if_else(delRelativeEnd == ERROR_NUMBER,
                                      0,
                                      delRelativeEnd),
             delRelativeEndTD = if_else(delRelativeEndTD == ERROR_NUMBER,
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
        

        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize == 0 & insSize != 0 & tandemDuplicationLength < TD_SIZE_CUTOFF ~ "INSERTION",
        Type != "WT" & Type != "SNV" & FAKE_DELIN_CHECK == FALSE & delSize == 0 & insSize != 0 & tandemDuplicationLength >= TD_SIZE_CUTOFF ~ "TANDEMDUPLICATION",
        
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
      TRIM_LEN
    ) 

   

  function_time("Step 7 took ")
  
  ###############################################################################
  #Process data: step 8
  ###############################################################################
  
  #calculate the fraction of reads with a certain outcome within a library
  FILE.data_improved9 = FILE.data_improved8 %>% group_by(FILE_NAME) %>% summarize(ReadCountTotal =
                                                                                  sum(ReadCount),
                                                                        .groups="drop")
  function_time("Step 8 took ")
  
  ###############################################################################
  #Process data: step 9
  ###############################################################################
  
  FILE.data_improved10 = left_join(FILE.data_improved8, FILE.data_improved9, by = "FILE_NAME") %>%
    mutate(fraction = ReadCount / ReadCountTotal) %>%
    #add columns for all the input options and software version
    mutate(countReadsTotal = NULL,
           FlankAUltEnd = FlankAUltEnd,
           Flank_A_isforward = FLANK_A_ISFORWARD,
           Focus_contig = FOCUS_CONTIG,
           Genotype = Genotype,
           DNASample = DNASample,
           Ecotype = Ecotype,
           RunID = RunID,
           Plasmid = PLASMID,
           Plasmid_alt = PLASMID_ALT,
           AgroGeno = AgroGeno,
           RemoveNonTranslocation = REMOVENONTRANS,
           GroupSamePosition = GROUPSAMEPOS,
           Primer_match_perfect = Primer_match_perfect,
           Alias = paste0(Library, "_", RunID),
           DSB_FW_END = DSB_FW_END,
           DSB_CONTIG = DSB_CONTIG,
           TDNA_LB_END = TDNA_LB_END,
           TDNA_RB_END = TDNA_RB_END,
           TDNA_IS_LBRB = TDNA_IS_LBRB,
           TDNA_ALT_LB_END = TDNA_ALT_LB_END,
           TDNA_ALT_RB_END = TDNA_ALT_RB_END,
           TDNA_ALT_IS_LBRB = TDNA_ALT_IS_LBRB,
           FLANKBEYONDDSB = FLANKBEYONDDSB,
           MinumumReadLength = MINLEN,
           MAXANCHORDIST = MAXANCHORDIST,
           program_version = hash,
           TD_SIZE_CUTOFF = TD_SIZE_CUTOFF,
           EXECUTE_DSB_AREA_CHECK = EXECUTE_DSB_AREA_CHECK)
  

  
  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, FILE.data_improved10)
  saveWorkbook(work_book, file = paste0(output_dir, Sample, "_", RunID, "_CISTRANSGUIDE_V2.xlsx"), overwrite = TRUE)
  
  
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
    filter(AnchorCount >= ANCHORCUTOFF,
           ANCHOR_DIST >= MINANCHORDIST,
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
    group_by(Alias, FLANK_B_CHROM, Plasmid, FLANK_B_ISFORWARD, DNASample, Subject, ID, Focus_contig, Genotype, Ecotype, Plasmid_alt, Family, FlankAUltEnd, AgroGeno, RemoveNonTranslocation, GroupSamePosition, Translocation, Translocation_del_resolved, TANDEM_DUPLICATION, Primer_match_perfect, DSB_FW_END, DSB_OVERHANG, DSB_CONTIG, TDNA_LB_END, TDNA_RB_END, TDNA_IS_LBRB, TDNA_ALT_LB_END, TDNA_ALT_RB_END, TDNA_ALT_IS_LBRB, FLANK_B_TDNA_SIDE)%>%
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
           ReadCount = ReadCountSum)

  #get a list of families
  GLOBAL.wb_family = GLOBAL.sample_info2 %>% select(Family) %>% distinct() %>% filter(Family!=0)
  if (nrow(GLOBAL.wb_family) == 0) {
    funlog("no family info detected")
    #if no families are indicated

    GLOBAL.wb_flag = GLOBAL.total_data_near_positioncombined %>%
        #then examine positions across samples and remove those that occur multiple times
        group_by(FLANK_B_START_POS) %>%
        mutate(duplicate_position = if_else(n() > 1,
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
        
        GLOBAL.wb_filter_current = total_data_near_positioncombined %>%
          filter(Family != i | Alias == j) %>% #events are either not of the current family, or they belong to the current alias
          group_by(FLANK_B_START_POS) %>%
          mutate(duplicate_position = if_else(n() > 1,
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
      mutate(duplicate_position = if_else(n() > 1,
                                          TRUE,
                                          FALSE)) %>%
      
      ungroup() %>%
      filter(duplicate_position == FALSE)
    funlog("combining surviving family and nonfamily events")
    wb_flag = rbind(wb_filter_total, wb_nonfamily)  
    
  }
} else{
  funlog("flagging problems only")
  GLOBAL.wb_flag = GLOBAL.wb %>%
    group_by(FLANK_B_START_POS) %>%
    mutate(duplicate_position = if_else(n() > 1,
                                        TRUE,
                                        FALSE)) %>%
    
    ungroup() 
}

GLOBAL.wb_flag1 = GLOBAL.wb_flag %>%
  mutate(RemoveProblematicEvents = REMOVEPROBLEMS)%>%
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
