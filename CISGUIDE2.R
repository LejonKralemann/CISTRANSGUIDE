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
input_dir= "./input/"
output_dir= "./output/"
MAX_DIST_FLANK_B_END = 10000 #distance from end of flank B to DSB, determines max deletion size and also affects maximum insertion size
LOCUS_WINDOW = 1000 #size of the window centered on the DSB, RB nick, or LB nick to determine locus info
GROUPSAMEPOS=TRUE #if true, it combines reads with the same genomic pos, which helps in removing artefacts. Typically used for TRANSGUIDE, but disabled for CISGUIDE.
REMOVENONTRANS=TRUE #if true, it only considers translocations. Typically used for TRANSGUIDE, but disabled for CISGUIDE. Note that some translocations on the same chromosome will also be removed thusly.
REMOVEPROBLEMS=TRUE #if true it removes all problematic reads from the combined datafile. Note if this is false, no duplicate filtering will be performed, because first reads due to barcode hopping need to be removed by removing events with few anchors.
CONVERTWT=TRUE #if true, it sets delRelativeStart, delRelativeEnd, insSize, delSize, homologyLength all to 0, and Translocation to FALSE when Type is WT. For troubleshooting this option should be set to FALSE. For SIQplotteR it should be TRUE.
ANCHORCUTOFF=3 #each event needs to have at least this number of anchors, otherwise it is marked as problematic (and potentially removed) 
MINANCHORDIST=150 #should be matching a situation where the mate is 100% flank B.
MAXTARGETLENGTH=10000 #this limits deletion calculation for when flank B is on the target chrom, but far away
MAXANCHORDIST=5000

###############################################################################
#set parameters - non-adjustable
###############################################################################
NF_NUMBER = as.integer(-99999999) #don't change
ERROR_NUMBER = as.integer(99999999) #don't change
hash=system("git rev-parse HEAD", intern=TRUE)
hash_little=substr(hash, 1, 8)
sample_info = read.csv(paste0(input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
TIME_START=round(as.numeric(Sys.time())*1000, digits=0)
FLANK_B_LEN_MIN = 30 #minimum length of flank B. Also determines the size of DSB_AREA_SEQ. do not change because the preprocessing program will still be set at 30.
MINMAPQUALA = 42 #minimum mapping quality (phred). 42 means a perfect, unambiguous match (well, should be)
MINMAPQUALB = 42 #same, but for flank B
PercentageDone = 0 #var for indicating progress
TotalFileSize = 0
CurrentFileSize = 0

###############################################################################
#Initial checks
###############################################################################
if (file.exists(paste0(input_dir, "Sample_information.txt"))==FALSE){
  message("Sample information sheet not present, aborting")
  quit()}
if (file.exists(paste0(input_dir, "read_numbers.txt"))==FALSE){
  message("Read numbers file not found, aborting")
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
function_time <-function(text){
  TIME_CURRENT=round(as.numeric(Sys.time())*1000, digits=0)
  message(paste0(text, (TIME_CURRENT - TIME_START), " milliseconds"))
  TIME_START<<-round(as.numeric(Sys.time())*1000, digits=0)
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###############################################################################
#Process data: step 0
#calculating total work
###############################################################################

message("calculating total work")
for (i in row.names(sample_info)){
  Sample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  RunID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(RunID))
  
  if (file.exists(paste0(input_dir, Sample, "_", RunID, "_A.txt"))==FALSE){
    next
  }else{
    TotalFileSize = TotalFileSize + (file.info((paste0(input_dir, Sample, "_", RunID, "_A.txt"))))$size
  }
}

message(paste0("Total file size to process: ", TotalFileSize, " bytes"))


###############################################################################
#Process data: step 1
#checking the file and reading metadata
###############################################################################

message("Checking file and reading metadata")

#check whether input file exists
for (i in row.names(sample_info)){
  
  ####################  general variables acquired from the information sheet  #####################
  
  Sample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  RunID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(RunID))
  
  ####################  check for existence of files  #####################
  
  if (file.exists(paste0(input_dir, Sample, "_", RunID, "_A.txt"))==FALSE){
    message("Primary processed file not found, moving to the next sample")
    next
  }else if (file.exists(paste0(output_dir, Sample, "_", RunID, "_CISTRANSGUIDE_V2.xlsx"))==TRUE){   #check whether file has already been processed
    message(paste0("File ", output_dir, Sample, "_", RunID, "_A.txt has already been processed, moving to the next sample"))
    #show progress
    CurrentFileSize = (file.info((paste0(input_dir, Sample, "_", RunID, "_A.txt"))))$size
    PercentageDone = PercentageDone + ((CurrentFileSize/TotalFileSize)*99)
    message(paste0("CISTRANSGUIDE analysis ", round(PercentageDone, digits=3), "% complete"))
    next
    }else{
      message(paste0("Processing ",input_dir, Sample, "_", RunID, "_A.txt"))
      CurrentFileSize = (file.info((paste0(input_dir, Sample, "_", RunID, "_A.txt"))))$size
    }
  
  data = read.csv(paste0(input_dir, Sample, "_", RunID, "_A.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
  if (nrow(data)==0){
    message("Primary processed file empty, moving to the next sample")
    next
  }
  
  ####################  continue general variables acquired from the information sheet  #####################
  
  DSB_CONTIG = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_CONTIG))#chromosome name or NA
  FOCUS_LOCUS = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Locus_name))#LB or RB if TRANSGUIDE, or a name of a locus if CISGUIDE
  Genotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Genotype))
  PLASMID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid))
  PLASMID_ALT = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid_alt))
  REF = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ref))
  DNASample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DNA))
  Ecotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ecotype))
  Library = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  AgroGeno = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(AgroGeno))
  FLANK_A_ORIENT = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FLANK_A_ORIENT))
  Primer_seq = str_replace_all(toupper(as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Primer))), "TCAGACGTGTGCTCTTCCGATCT", "")
  DSB_FW_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_FW_END)) #end of left flank before DSB
  DSB_OVERHANG = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_OVERHANG)) #e.g. 1 if cas9, 5 if cas12a
  TDNA_LB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_LB_END))
  TDNA_RB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_RB_END))
  TDNA_ALT_LB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_ALT_LB_END))
  TDNA_ALT_RB_END = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(TDNA_ALT_RB_END))
 

  ####################  calculated general variables  #####################
  
  if (TDNA_LB_END < TDNA_RB_END){
    TDNA_IS_LBRB = TRUE
  }else{
    TDNA_IS_LBRB = FALSE
  }
  if (is.na(TDNA_ALT_LB_END)==FALSE & is.na(TDNA_ALT_RB_END)==FALSE){
  if (TDNA_ALT_LB_END < TDNA_ALT_RB_END){
    TDNA_ALT_IS_LBRB = TRUE
  }else{
    TDNA_ALT_IS_LBRB = FALSE
  }
  }else{
    TDNA_ALT_IS_LBRB = NA
  }
  
  if (FOCUS_LOCUS == "LB"){
    FlankAUltEnd = TDNA_LB_END
    FOCUS_CONTIG = PLASMID
    if (FLANK_A_ORIENT == "FW"){
      if (TDNA_IS_LBRB == TRUE){
        message("T-DNA orientation conflict. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
        next 
      }
      FlankBUltStart = TDNA_LB_END + 1
    }else{
      if (TDNA_IS_LBRB == FALSE){
        message("T-DNA orientation conflict. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
        next 
      }
      FlankBUltStart = TDNA_LB_END - 1
    }
  }else if (FOCUS_LOCUS == "RB"){
    FlankAUltEnd = TDNA_RB_END  
    FOCUS_CONTIG = PLASMID
    if (FLANK_A_ORIENT == "FW"){
      if (TDNA_IS_LBRB == FALSE){
        message("T-DNA orientation conflict. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
        next 
      }
      FlankBUltStart = TDNA_RB_END + 1
    }else{
      if (TDNA_IS_LBRB == TRUE){
        message("T-DNA orientation conflict. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
        next 
      }
      FlankBUltStart = TDNA_RB_END - 1
    }
  }else {
    FOCUS_CONTIG = DSB_CONTIG
    if (FLANK_A_ORIENT == "FW"){
      FlankAUltEnd = DSB_FW_END
      FlankBUltStart = (DSB_FW_END+1)-DSB_OVERHANG
    }else{
      FlankAUltEnd = (DSB_FW_END+1)-DSB_OVERHANG
      FlankBUltStart = DSB_FW_END
    }
    }
 
  Primer_seq_len = nchar(Primer_seq)
  
  ####################  REF checks  #####################
  
  if (file.exists(paste0(input_dir, REF))==FALSE){
    message("Reference fasta not found. Moving to next sample.")
    next
  }else{
    message(paste0("Using ref ", input_dir, REF))
  }
  
  genomeseq = readDNAStringSet(paste0(input_dir, REF) , format="fasta")
  contig_seq = as.character(eval(parse(text = paste0("genomeseq$", FOCUS_CONTIG))))
  if (length(contig_seq)==0){
    message("Focus contig not found in reference fasta. Did you fill in the Sample_information sheet correctly? Moving to next sample.")
    next 
  }
  
  ####################  continue calculated general variables  #####################
  
  Primer_match = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  
  DSB_AREA_SEQ = (if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= FlankAUltEnd - ((FLANK_B_LEN_MIN/2)-1), stop= FlankAUltEnd + (FLANK_B_LEN_MIN/2))
  }else if (FLANK_A_ORIENT == "RV"){
    as.character(reverseComplement(DNAString(substr(contig_seq, start= FlankAUltEnd - (FLANK_B_LEN_MIN/2), stop= FlankAUltEnd + ((FLANK_B_LEN_MIN/2)-1)))))
  }else{
    ""
  })
  
  DSB_AREA_SEQ_RC = as.character(reverseComplement(DNAString(DSB_AREA_SEQ)))
  
  ####################  Primer checks  #####################
  
  if (Primer_seq != "NA"){
  FASTA_MODE = FALSE
  }else{
  FASTA_MODE = TRUE
  message("Primer seq not found, running in fasta mode.")}
  
  if (FASTA_MODE == FALSE){
  if (FLANK_A_ORIENT == "FW"){
    Primer_match = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))  
    Primer_match_3 = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 3, fixed=TRUE)) 
    if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
      Primer_pos = as.numeric(Primer_match$start)
      Primer_match_perfect=TRUE
    }else if (nrow(Primer_match) >1){
      message("Primer found several times in the genome")
      next
    }else if (nrow(Primer_match_3) == 1){
      message("Note! Primer does not match fully. Continuing anyway.")
      Primer_pos = as.numeric(Primer_match_3$start)
      Primer_match_perfect=FALSE
    }else{
      message("Primer not found")
      next
    }
  }else if (FLANK_A_ORIENT == "RV"){
    Primer_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
    Primer_match_3 = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 3, fixed=TRUE))
    if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
      Primer_pos = as.numeric(Primer_match$end)
      Primer_match_perfect=TRUE
    }else if (nrow(Primer_match) >1){
      message("Primer found several times in the genome")
      next
    }else if (nrow(Primer_match_3) == 1){
      message("Note! Primer does not match fully. Continuing anyway.")
      Primer_pos = as.numeric(Primer_match_3$end)
      Primer_match_perfect=FALSE
    }else{
      message("Primer not found")
      next
    }
  }else{
    next
  }
    
  ####################  continue calculated general variables  #####################  
    
  #calculate the length from primer to DSB
  PRIMER_TO_DSB_GLOBAL = if (FLANK_A_ORIENT == "FW"){
    FlankAUltEnd - (Primer_pos -1)
  }else if (FLANK_A_ORIENT=="RV"){
    Primer_pos - (FlankAUltEnd -1)
  }else{
    ERROR_NUMBER
  }

  if (PRIMER_TO_DSB_GLOBAL>300){
    message("primer too far away from the DSB. Did you fill in the Sample_information sheet correctly?")
    next
  }

  #get the REF seq for flank A. from primer start to DSB +3 if RV primer, not if FW primer. Because CAS9 can cut further away from the PAM, but not closer. So the FLANK_A_REF is going as far as FLANK A is allowed to go.

  FLANK_A_REF_GLOBAL = if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= Primer_pos, stop= FlankAUltEnd)
  }else if (FLANK_A_ORIENT=="RV"){
    as.character(reverseComplement(DNAString(
      substr(contig_seq, start= FlankAUltEnd, stop= Primer_pos))
    ))
  }else{
    ""
  }

  GLOBAL_TOTAL_REF = if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= Primer_pos, stop= Primer_pos+MAX_DIST_FLANK_B_END+PRIMER_TO_DSB_GLOBAL)
  }else if (FLANK_A_ORIENT=="RV"){
    as.character(reverseComplement(DNAString(
      substr(contig_seq, start= Primer_pos-(MAX_DIST_FLANK_B_END+PRIMER_TO_DSB_GLOBAL), stop= Primer_pos))
    ))
  }else{
    ""
  }}
  
  #if there was a problem with the primer (some mismatches) then the refs need to be changed so they match the reads.
  #also the minimum mapping quality will need to be adjusted because of the mismatches
  if (Primer_match_perfect==FALSE){
    FLANK_A_REF_GLOBAL = paste0(Primer_seq, substr(FLANK_A_REF_GLOBAL, start=Primer_seq_len+1, stop=nchar(FLANK_A_REF_GLOBAL)))
    GLOBAL_TOTAL_REF = paste0(Primer_seq, substr(GLOBAL_TOTAL_REF, start=Primer_seq_len+1, stop=nchar(GLOBAL_TOTAL_REF)))
    MINMAPQUALA = 8
  }
  
  #set the minimum length of a read
  if (FASTA_MODE == FALSE){
  MINLEN = PRIMER_TO_DSB_GLOBAL+FLANK_B_LEN_MIN
  }else{
  MINLEN = 60
  message("MINLEN set to 60 because no primer seq available")
  }
  
  function_time("Step 1 took ")

  ###############################################################################
  #Process data: step 2
  ###############################################################################
  
  data_improved_a  = data %>%
    
    #filter(QNAME == "A01685:194:HHCVJDSX7:3:2532:20690:21449")%>%
    
    #Count number of Ns and remove any reads with Ns
    mutate(NrN = str_count(SEQ_1, pattern = "N"),
           SEQ_1_LEN = nchar(SEQ_1)) %>%
    filter(NrN < 1) %>%
    
    #remove reads that do not have perfectly mapped ends
    filter(A_MAPQ >= MINMAPQUALA,
           B_MAPQ >= MINMAPQUALB)%>%
    
    #filter away reads that are too short
    filter(!(SEQ_1_LEN < MINLEN)) %>%
    
    #calculate the average base quality
    rowwise() %>%
    mutate(AvgBaseQual_1 = mean(utf8ToInt(QUAL_1)-33)) %>%
    mutate(AvgBaseQual_2 = mean(utf8ToInt(QUAL_2)-33)) %>%
    ungroup()
  
  if (REMOVENONTRANS==TRUE){
    data_improved_b = data_improved_a %>%
      filter(FLANK_B_CHROM != FOCUS_CONTIG)
  }else{
    data_improved_b = data_improved_a
  }
  
  
  data_improved_c = data_improved_b %>%
    
    #calculate the end position of the B flank in the mate. This is the end the furthest away from the junction.
    mutate(MATE_B_END_POS = case_when(MATE_B_ORIENT=="FW" ~ as.integer(MATE_B_POS),
                                      MATE_B_ORIENT=="RV" ~ as.integer(MATE_B_POS + 29),
                                      TRUE ~ ERROR_NUMBER))  %>%
    
    #test whether B flanks from read and mate are mapped to same chromosome, and agree on orientation.
    #note that the two reads are in two different orientations, so the B orientation should not agree
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_ORIENT != MATE_B_ORIENT,
                                              TRUE,
                                              FALSE)) 
  
  data_improved1 = data_improved_c %>%
    
    #count reads, taking the highest quals
    #acts as basic dupfilter
    #and also summarizes anchors, but puts them in a list
    group_by(
      A_CHROM,
      A_POS,
      FLANK_B_CHROM,
      B_POS,
      FLANK_B_ORIENT,
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
    mutate(SEQ_1_start = substr(SEQ_1, 1, Primer_seq_len)) %>%
    filter(SEQ_1_start == Primer_seq) %>%
    
    #the lines below until the grouping are under investigation. Maybe doesn't work well yet. First make the code so that the output is identical to the output before the changes of a non fasta sample. Then test the fasta data.
    mutate(PRIMER_SEQ = if (FASTA_MODE == FALSE){
      PRIMER_SEQ
      }else{
        substr(SEQ_1, 1, 30)
      }) %>%
    mutate(PRIMER_POS_FAKE_match = if (FASTA_MODE == FALSE){
      NA
    }else{
      if (FLANK_A_ORIENT == "FW"){
        list(matchPattern(DNAString(PRIMER_SEQ), DNAString(contig_seq), max.mismatch = 0))
      }else{
        list(matchPattern(as.character(reverseComplement(DNAString(PRIMER_SEQ))), DNAString(contig_seq), max.mismatch = 0))
      }}) %>%
    mutate(PRIMER_POS_FAKE = if (FASTA_MODE == FALSE){
      as.integer(ERROR_NUMBER)
    }else{
      if (FLANK_A_ORIENT == "FW"){
      as.integer(PRIMER_POS_FAKE_match@ranges@start)
      }else{
      as.integer(PRIMER_POS_FAKE_match@ranges@start)+29
      }
    }) %>%
    mutate(PRIMER_TO_DSB = if (FASTA_MODE == FALSE){
      as.integer(PRIMER_TO_DSB_GLOBAL)
    }else{
      if (FLANK_A_ORIENT == "FW"){
        FlankAUltEnd - (PRIMER_POS_FAKE -1)
      }else{
        PRIMER_POS_FAKE - (FlankAUltEnd -1)
      }}) %>%
    mutate(FLANK_A_REF = if (FASTA_MODE == FALSE){
      FLANK_A_REF_GLOBAL
    }else{
      if (FLANK_A_ORIENT == "FW"){
        substr(contig_seq, start= PRIMER_POS_FAKE, stop= FlankAUltEnd)
      }else{
        as.character(reverseComplement(DNAString(substr(contig_seq, start= FlankAUltEnd, stop= PRIMER_POS_FAKE))))
      }}) %>%
    mutate(TOTAL_REF = if (FASTA_MODE == FALSE){
      GLOBAL_TOTAL_REF
    }else{
      if (FLANK_A_ORIENT == "FW"){
        substr(contig_seq, start= PRIMER_POS_FAKE, stop= PRIMER_POS_FAKE+MAX_DIST_FLANK_B_END+PRIMER_TO_DSB)
      }else{
        as.character(reverseComplement(DNAString(
          substr(contig_seq, start= PRIMER_POS_FAKE-(MAX_DIST_FLANK_B_END+PRIMER_TO_DSB), stop= PRIMER_POS_FAKE))))
      }}) %>%
    ungroup() 
  
  #check if any reads have survived
  if (nrow(data_improved1)==0){
    message(paste0("No reads surviving for sample ", DNASample))
    #show progress
    PercentageDone = PercentageDone + ((CurrentFileSize/TotalFileSize)*100)
    message(paste0("CISTRANSGUIDE analysis ", round(PercentageDone, digits=3), "% complete"))
    next
  }
  
  function_time("Step 2 took ")
  
  ###############################################################################
  #Process data: step 3
  ###############################################################################
  
  data_improved2 = data_improved1 %>%
    
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
    mutate(FLANK_B_END_POS = case_when(FLANK_B_ORIENT=="RV" ~ as.integer(B_POS),
                                       FLANK_B_ORIENT=="FW" ~ as.integer(B_POS + 29),
                                       TRUE ~ ERROR_NUMBER))
    

  function_time("Step 3 took ")
  
  ###############################################################################
  #Process data: step 4
  ###############################################################################
  
  data_improved3 = data_improved2 %>%

    #Find how much SEQ_1 matches with FLANK_A_REF. Allow 1 bp mismatch somewhere, if the alignment after that continues for at least another 10 bp.
    mutate(FLANK_A_REF_LEN = as.integer(nchar(FLANK_A_REF))) %>%
    rowwise() %>%
    mutate(FLANK_A_MATCH = matcher_skipper(FLANK_A_REF, SEQ_1)) %>%
    ungroup() %>%
    mutate(FLANK_A_LEN = nchar(FLANK_A_MATCH)) %>%
    mutate(FLANK_A_END_POS = case_when(FLANK_A_ORIENT == "FW" ~ as.integer(FLANK_A_START_POS + (FLANK_A_LEN -1)),
                                       FLANK_A_ORIENT == "RV" ~ as.integer(FLANK_A_START_POS - (FLANK_A_LEN -1)),
                                       TRUE ~ ERROR_NUMBER)) %>%
    
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    mutate(SEQ_1_WO_A_LEN = nchar(SEQ_1_WO_A))%>%
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        FLANK_A_LEN != ERROR_NUMBER & FLANK_A_LEN != NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
        TRUE ~ ERROR_NUMBER
      )
    ) 
  
  function_time("Step 4 took ")
  
  ###############################################################################
  #Process data: step 5
  ###############################################################################
  
  data_improved3b = data_improved3 %>%
    mutate(FLANK_B_CLOSE_BY = case_when(abs(FlankAUltEnd - FLANK_B_END_POS) < MAX_DIST_FLANK_B_END ~ TRUE,
                                        TRUE ~ FALSE)) %>%
    #FLANK_B_REF. This ref includes homology.
    rowwise() %>%
    mutate(FLANK_B_REF =
             if (FLANK_B_CHROM == FOCUS_CONTIG & 
                       FLANK_B_CLOSE_BY==TRUE & 
                       FLANK_A_ORIENT == FLANK_B_ORIENT){
               if (FLANK_B_ORIENT == "FW" & 
                   FLANK_A_START_POS < (FLANK_B_END_POS-(SEQ_1_LEN-1)) & 
                   FLANK_A_END_POS < FLANK_B_END_POS) {
                 substr(TOTAL_REF, start = 1, stop = (FLANK_B_END_POS-(FLANK_A_START_POS-1)))
               }else if (FLANK_B_ORIENT == "RV" & 
                         FLANK_A_START_POS > (FLANK_B_END_POS+(SEQ_1_LEN-1)) & 
                         FLANK_A_END_POS > FLANK_B_END_POS){
                 substr(TOTAL_REF, start = 1, stop = (FLANK_A_START_POS-(FLANK_B_END_POS-1)))
               }else{
                 if (FLANK_B_ORIENT == "FW"){
                 substr(as.character(eval(parse(text = paste0("genomeseq$", FLANK_B_CHROM)))), FLANK_B_END_POS-(SEQ_1_LEN-1), FLANK_B_END_POS)
                 }else{
                   as.character(reverseComplement(DNAString(substr(as.character(eval(parse(text = paste0("genomeseq$", FLANK_B_CHROM)))), FLANK_B_END_POS, FLANK_B_END_POS+(SEQ_1_LEN-1)))))
                 }
               }
             }else{
               if (FLANK_B_ORIENT == "FW"){
                 substr(as.character(eval(parse(text = paste0("genomeseq$", FLANK_B_CHROM)))), FLANK_B_END_POS-(SEQ_1_LEN-1), FLANK_B_END_POS)
               }else{
                 as.character(reverseComplement(DNAString(substr(as.character(eval(parse(text = paste0("genomeseq$", FLANK_B_CHROM)))), FLANK_B_END_POS, FLANK_B_END_POS+(SEQ_1_LEN-1)))))
               }
             }) %>%
    mutate(FLANK_B_MATCH = stri_reverse(matcher_skipper(stri_reverse(FLANK_B_REF), stri_reverse(SEQ_1)))) %>%     #find the full flank b match, while skipping over seq errors. But causes a problem on focus contig when del=1 en ins=1
    ungroup() %>%
    mutate(FLANK_B_MATCH_LEN = nchar(FLANK_B_MATCH)) %>%
  
  #then do some checks to find wt events that because of mutations did not get called as such
  #first find a sequence around the DSB that would indicate no DSB has been made or is repaired perfectly
  #check in both seq1 and seq2
  rowwise() %>%
    mutate(DSB_AREA_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1))) %>%
    mutate(DSB_AREA_COUNT = length(DSB_AREA_CHECK@ranges))%>%
    mutate(DSB_AREA2_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ_RC), DNAString(SEQ_2_first), max.mismatch = 1))) %>%
    mutate(DSB_AREA2_COUNT = length(DSB_AREA2_CHECK@ranges))%>%
    ungroup()%>%
    mutate(DSB_AREA_HIT = "",
           DSB_AREA2_HIT = "",
           DSB_AREA_INTACT = FALSE,
           DSB_AREA_1MM = FALSE)
  
  for (j in row.names(data_improved3b)){
    j_int=as.integer(j)
    if (data_improved3b$DSB_AREA_COUNT[[j_int]]==1){
      if ((data_improved3b$SEQ_1_LEN[[j_int]] >= (data_improved3b$DSB_AREA_CHECK[[j_int]]@ranges@start[1]+data_improved3b$DSB_AREA_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved3b$DSB_AREA_CHECK[[j_int]]@ranges@start[1]>0)){
        data_improved3b$DSB_AREA_HIT[[j_int]] = as.character(data_improved3b$DSB_AREA_CHECK[[j_int]])
      }else{
        data_improved3b$DSB_AREA_HIT[[j_int]] = ""
      }
    }else{
      data_improved3b$DSB_AREA_HIT[[j_int]] = ""
    }
    if (data_improved3b$DSB_AREA2_COUNT[[j_int]]==1){
      if ((nchar(data_improved3b$SEQ_2_first[[j_int]]) >= (data_improved3b$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]+data_improved3b$DSB_AREA2_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved3b$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]>0)){
        data_improved3b$DSB_AREA2_HIT[[j_int]] = as.character(data_improved3b$DSB_AREA2_CHECK[[j_int]])
      }else{
        data_improved3b$DSB_AREA2_HIT[[j_int]] = ""
      }
    }else{
      data_improved3b$DSB_AREA2_HIT[[j_int]] = ""
    }
  }
  
  
  
  
  data_improved4 =   data_improved3b %>%
    mutate(DSB_AREA_INTACT = if_else(DSB_AREA_HIT == DSB_AREA_SEQ | DSB_AREA2_HIT == DSB_AREA_SEQ_RC,
                                     "TRUE",
                                     "FALSE"))%>%
    mutate(DSB_AREA_1MM = if_else(
      DSB_AREA_INTACT == FALSE & (DSB_AREA_COUNT>0 | DSB_AREA2_COUNT>0),
      "TRUE",
      "FALSE")) %>%
    mutate(CASE_WT = if_else((
      DSB_AREA_INTACT==TRUE ),
      TRUE,
      FALSE)) %>%
    mutate(DSB_HIT_MULTI = if_else(DSB_AREA_COUNT>1 | DSB_AREA2_COUNT>1,
                                   "TRUE",
                                   "FALSE")) %>%
    rowwise()%>%
    #fix the FLANK_B_match in case focus contig has 1bp del and 1bp ins
    mutate(FLANK_B_MATCH_LEN = if_else(
      FLANK_B_MATCH_LEN == SEQ_1_LEN & DSB_AREA_INTACT == FALSE,
      lcprefix(stri_reverse(FLANK_B_REF), stri_reverse(SEQ_1)),
               FLANK_B_MATCH_LEN))%>%
    ungroup() %>%
   
    #flank b start position including MH
    mutate(FLANK_B_START_POS_MH = case_when(FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_END_POS-(FLANK_B_MATCH_LEN-1)),
                                            FLANK_B_ORIENT == "RV" ~ as.integer(FLANK_B_END_POS+(FLANK_B_MATCH_LEN-1)),
                                            TRUE ~ ERROR_NUMBER)) %>%

    rowwise() %>%
    mutate(FLANK_B_LEN_MH = if_else(FLANK_B_ORIENT == "FW",
                                    FLANK_B_END_POS-(FLANK_B_START_POS_MH-1),
                                    FLANK_B_START_POS_MH-(FLANK_B_END_POS-1))) %>%
    #Extract the MH sequence
    mutate(MH = 
      if (FLANK_B_CHROM == FOCUS_CONTIG & FLANK_B_CLOSE_BY==TRUE & FLANK_A_ORIENT == FLANK_B_ORIENT & (FLANK_A_LEN > (SEQ_1_LEN - FLANK_B_LEN_MH))) {
        if (FLANK_B_ORIENT == "FW") {
          if (FLANK_A_END_POS >= FLANK_B_START_POS_MH) {
            ""
          } else{
            substr(SEQ_1, start = (1 + SEQ_1_LEN - FLANK_B_LEN_MH), stop = FLANK_A_LEN)
          }
        } else{
          if (FLANK_A_END_POS <= FLANK_B_START_POS_MH) {
            ""
          } else{
            substr(SEQ_1, start = (1 + SEQ_1_LEN - FLANK_B_LEN_MH), stop = FLANK_A_LEN)
          }
        }
      } else{
        substr(SEQ_1, start = (1 + SEQ_1_LEN - FLANK_B_LEN_MH), stop = FLANK_A_LEN)
      }
    ) %>%   
    #flank B start pos excluding MH (for del calculation)
      mutate(FLANK_B_START_POS_DEL = 
        if (FLANK_B_ORIENT == "FW"){
          if (FLANK_B_CHROM == FOCUS_CONTIG & FLANK_A_END_POS >= (FLANK_B_START_POS_MH + nchar(MH)) & FLANK_A_END_POS < FLANK_B_END_POS){
            FLANK_A_END_POS +1 #if filler is b continuation (filler plus homology)
          }else{
            FLANK_B_START_POS_MH + nchar(MH)
          }
          
        }else{#if RV
          if (FLANK_B_CHROM == FOCUS_CONTIG & FLANK_A_END_POS <= (FLANK_B_START_POS_MH - nchar(MH)) & FLANK_A_END_POS > FLANK_B_END_POS){
            FLANK_A_END_POS -1 #if filler is b continuation (filler plus homology)
          }else{
          FLANK_B_START_POS_MH - nchar(MH)}
        })%>%
    #correct the flank b start pos in case there is a deletion in A, and an insertion, and flank b continues beyond the FLANK_B_ULTSTART
    mutate(FLANK_B_START_POS_DEL = case_when(FLANK_B_CHROM == FOCUS_CONTIG & FLANK_B_ORIENT == "FW" & FLANK_B_START_POS_DEL < FlankBUltStart & FLANK_B_END_POS > FlankBUltStart ~ FlankBUltStart,
                                             FLANK_B_CHROM == FOCUS_CONTIG & FLANK_B_ORIENT == "RV" & FLANK_B_START_POS_DEL > FlankBUltStart & FLANK_B_END_POS < FlankBUltStart ~ FlankBUltStart,
                                             TRUE ~ FLANK_B_START_POS_DEL)) %>%
    
      #calculate length of flank B minus the MH
    mutate(FLANK_B_LEN_DEL = if_else(FLANK_B_ORIENT == "FW",
                                     FLANK_B_END_POS-(FLANK_B_START_POS_DEL-1),
                                     FLANK_B_START_POS_DEL-(FLANK_B_END_POS-1))) %>%

    ungroup()
  

  function_time("Step 5 took ")
  
  ###############################################################################
  #Process data: step 6
  ###############################################################################
  
  data_improved5 = data_improved4 %>%
    
    #determine the filler sequence 
    mutate(
      FILLER = if_else(
        FLANK_B_LEN_DEL < SEQ_1_WO_A_LEN,
        substr(SEQ_1_WO_A, start = 1, stop = SEQ_1_WO_A_LEN - FLANK_B_LEN_DEL),	 
        "" #no filler
      )) %>%
    mutate(insSize = nchar(FILLER)) %>%
    

    #determine how much from flank B has been deleted
    mutate(FLANK_B_DEL = case_when(FLANK_B_CHROM == FOCUS_CONTIG &                                                  #only report deletion when ends are known
                                     FLANK_B_CLOSE_BY == TRUE &                                                     #long distances likely represent joining two distant break ends, not long deletions
                                     FLANK_A_ORIENT == FLANK_B_ORIENT &                                             #in case inverted repeat
                                     FlankBUltStart <= FLANK_B_START_POS_DEL &                                      #this is required in case there is a tandem repeat
                                     FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_START_POS_DEL - (FlankAUltEnd+1)),
                                   FLANK_B_CHROM == FOCUS_CONTIG &                                                  #only report deletion when ends are known
                                     FLANK_B_CLOSE_BY == TRUE &                                                     #long distances likely represent joining two distant break ends, not long deletions
                                     FLANK_A_ORIENT == FLANK_B_ORIENT &                                             #in case inverted repeat
                                     FlankBUltStart >= FLANK_B_START_POS_DEL &                                      #this is required in case there is a tandem repeat
                                     FLANK_B_ORIENT == "RV" ~ as.integer((FlankAUltEnd-1) - FLANK_B_START_POS_DEL),
                                   FLANK_B_START_POS_DEL == FLANK_A_START_POS ~ as.integer(0),                      #in wt case
                                   TRUE ~ ERROR_NUMBER
                                   ))%>%

    #calculate total deletion length
    mutate(delSize = if_else(
      FLANK_B_DEL != ERROR_NUMBER,
      as.integer(FLANK_A_DEL + FLANK_B_DEL),
      ERROR_NUMBER
    )) %>%
    #correct MH if delSize == 0
    mutate(MH = if_else(delSize <= 0,
                        "",
                        MH)) %>%


    
    mutate(homologyLength = if_else(insSize == 0,
                                    nchar(MH),
                                    as.integer(-1))) %>% 
    
    #new logical column indicating orientation of flank B
    mutate(FLANK_B_ISFORWARD = if_else(FLANK_B_ORIENT=="FW",
                                TRUE,
                                FALSE)) %>%
    
    #fix the FLANK_B_START_POS again if the event is wt
    rowwise()%>%
    mutate(FLANK_B_START_POS = if (delSize == 0 & insSize == 0){
      if (FLANK_A_ORIENT == "FW"){
        FlankAUltEnd + 1
      }else{
        FlankAUltEnd - 1
      }
    }else{
      if (insSize > 0 ){      #for cases on the focus contig where the filler is a flank B continuation
        FLANK_B_START_POS_DEL
      }else{
      FLANK_B_START_POS_MH}
    }) %>%
    #remove reads with FLANK_B 's shorter than the minimum.
    filter(FLANK_B_LEN_MH >= FLANK_B_LEN_MIN)
  
  function_time("Step 6 took ")
  
  ###############################################################################
  #Process data: step 7
  ###############################################################################
  
  data_improved5b = data_improved5 %>%
  
    #for SIQPlotter delRelativeStart is FLANK_A_DEL (but negative) and delRelativeEnd is FLANK_B_DEL.
    mutate(
      delRelativeStart = as.integer(FLANK_A_DEL*-1),
      delRelativeEnd = as.integer(FLANK_B_DEL)) %>%
    
    #then do some checks to find wt events that because of mutations did not get called as such
    #first find a sequence around the DSB that would indicate no DSB has been made or is repaired perfectly
    #check in both seq1 and seq2
    rowwise() %>%
    mutate(DSB_AREA_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1))) %>%
    mutate(DSB_AREA_COUNT = length(DSB_AREA_CHECK@ranges))%>%
    mutate(DSB_AREA2_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ_RC), DNAString(SEQ_2_first), max.mismatch = 1))) %>%
    mutate(DSB_AREA2_COUNT = length(DSB_AREA2_CHECK@ranges))%>%
    ungroup()%>%
    mutate(DSB_AREA_HIT = "",
           DSB_AREA2_HIT = "",
           DSB_AREA_INTACT = FALSE,
           DSB_AREA_1MM = FALSE)
  
  for (j in row.names(data_improved5b)){
    j_int=as.integer(j)
    if (data_improved5b$DSB_AREA_COUNT[[j_int]]==1){
      if ((data_improved5b$SEQ_1_LEN[[j_int]] >= (data_improved5b$DSB_AREA_CHECK[[j_int]]@ranges@start[1]+data_improved5b$DSB_AREA_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved5b$DSB_AREA_CHECK[[j_int]]@ranges@start[1]>0)){
    data_improved5b$DSB_AREA_HIT[[j_int]] = as.character(data_improved5b$DSB_AREA_CHECK[[j_int]])
      }else{
        data_improved5b$DSB_AREA_HIT[[j_int]] = ""
      }
    }else{
      data_improved5b$DSB_AREA_HIT[[j_int]] = ""
    }
    if (data_improved5b$DSB_AREA2_COUNT[[j_int]]==1){
      if ((nchar(data_improved5b$SEQ_2_first[[j_int]]) >= (data_improved5b$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]+data_improved5b$DSB_AREA2_CHECK[[j_int]]@ranges@width[1]-1)) & (data_improved5b$DSB_AREA2_CHECK[[j_int]]@ranges@start[1]>0)){
        data_improved5b$DSB_AREA2_HIT[[j_int]] = as.character(data_improved5b$DSB_AREA2_CHECK[[j_int]])
      }else{
        data_improved5b$DSB_AREA2_HIT[[j_int]] = ""
      }
    }else{
      data_improved5b$DSB_AREA2_HIT[[j_int]] = ""
    }
  }

  
  
  
  data_improved6 =   data_improved5b %>%
    mutate(DSB_AREA_INTACT = if_else(DSB_AREA_HIT == DSB_AREA_SEQ | DSB_AREA2_HIT == DSB_AREA_SEQ_RC,
                                     "TRUE",
                                     "FALSE"))%>%
    mutate(DSB_AREA_1MM = if_else(
      DSB_AREA_INTACT == FALSE & (DSB_AREA_COUNT>0 | DSB_AREA2_COUNT>0),
      "TRUE",
      "FALSE")) %>%
    mutate(CASE_WT = if_else((
      DSB_AREA_INTACT==TRUE ),
      TRUE,
      FALSE)) %>%
    mutate(DSB_HIT_MULTI = if_else(DSB_AREA_COUNT>1 | DSB_AREA2_COUNT>1,
                                   "TRUE",
                                   "FALSE")) %>%
    
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
    mutate(read_minimum_length = if_else((DSB_AREA_COUNT==0 & FAKE_DELIN_CHECK == FALSE),
                                         FLANK_A_LEN + insSize + 30,
                                         FLANK_A_REF_LEN+30)) %>%
    filter(SEQ_1_LEN >= read_minimum_length)%>%
    mutate(SEQ_1_trimmed = substr(SEQ_1, 1, read_minimum_length))
  
  #here grouping will occur to determine the number of anchors.
  #for TRANSGUIDE a consensus outcome will be determined
         if (GROUPSAMEPOS == FALSE){
      
        data_improved8pre2 =data_improved6 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
      group_by(
        FILE_NAME,
        PRIMER_SEQ,
        delRelativeStart,
        delRelativeEnd,
        FLANK_A_LEN,
        FLANK_A_END_POS,
        FLANK_B_CHROM,
        FLANK_B_START_POS,
        MATE_FLANK_B_CHROM_AGREE,
        MATE_FLANK_B_CHROM,
        FLANK_B_ISFORWARD,
        FLANK_B_ORIENT,
        FILLER,
        MH,
        insSize,
        delSize,
        homologyLength,
        FAKE_DELIN_CHECK,
        DSB_AREA_INTACT,
        DSB_AREA_1MM,
        DSB_HIT_MULTI,
        TRIM_LEN
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
        
        data_improved8pre2 =data_improved6 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
          group_by(
            FILE_NAME,
            PRIMER_SEQ,
            FLANK_B_CHROM,
            FLANK_B_START_POS,
            FLANK_B_ISFORWARD,
            FLANK_B_ORIENT,
            TRIM_LEN,
            MATE_FLANK_B_CHROM_AGREE,
            MATE_FLANK_B_CHROM
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
            FLANK_A_END_POS_con = as.integer(names(which.max(table(FLANK_A_END_POS)))),
            FILLER_con = names(which.max(table(FILLER))),
            MH_con =names(which.max(table(MH))),
            insSize_con = as.integer(names(which.max(table(insSize)))),
            delSize_con = as.integer(names(which.max(table(delSize)))),
            homologyLength_con = as.integer(names(which.max(table(homologyLength)))),
            FAKE_DELIN_CHECK_con = as.logical(names(which.max(table(FAKE_DELIN_CHECK)))),
            DSB_AREA_INTACT_con = as.logical(names(which.max(table(DSB_AREA_INTACT)))),
            DSB_AREA_1MM_con = as.logical(names(which.max(table(DSB_AREA_1MM)))),
            DSB_HIT_MULTI_con = as.logical(names(which.max(table(DSB_HIT_MULTI)))),
            .groups="drop"
          )%>%
          mutate(Consensus_freq = Count_consensus/ReadCount)%>%
          #rename columns to that code below is compatible with TRANS and CISGUIDE
          rename(FILLER = FILLER_con,
                 MH = MH_con,
                 FAKE_DELIN_CHECK = FAKE_DELIN_CHECK_con,
                 DSB_AREA_INTACT = DSB_AREA_INTACT_con,
                 DSB_AREA_1MM = DSB_AREA_1MM_con,
                 DSB_HIT_MULTI = DSB_HIT_MULTI_con,
                 delRelativeStart = delRelativeStart_con,
                 delRelativeEnd = delRelativeEnd_con,
                 insSize = insSize_con,
                 delSize = delSize_con,
                 homologyLength = homologyLength_con,
                 FLANK_A_END_POS = FLANK_A_END_POS_con)
      }
        
    data_improved8 =data_improved8pre2 %>%
      
      #then calculate the difference between start of flank B (in the read) and the position of of flank B at the end of the mate.
      mutate(ANCHOR_DIST = case_when(FLANK_B_ORIENT=="FW" & FLANK_B_CHROM == MATE_FLANK_B_CHROM & MATE_B_END_POS_max > FLANK_B_START_POS ~ as.integer(1+MATE_B_END_POS_max - FLANK_B_START_POS),
                                     FLANK_B_ORIENT=="RV" & FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_START_POS > MATE_B_END_POS_min ~ as.integer(1+FLANK_B_START_POS - MATE_B_END_POS_min),
                                     TRUE ~ NF_NUMBER)) %>%

      #determine whether there is a translocation, and if so, change some other variables to support translocation plotting in SIQplotteR.
      #in the case of T-DNA or other transfected DNA, the real end of the molecule is known, and in that case the deletion length can be determined
      #in other translocation cases the delRelativeEnd becomes 0 for practical plotting reasons, not because it is really 0 (it is unknowable).
      #if the T-DNA backbone, it will also be 0. Firstly because it is unclear whether the backbone got transferred via border skipping or as a separate "T-DNA", secondly because the map is linear and breaks in the backbone. The calculation therefore is more difficult.
      mutate(Translocation = case_when(delSize == ERROR_NUMBER ~ TRUE,
                                       TRUE ~ FALSE))%>%
      mutate(Translocation_del_resolved = case_when(
                (FOCUS_LOCUS == "LB" | FOCUS_LOCUS == "RB") & FLANK_B_CHROM == DSB_CONTIG ~ TRUE, #if transguide, del can be calculated on the induced genomic break locus
                Translocation == TRUE & FOCUS_LOCUS != "LB" & FOCUS_LOCUS != "RB" & FLANK_B_CHROM == DSB_CONTIG ~ TRUE, #translocation with the same chromosome
                Translocation == TRUE & FLANK_B_CHROM == PLASMID & TDNA_IS_LBRB == TRUE & FLANK_B_START_POS >= TDNA_LB_END & FLANK_B_START_POS <= TDNA_RB_END ~ TRUE, #if cisguide or transguide with translocation, and within the T-DNA
                Translocation == TRUE & FLANK_B_CHROM == PLASMID & TDNA_IS_LBRB == FALSE & FLANK_B_START_POS >= TDNA_RB_END & FLANK_B_START_POS <= TDNA_LB_END ~ TRUE, #if cisguide or transguide with translocation, and within the T-DNA
                Translocation == TRUE & FLANK_B_CHROM == PLASMID_ALT & TDNA_ALT_IS_LBRB == TRUE & FLANK_B_START_POS >= TDNA_ALT_LB_END & FLANK_B_START_POS <= TDNA_ALT_RB_END ~ TRUE, #if cisguide or transguide with translocation, and within the T-DNA_ALT
                Translocation == TRUE & FLANK_B_CHROM == PLASMID_ALT & TDNA_ALT_IS_LBRB == FALSE & FLANK_B_START_POS >= TDNA_ALT_RB_END & FLANK_B_START_POS <= TDNA_ALT_LB_END ~ TRUE, #if cisguide or transguide with translocation, and within the T-DNA_ALT 
                                         TRUE ~ FALSE))%>% #else is some other genomic location or plasmid backbone
    mutate(delRelativeEnd = case_when(
               (FOCUS_LOCUS == "LB" | FOCUS_LOCUS == "RB") & FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ISFORWARD == TRUE & FLANK_B_START_POS >= (1+DSB_FW_END-DSB_OVERHANG) ~ FLANK_B_START_POS - (1+DSB_FW_END-DSB_OVERHANG),
               (FOCUS_LOCUS == "LB" | FOCUS_LOCUS == "RB") & FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ISFORWARD == FALSE & DSB_FW_END >= FLANK_B_START_POS  ~ DSB_FW_END - FLANK_B_START_POS,
               
               Translocation == TRUE & FOCUS_LOCUS != "LB" & FOCUS_LOCUS != "RB" & FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ORIENT == "FW" & FLANK_B_START_POS > (1+DSB_FW_END-DSB_OVERHANG) ~ (FLANK_B_START_POS - (1+DSB_FW_END-DSB_OVERHANG)), #if flank B is on right side of break
               Translocation == TRUE & FOCUS_LOCUS != "LB" & FOCUS_LOCUS != "RB" & FLANK_B_CHROM == DSB_CONTIG & FLANK_B_ORIENT == "RV" & FLANK_B_START_POS < (DSB_FW_END+DSB_OVERHANG) ~ (DSB_FW_END+DSB_OVERHANG) - FLANK_B_START_POS, #if flank b is on left side of break
               
               Translocation == TRUE & FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == TRUE & FLANK_B_START_POS >= TDNA_LB_END & FLANK_B_START_POS < TDNA_RB_END ~ FLANK_B_START_POS - TDNA_LB_END, #deletion from the LB
               Translocation == TRUE & FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == TRUE & TDNA_IS_LBRB == FALSE & FLANK_B_START_POS >= TDNA_RB_END & FLANK_B_START_POS < TDNA_LB_END ~ FLANK_B_START_POS - TDNA_RB_END, #deletion from the RB
               Translocation == TRUE & FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == TRUE & TDNA_RB_END >= FLANK_B_START_POS & TDNA_LB_END < FLANK_B_START_POS ~ TDNA_RB_END - FLANK_B_START_POS, #deletion from the RB
               Translocation == TRUE & FLANK_B_CHROM == PLASMID & FLANK_B_ISFORWARD == FALSE & TDNA_IS_LBRB == FALSE & TDNA_LB_END >= FLANK_B_START_POS & TDNA_RB_END < FLANK_B_START_POS ~ TDNA_LB_END - FLANK_B_START_POS, #deletion from the LB
                                        TRUE ~ delRelativeEnd))%>%
             #correct in the case that the delsize is larger than the maximum allowed delsize.
        mutate(Translocation_del_resolved = case_when(delRelativeEnd > MAXTARGETLENGTH ~ FALSE,
                                                    TRUE ~ Translocation_del_resolved),
               delRelativeEnd = case_when(delRelativeEnd > MAXTARGETLENGTH ~ 0,
                                          TRUE ~ delRelativeEnd),
               
             delSize = case_when(delSize == ERROR_NUMBER ~ (delRelativeStart*-1)+delRelativeEnd,
                                 TRUE ~ delSize)) %>%
      
      
      
      #Add a Subject and Type column. Also add two extra columns for SIQplotteR that are equal to the other del columns.
    mutate(
      Subject = FOCUS_LOCUS,
      Type = case_when(
        (delSize == 0 & insSize == 0) | DSB_AREA_INTACT==TRUE | DSB_AREA_1MM==TRUE | FAKE_DELIN_CHECK == TRUE ~ "WT",
        delSize != 0 & delSize != ERROR_NUMBER & insSize == 0 ~ "DELETION",
        insSize != 0 & delSize == 0 ~ "INSERTION",
        delSize != 0 & delSize != ERROR_NUMBER & insSize != 0 ~ "DELINS",
        TRUE ~ "OTHER"),
      delRelativeStartTD = delRelativeStart,
      delRelativeEndTD = delRelativeEnd
    ) %>%
      
    
    
    #select the most important columns
    select(
      Name,
      FILE_NAME,
      PRIMER_SEQ,
      delRelativeStart,
      delRelativeStartTD,
      delRelativeEnd,
      delRelativeEndTD,
      AnchorCount,
      ANCHOR_DIST,
      ReadCount,
      FLANK_B_CHROM,
      FLANK_B_START_POS,
      FLANK_B_ISFORWARD,
      FLANK_B_ORIENT,
      FILLER,
      MH,
      insSize,
      homologyLength,
      delSize,
      Type,
      Subject,
      SEQ_1_con,
      SEQ_2_con,
      DSB_AREA_INTACT,
      DSB_AREA_1MM,
      DSB_HIT_MULTI,
      TRIM_LEN,
      Consensus_freq,
      FAKE_DELIN_CHECK,
      MATE_FLANK_B_CHROM_AGREE,
      MATE_FLANK_B_CHROM,
      FLANK_A_END_POS,
      Translocation,
      Translocation_del_resolved
    ) 
    

    
  
  function_time("Step 7 took ")
  
  ###############################################################################
  #Process data: step 8
  ###############################################################################
  
  #calculate the fraction of reads with a certain outcome within a library
  data_improved9 = data_improved8 %>% group_by(FILE_NAME) %>% summarize(ReadCountTotal =
                                                                                  sum(ReadCount),
                                                                        .groups="drop")
  function_time("Step 8 took ")
  
  ###############################################################################
  #Process data: step 9
  ###############################################################################
  
  data_improved10 = left_join(data_improved8, data_improved9, by = "FILE_NAME") %>%
    mutate(fraction = ReadCount / ReadCountTotal) %>%
    #add columns for all the input options and software version
    mutate(countReadsTotal = NULL,
           FlankAUltEnd = FlankAUltEnd,
           FlankBUltStart = FlankBUltStart,
           Flank_A_orient = FLANK_A_ORIENT,
           Focus_contig = FOCUS_CONTIG,
           Genotype = Genotype,
           DNASample = DNASample,
           Ecotype = Ecotype,
           RunID = RunID,
           program_version = hash,
           Plasmid = PLASMID,
           Plasmid_alt = PLASMID_ALT,
           AgroGeno = AgroGeno,
           RemoveNonTranslocation = REMOVENONTRANS,
           GroupSamePosition = GROUPSAMEPOS,
           Primer_match_perfect = Primer_match_perfect,
           Alias = paste0(Library, "_", RunID),
           DSB_FW_END = DSB_FW_END,
           DSB_OVERHANG = DSB_OVERHANG,
           DSB_CONTIG = DSB_CONTIG,
           TDNA_LB_END = TDNA_LB_END,
           TDNA_RB_END = TDNA_RB_END,
           TDNA_IS_LBRB = TDNA_IS_LBRB,
           TDNA_ALT_LB_END = TDNA_ALT_LB_END,
           TDNA_ALT_RB_END = TDNA_ALT_RB_END,
           TDNA_ALT_IS_LBRB = TDNA_ALT_IS_LBRB,
           MAXTARGETLENGTH = MAXTARGETLENGTH)
  

  
  
  
  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, data_improved10)
  saveWorkbook(work_book, file = paste0(output_dir, Sample, "_", RunID, "_CISTRANSGUIDE_V2.xlsx"), overwrite = TRUE)
  
  function_time("Step 9 took ")
  

  #show progress
  PercentageDone = PercentageDone + ((CurrentFileSize/TotalFileSize)*100)
  message(paste0("CISTRANSGUIDE analysis ", round(PercentageDone, digits=3), "% complete"))
  
}


###############################################################################
#Combine data: step 10
###############################################################################

sample_list = list.files(path=output_dir, pattern = "CISTRANSGUIDE_V2.xlsx")
wb_pre = tibble()

for (i in sample_list){
  wb_pre=bind_rows(wb_pre, read.xlsx(paste0(output_dir, i)) %>% select_if(function(x) !(all(is.na(x)) | all(x==""))))

}
#in case the following columns are NA, they will not have been imported from the excel. I therefore need to add them. But also allowing the possibility that they are already there.
missing_columns = wb_pre %>%
  select(Name) %>%
  mutate(Plasmid_alt = NA,
         DSB_FW_END = NA,
         DSB_OVERHANG = NA,
         DSB_CONTIG = NA,
         TDNA_ALT_LB_END = NA,
         TDNA_ALT_RB_END = NA,
         TDNA_ALT_IS_LBRB = NA)

wb = left_join(wb_pre, missing_columns)

#remove previously marked problematic events as well as duplicate positions
if (REMOVEPROBLEMS == TRUE) {
  message("removing problematic events")
  total_data_positioncompare_pre = wb %>%
    #first remove problematic events based on characteristics of the events themselves
    filter(AnchorCount >= ANCHORCUTOFF,
           ANCHOR_DIST >= MINANCHORDIST,
           ANCHOR_DIST <= MAXANCHORDIST,
           Consensus_freq >= 0.75,
           MATE_FLANK_B_CHROM_AGREE == TRUE,
           FAKE_DELIN_CHECK == FALSE,
           !(DSB_AREA_INTACT==TRUE & (delSize !=0 | insSize != 0 | homologyLength != 0))) %>% 
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
  ID_prev = 0
  for (i in 2:length(total_data_positioncompare_pre$Alias)){
    if (total_data_positioncompare_pre$same_as_prev[i] == TRUE){
      total_data_positioncompare_pre$ID[i] = ID_prev
      ID_prev = total_data_positioncompare_pre$ID[i]
    }else{
      total_data_positioncompare_pre$ID[i] = ID_prev + 1
      ID_prev = total_data_positioncompare_pre$ID[i]
    }
  }
  
  #add family info
  sample_info2  = sample_info %>%
    select(Family, RunID, Sample)%>%
    rowwise()%>%
    mutate(Alias = paste(Sample, RunID, sep="_"))%>%
    select(Alias, Family)%>%
    ungroup()
  
  total_data_positioncompare = left_join(sample_info2, total_data_positioncompare_pre, by=c("Alias"))%>%
    filter(!is.na(Family) & !is.na(Name))
  
  message("combining junctions with similar positions")
  #combine junctions with similar positions and get the characteristics of the consensus event from the event the most anchors 
  total_data_near_positioncombined = total_data_positioncompare %>%
    group_by(Alias, FLANK_B_CHROM, Plasmid, FLANK_B_ISFORWARD, DNASample, Subject, ID, Focus_contig, Genotype, Ecotype, Plasmid_alt, Family, FlankAUltEnd, FlankBUltStart, AgroGeno, RemoveNonTranslocation, GroupSamePosition, Translocation, Translocation_del_resolved, Primer_match_perfect, DSB_FW_END, DSB_OVERHANG, DSB_CONTIG, TDNA_LB_END, TDNA_RB_END, TDNA_IS_LBRB, TDNA_ALT_LB_END, TDNA_ALT_RB_END, TDNA_ALT_IS_LBRB, MAXTARGETLENGTH)%>%
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
  wb_family = sample_info2 %>% select(Family) %>% distinct() %>% filter(Family!=0)
  if (nrow(wb_family) == 0) {
    message("no family info detected")
    #if no families are indicated

    wb_flag = total_data_near_positioncombined %>%
        #then examine positions across samples and remove those that occur multiple times
        group_by(FLANK_B_START_POS) %>%
        mutate(duplicate_position = if_else(n() > 1,
                                            TRUE,
                                            FALSE)) %>%
        ungroup() %>%
        filter(duplicate_position == FALSE)
      
  
  } else{
    message("taking family into consideration")
    #take families into account
    wb_filter_total = total_data_near_positioncombined %>% filter(Family == 99999999) #make an empty file
    
    for (i in wb_family$Family) {
      message(paste0("Cleanup family ", i))
      #cleanup per family
      wb_current_family = total_data_near_positioncombined %>% filter(Family == i) %>% select(Alias) %>% distinct() #make a list of aliases within the current family
      wb_filter_subtotal = total_data_near_positioncombined %>% filter(Family == 99999999) #make an empty file
      
      for (j in wb_current_family$Alias) {
        #per alias in that family
        
        wb_filter_current = total_data_near_positioncombined %>%
          filter(Family != i | Alias == j) %>% #events are either not of the current family, or they belong to the current alias
          group_by(FLANK_B_START_POS) %>%
          mutate(duplicate_position = if_else(n() > 1,
                                              TRUE,
                                              FALSE)) %>%
          
          ungroup() %>%
          filter(duplicate_position == FALSE &
                   Family == i) #keep only events belonging to the current file, and remove duplicate positions
        
        wb_filter_subtotal = rbind(wb_filter_subtotal, wb_filter_current) #combine surviving events from the current family
        
      }
      wb_filter_total = rbind(wb_filter_total, wb_filter_subtotal) #combining surviving events from all families
      
    }
    message("Fetching non-duplicate position events not belonging to a family")
    wb_nonfamily = total_data_near_positioncombined %>%
      filter(Family == 0) %>%
      group_by(FLANK_B_START_POS) %>%
      mutate(duplicate_position = if_else(n() > 1,
                                          TRUE,
                                          FALSE)) %>%
      
      ungroup() %>%
      filter(duplicate_position == FALSE)
    message("combining surviving family and nonfamily events")
    wb_flag = rbind(wb_filter_total, wb_nonfamily)  
    
  }
} else{
  message("flagging problems only")
  wb_flag = wb %>%
    group_by(FLANK_B_START_POS) %>%
    mutate(duplicate_position = if_else(n() > 1,
                                        TRUE,
                                        FALSE)) %>%
    
    ungroup() 
}

wb_flag1 = wb_flag %>%
  mutate(RemoveProblematicEvents = REMOVEPROBLEMS)%>%
  #add/ change several things for compatibility with SIQplotteR
  mutate(getHomologyColor = "dummy",
         Barcode = "dummy")%>%
  rename(countEvents = ReadCount,
         insertion = FILLER,
         homology = MH)

#correct all WT events 
if (CONVERTWT == TRUE){
  wb_flag2 = wb_flag1 %>%
    mutate(delRelativeStart = if_else(Type=="WT",
                                    0,
                                    as.integer(delRelativeStart)),
           delRelativeStartTD = if_else(Type=="WT",
                                      0,
                                      as.integer(delRelativeStart)),
           delRelativeEnd = if_else(Type=="WT",
                                    0,
                                    as.integer(delRelativeEnd)),
           delRelativeEndTD = if_else(Type=="WT",
                                    0,
                                    as.integer(delRelativeEnd)),
           delSize = if_else(Type=="WT",
                                    0,
                                    as.integer(delSize)),
           insSize = if_else(Type=="WT",
                                    0,
                                    as.integer(insSize)),
           homologyLength = if_else(Type=="WT",
                                    0,
                                    as.integer(homologyLength)),
           Translocation = if_else(Type=="WT",
                                   FALSE,
                                   as.logical(Translocation))
           )
}else{
  wb_flag2 = wb_flag1
}


message("Writing output")
work_book2 <- createWorkbook()
addWorksheet(work_book2, "rawData")
writeData(work_book2, sheet = 1, wb_flag2)

#Write an additional sheet with read number info
read_numbers_info = read.csv(paste0(input_dir, "read_numbers.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
wb_numbers = read_numbers_info %>% 
  mutate(Alias = paste0(Sample, "_", RunID))%>%
  mutate(Sample = NULL,
         RunID = NULL)

addWorksheet(work_book2, "Information")
writeData(work_book2, sheet = 2, wb_numbers)
saveWorkbook(work_book2, file = paste0(output_dir, "Data_combined_CISTRANSGUIDE_V2_", as.integer(Sys.time()), ".xlsx"), overwrite = TRUE)

message("CISTRANSGUIDE analysis has completed")
