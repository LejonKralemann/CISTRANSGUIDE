###############################################################################
#install and load packages
###############################################################################
if (require(BSgenome)==FALSE){install.packages("BSgenome")}
if (require(Biostrings)==FALSE){install.packages("Biostrings")}
if (require(stringi)==FALSE){install.packages("stringi")}
if (require(stringdist)==FALSE){install.packages("stringdist")}
if (require(tidyverse)==FALSE){install.packages("tidyverse")}
if (require(openxlsx)==FALSE){install.packages("openxlsx")}



###############################################################################
#set parameters
###############################################################################
input_dir= "./input/"
output_dir= "./output/"
hash=system("git rev-parse HEAD", intern=TRUE)
hash_little=substr(hash, 1, 8)
NF_NUMBER = as.integer(-99999999) #don't change
ERROR_NUMBER = as.integer(99999999) #don't change
MINBASEQUAL = 0.75 #minimum base quality
MAX_DIST_FLANK_B_END = 10000 #distance from end of flank B to DSB, determines max deletion size and also affects maximum insertion size
FLANK_B_LEN_MIN = 15 #minimum length of flank B. Also affects size of DSB_AREA_SEQ (2x FLANK_B_LEN_MIN)
LOCUS_WINDOW = 1000 #size of the window centered on the DSB, RB nick, or LB nick to determine locus info
sample_info = read.csv(paste0(input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
TIME_START=as.numeric(Sys.time())*1000
DEBUG=FALSE #if on, reads will not be discarded when a problem has been detected, but flagged.

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
    compare_len_1 = lcprefix(ref, seq1)
    compare_flank_ref = if ((compare_len_1 + 2)< nchar(ref)){
      substr(ref, start= compare_len_1+2, stop=nchar(ref))
    }else{
      ""}
    compare_seq_1 = if (compare_flank_ref != ""){
      substr(seq1, start = compare_len_1 +2, stop = nchar(seq1))
    }else{
      ""}
    compare_len_2 = if (compare_flank_ref != ""){
      as.integer(lcprefix(compare_flank_ref, compare_seq_1))
    }else{
      as.integer(0)}
    flank_match = if (compare_len_2 > 9){
      substr(ref, start = 1, stop = compare_len_1 + compare_len_2 + 1)
    }else{
      substr(ref, start = 1, stop = compare_len_1)}
    return(flank_match)
  }
}
function_time <-function(text){
  TIME_CURRENT=as.numeric(Sys.time())*1000
  message(paste0(text, (TIME_CURRENT - TIME_START), " milliseconds"))
  TIME_START<<-as.numeric(Sys.time())*1000
}

###############################################################################
#Process data: step 1
###############################################################################

for (i in row.names(sample_info)){
  Sample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  RunID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(RunID))
  
  if (file.exists(paste0(input_dir, Sample, "_", RunID, "_A.txt"))==FALSE){
    message("Primary processed file not found, moving to next sample")
    next
  }else{
    message(paste0("Processing ",input_dir, Sample, "_", RunID, "_A.txt"))
  }
  
  data = read.csv(paste0(input_dir, Sample, "_", RunID, "_A.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
  if (nrow(data)==0){
    message("Primary processed file empty, moving to next sample")
    next
  }
  
  FOCUS_CONTIG = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DSB_chrom))
  FOCUS_LOCUS = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Locus_name))
  Genotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Genotype))
  PLASMID = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid))
  PLASMID_ALT = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Plasmid_alt))
  FlankAUltEnd = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FlankAUltEnd))
  FlankBUltStart = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FlankBUltStart))
  DNASample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DNA))
  Ecotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ecotype))
  Library = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  
  if (file.exists(paste0(input_dir, PLASMID, ".fa"))==FALSE){
    message("Reference fasta not found")
    next
  }
  
  genomeseq = readDNAStringSet(paste0(input_dir, PLASMID,".fa") , format="fasta")
  contig_seq = as.character(eval(parse(text = paste0("genomeseq$", FOCUS_CONTIG))))
  Primer_seq = str_replace_all(as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Primer)), "TCAGACGTGTGCTCTTCCGATCT", "")
  Primer_match = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  FLANK_A_ORIENT = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FLANK_A_ORIENT))

  
  DSB_AREA_SEQ = (if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= FlankAUltEnd - (FLANK_B_LEN_MIN-1), stop= FlankAUltEnd + FLANK_B_LEN_MIN)
  }else if (FLANK_A_ORIENT == "RV"){
    as.character(reverseComplement(DNAString(substr(contig_seq, start= FlankAUltEnd - FLANK_B_LEN_MIN, stop= FlankAUltEnd + (FLANK_B_LEN_MIN-1)))))
  }else{
    ""
  })
  
  if (FLANK_A_ORIENT == "FW"){
    Primer_match = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))  
    if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
      Primer_pos = as.numeric(Primer_match$start)
    }else if (nrow(Primer_match) >1){
      message("Primer found several times in the genome")
      next
    }else{
      message("Primer not found")
      next
    }
  }else if (FLANK_A_ORIENT == "RV"){
    Primer_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
    if (nrow(Primer_match) > 0 & nrow(Primer_match) < 2){
      Primer_pos = as.numeric(Primer_match$end)
    }else if (nrow(Primer_match) >1){
      message("Primer found several times in the genome")
      next
    }else{
      message("Primer not found")
      next
    }
  }else{
    next
  }
  #calculate the length from primer to DSB
  PRIMER_TO_DSB = if (FLANK_A_ORIENT == "FW"){
    FlankAUltEnd - (Primer_pos -1)
  }else if (FLANK_A_ORIENT=="RV"){
    Primer_pos - (FlankAUltEnd -1)
  }else{
    ERROR_NUMBER
  }
  #set the minimum length of a read
  MINLEN = PRIMER_TO_DSB+FLANK_B_LEN_MIN

  #get the REF seq for flank A. from primer start to DSB +3 if RV primer, not if FW primer. Because CAS9 can cut further away from the PAM, but not closer. So the FLANK_A_REF is going as far as FLANK A is allowed to go.
  FLANK_A_REF = if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= Primer_pos, stop= FlankAUltEnd)
  }else if (FLANK_A_ORIENT=="RV"){
    as.character(reverseComplement(DNAString(
      substr(contig_seq, start= FlankAUltEnd, stop= Primer_pos))
    ))
  }else{
    ""
  }
  GLOBAL_TOTAL_REF = if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= Primer_pos, stop= Primer_pos+MAX_DIST_FLANK_B_END+PRIMER_TO_DSB)
  }else if (FLANK_A_ORIENT=="RV"){
    as.character(reverseComplement(DNAString(
      substr(contig_seq, start= Primer_pos-(MAX_DIST_FLANK_B_END+PRIMER_TO_DSB), stop= Primer_pos))
    ))
  }else{
    ""
  } 
  
  function_time("Step 1 took ")

  ###############################################################################
  #Process data: step 2
  ###############################################################################
  
  data_improved  = data %>%
    
    #filter(QNAME == "A00379:436:H3CHWDMXY:1:2448:17381:24674") %>%
    
    #Count number of Ns and remove any reads with Ns
    mutate(NrN = str_count(SEQ_1, pattern = "N"),
           SEQ_1_LEN = nchar(SEQ_1)) %>%
    filter(NrN < 1) %>%
    
    #filter away reads that are too short
    filter(!(SEQ_1_LEN < MINLEN)) %>%
    
    #calculate the average base quality, and filter away low quality reads
    mutate(AvgBaseQual_1 =   ((
      (str_count(QUAL_1, pattern = '\\"') * 1) +
        (str_count(QUAL_1, pattern = "#") * 2) +
        (str_count(QUAL_1, pattern = "\\$") * 3) +
        (str_count(QUAL_1, pattern = "%") * 4) +
        (str_count(QUAL_1, pattern = "&") * 5) +
        (str_count(QUAL_1, pattern = "'") * 6) +
        (str_count(QUAL_1, pattern = "\\(") * 7) +
        (str_count(QUAL_1, pattern = "\\)") * 8) +
        (str_count(QUAL_1, pattern = "\\*") * 9) +
        (str_count(QUAL_1, pattern = "\\+") * 10) +
        (str_count(QUAL_1, pattern = ",") * 11) +
        (str_count(QUAL_1, pattern = "\\-") * 12) +
        (str_count(QUAL_1, pattern = "\\.") * 13) +
        (str_count(QUAL_1, pattern = "/") * 14) +
        (str_count(QUAL_1, pattern = "0") * 15) +
        (str_count(QUAL_1, pattern = "1") * 16) +
        (str_count(QUAL_1, pattern = "2") * 17) +
        (str_count(QUAL_1, pattern = "3") * 18) +
        (str_count(QUAL_1, pattern = "4") * 19) +
        (str_count(QUAL_1, pattern = "5") * 20) +
        (str_count(QUAL_1, pattern = "6") * 21) +
        (str_count(QUAL_1, pattern = "7") * 22) +
        (str_count(QUAL_1, pattern = "8") * 23) +
        (str_count(QUAL_1, pattern = "9") * 24) +
        (str_count(QUAL_1, pattern = ":") * 25) +
        (str_count(QUAL_1, pattern = ";") * 26) +
        (str_count(QUAL_1, pattern = "<") * 27) +
        (str_count(QUAL_1, pattern = "=") * 28) +
        (str_count(QUAL_1, pattern = ">") * 29) +
        (str_count(QUAL_1, pattern = "\\?") * 30) +
        (str_count(QUAL_1, pattern = "@") * 31) +
        (str_count(QUAL_1, pattern = "A") * 32) +
        (str_count(QUAL_1, pattern = "B") * 33) +
        (str_count(QUAL_1, pattern = "C") * 34) +
        (str_count(QUAL_1, pattern = "D") * 35) +
        (str_count(QUAL_1, pattern = "E") * 36) +
        (str_count(QUAL_1, pattern = "F") * 37) +
        (str_count(QUAL_1, pattern = "G") * 38) +
        (str_count(QUAL_1, pattern = "H") * 39) +
        (str_count(QUAL_1, pattern = "I") * 40)
    ) / (40 * nchar(QUAL_1))
    )) %>%
    filter(AvgBaseQual_1 > MINBASEQUAL)%>%
    group_by(
      RNAME_1,
      POS_1,
      CIGAR_1,
      SATAG_1,
      SEQ_RCed_1,
      SEQ_1,
      RNAME_2,
      POS_2,
      CIGAR_2,
      SATAG_2,
      SEQ_RCed_2,
      FILE_NAME,
      PRIMER_SEQ,
      SEQ_1_LEN,
      DEDUP_METHOD,
      TRIM_LEN) %>%
    summarize(
      countEvents_init = n(),
      Max_BaseQual_1 = max(AvgBaseQual_1),
      QNAME_first = first(QNAME),
      SEQ_2_first = first(SEQ_2))
  
  function_time("Step 2 took ")
  
  ###############################################################################
  #Process data: step 3
  ###############################################################################
  
  data_improved1 = data_improved %>%
    
    #check for expected position and orientation base on primer seqs
    mutate(
      FLANK_A_START_POS = Primer_pos
    )%>%
    
    #calculate the number of alignments
    mutate(CIGAR_1_LEN = nchar(str_remove_all(CIGAR_1, "[1234567890]"))) %>%
    mutate(CIGAR_2_LEN = nchar(str_remove_all(CIGAR_2, "[1234567890]"))) %>%
    
    #split all the CIGARs and SATags for later
    rowwise() %>%
    mutate(
      CIGAR_1_N = str_split(CIGAR_1, pattern = "([MSHID])"),
      CIGAR_1_rm = str_remove_all(CIGAR_1, "[1234567890]"),
      CIGAR_1_L = str_split(CIGAR_1_rm, pattern = ""),
      CIGAR_1_rm = NULL,
      CIGAR_1_N = list(stri_remove_empty(CIGAR_1_N))
    ) %>%
    ungroup() %>%
    separate(
      col = SATAG_1,
      remove = TRUE,
      into = c(
        "SATAG_1_1",
        "SATAG_1_2",
        "SATAG_1_3",
        "SATAG_1_4",
        "SATAG_1_5",
        "SATAG_1_6",
        "SATAG_1_7",
        "SATAG_1_8",
        "SATAG_1_9",
        "SATAG_1_10",
        "SATAG_1_11",
        "SATAG_1_12",
        "SATAG_1_13",
        "SATAG_1_14"
      ),
      sep = "([:,;])"
    ) %>%
    mutate(
      SATAG_1_2 = NULL,
      SATAG_1_7 = NULL,
      SATAG_1_8 = NULL,
      SATAG_1_13 = NULL,
      SATAG_1_14 = NULL,
      SATAG_1_6_CIGARLEN = nchar(str_remove_all(SATAG_1_6, "[1234567890]")),
      SATAG_1_12_CIGARLEN = nchar(str_remove_all(SATAG_1_12, "[1234567890]"))
    ) %>%
    
    rowwise() %>%
    mutate(
      SATAG_1_6_N = str_split(SATAG_1_6, pattern = "([MSHID])"),
      SATAG_1_6_rm = str_remove_all(SATAG_1_6, "[1234567890]"),
      SATAG_1_6_L = str_split(SATAG_1_6_rm, pattern = ""),
      SATAG_1_6_rm = NULL,
      SATAG_1_6_N = list(stri_remove_empty(SATAG_1_6_N)),
      SATAG_1_12_N = str_split(SATAG_1_12, pattern = "([MSHID])"),
      SATAG_1_12_rm = str_remove_all(SATAG_1_12, "[1234567890]"),
      SATAG_1_12_L = str_split(SATAG_1_12_rm, pattern = ""),
      SATAG_1_12_rm = NULL,
      SATAG_1_12_N = list(stri_remove_empty(SATAG_1_12_N))
    )  %>%
    ungroup() %>%
    
    #make certain columns integers, for calculating
    mutate(
      POS_1 = as.integer(POS_1),
      SATAG_1_4 = as.integer(SATAG_1_4),
      SATAG_1_10 = as.integer(SATAG_1_10)
    ) %>%
    
    #do the same stuff, but now for the mate
    rowwise() %>%
    mutate(
      CIGAR_2_N = str_split(CIGAR_2, pattern = "([MSHID])"),
      CIGAR_2_rm = str_remove_all(CIGAR_2, "[1234567890]"),
      CIGAR_2_L = str_split(CIGAR_2_rm, pattern = ""),
      CIGAR_2_rm = NULL,
      CIGAR_2_N = list(stri_remove_empty(CIGAR_2_N))
    ) %>%
    ungroup() %>%
    separate(
      col = SATAG_2,
      remove = TRUE,
      into = c(
        "SATAG_2_1",
        "SATAG_2_2",
        "SATAG_2_3",
        "SATAG_2_4",
        "SATAG_2_5",
        "SATAG_2_6",
        "SATAG_2_7",
        "SATAG_2_8",
        "SATAG_2_9",
        "SATAG_2_10",
        "SATAG_2_11",
        "SATAG_2_12",
        "SATAG_2_13",
        "SATAG_2_14"
      ),
      sep = "([:,;])"
    ) %>%
    mutate(
      SATAG_2_2 = NULL,
      SATAG_2_7 = NULL,
      SATAG_2_8 = NULL,
      SATAG_2_13 = NULL,
      SATAG_2_14 = NULL,
      SATAG_2_6_CIGARLEN = nchar(str_remove_all(SATAG_2_6, "[1234567890]")),
      SATAG_2_12_CIGARLEN = nchar(str_remove_all(SATAG_2_12, "[1234567890]"))
    ) %>%
    
    rowwise() %>%
    mutate(
      SATAG_2_6_N = str_split(SATAG_2_6, pattern = "([MSHID])"),
      SATAG_2_6_rm = str_remove_all(SATAG_2_6, "[1234567890]"),
      SATAG_2_6_L = str_split(SATAG_2_6_rm, pattern = ""),
      SATAG_2_6_rm = NULL,
      SATAG_2_6_N = list(stri_remove_empty(SATAG_2_6_N)),
      SATAG_2_12_N = str_split(SATAG_2_12, pattern = "([MSHID])"),
      SATAG_2_12_rm = str_remove_all(SATAG_2_12, "[1234567890]"),
      SATAG_2_12_L = str_split(SATAG_2_12_rm, pattern = ""),
      SATAG_2_12_rm = NULL,
      SATAG_2_12_N = list(stri_remove_empty(SATAG_2_12_N))
    )  %>%
    ungroup() %>%
    
    #make certain columns integers, for calculating
    mutate(
      POS_2 = as.integer(POS_2),
      SATAG_2_4 = as.integer(SATAG_2_4),
      SATAG_2_10 = as.integer(SATAG_2_10)
    ) %>%
  
  #split some columns, keep original columns for now
    rowwise()%>%
    mutate(
    CIGAR_1_L_first_element = head(CIGAR_1_L, 1),
    CIGAR_1_L_last_element = tail(CIGAR_1_L, 1),
    SATAG_1_6_L_last_element = tail(SATAG_1_6_L, 1),
    SATAG_1_6_L_first_element = head(SATAG_1_6_L, 1),
    CIGAR_1_N_first_element = as.integer(head(CIGAR_1_N, 1)),
    CIGAR_1_N_last_element = as.integer(tail(CIGAR_1_N, 1)),
    SATAG_1_6_N_last_element = as.integer(tail(SATAG_1_6_N, 1)),
    SATAG_1_6_N_first_element = as.integer(head(SATAG_1_6_N, 1)),
    CIGAR_1_N_second_element = as.integer(tail(head(CIGAR_1_N, 2), 1)),
    CIGAR_1_L_second_element = tail(head(CIGAR_1_L, 2), 1),
    CIGAR_2_L_first_element = head(CIGAR_2_L, 1),
    CIGAR_2_L_last_element = tail(CIGAR_2_L, 1),
    SATAG_2_6_L_first_element = head(SATAG_2_6_L, 1),
    SATAG_2_6_L_last_element = tail(SATAG_2_6_L, 1),
    SATAG_2_12_L_first_element = head(SATAG_2_12_L, n = 1),
    SATAG_2_12_L_last_element = tail(SATAG_2_12_L, n = 1)
  )%>%
    ungroup()
  
  
  function_time("Step 3 took ")
  
  ###############################################################################
  #Process data: step 4
  ###############################################################################
  
  data_improved2 = data_improved1 %>%
    
    #temp filters
    filter(SATAG_1_1 != "XA") %>%
    filter(SATAG_2_1 != "XA") %>%
    
    #calculate read span length. This is basically equal to read length, plus deletion length within an alignment.
    rowwise() %>%
    mutate(READ_SPAN = sum(as.integer(CIGAR_1_N))) %>%
    ungroup() %>%
    
    #determine which chromosome the B flank is aligning to
    mutate(
      FLANK_B_CHROM = (case_when(
        (CIGAR_1_LEN == 1) ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == TRUE & CIGAR_1_L_first_element == "M") ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == TRUE  & CIGAR_1_L_first_element == "S" & SATAG_1_1 != "SA") ~ "NOT_FOUND",
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == FALSE & CIGAR_1_L_last_element == "M") ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == FALSE  & CIGAR_1_L_last_element == "S" & SATAG_1_1 != "SA") ~ "NOT_FOUND",
        (SEQ_RCed_1 == TRUE & SATAG_1_1 == "SA" & CIGAR_1_L_first_element == "S" & SATAG_1_6_L_last_element == "M" & SATAG_1_5 == "+") ~ SATAG_1_3,
        (SEQ_RCed_1 == TRUE & SATAG_1_1 == "SA" & CIGAR_1_L_first_element == "S" & SATAG_1_6_L_last_element == "S" & SATAG_1_5 == "+" & is.na(SATAG_1_12)) ~ "NOT_FOUND",
        (
          SEQ_RCed_1 == TRUE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_first_element == "S" &
            SATAG_1_6_L_first_element == "M" & SATAG_1_5 == "-"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == TRUE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_first_element == "S" &
            SATAG_1_6_L_first_element == "S" &
            SATAG_1_5 == "-" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_last_element == "S" &
            SATAG_1_6_L_last_element == "M" & SATAG_1_5 == "+"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_last_element == "S" &
            SATAG_1_6_L_last_element == "S" &
            SATAG_1_5 == "+" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_last_element == "S" &
            SATAG_1_6_L_first_element == "M" & SATAG_1_5 == "-"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            CIGAR_1_L_last_element == "S" &
            SATAG_1_6_L_first_element == "S" &
            SATAG_1_5 == "-" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        TRUE ~ "ERROR"
      ))
    ) %>%
    
    #mark cases where flank B chrom is not found
    mutate(
      hasPROBLEM = if_else(FLANK_B_CHROM != "NOT_FOUND",
                           FALSE,
                           TRUE)) %>%
    
    #write down the various potential insertions
    rowwise() %>%
    mutate(
      POS2_I = case_when(
        CIGAR_1_LEN > 2  &
          CIGAR_1_L_second_element == "I" ~ CIGAR_1_N_second_element,
        TRUE ~ as.integer(0)
      ),
      POS3_I = case_when(
        CIGAR_1_LEN > 3  &
          tail(head(CIGAR_1_L, 3), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 3), 1)),
        TRUE ~ as.integer(0)
      ),
      POS4_I = case_when(
        CIGAR_1_LEN > 4  &
          tail(head(CIGAR_1_L, 4), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 4), 1)),
        TRUE ~ as.integer(0)
      ),
      POS5_I = case_when(
        CIGAR_1_LEN > 5  &
          tail(head(CIGAR_1_L, 5), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 5), 1)),
        TRUE ~ as.integer(0)
      ),
      POS6_I = case_when(
        CIGAR_1_LEN > 6  &
          tail(head(CIGAR_1_L, 6), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 6), 1)),
        TRUE ~ as.integer(0)
      ),
      POS7_I = case_when(
        CIGAR_1_LEN > 7  &
          tail(head(CIGAR_1_L, 7), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 7), 1)),
        TRUE ~ as.integer(0)
      )
    ) %>%
    mutate(READ_SPAN_MINUS_I = as.integer(READ_SPAN - (POS2_I + POS3_I + POS4_I + POS5_I + POS6_I + POS7_I)))%>%
    ungroup()
  
  function_time("Step 4 took ")
  
  ###############################################################################
  #Process data: step 5
  ###############################################################################
  
  data_improved3 = data_improved2 %>%
    
    #calculate the end position of the B flank, in order to get the ref seq. This is the end the furthest away from the junction.
    
    mutate(
      FLANK_B_END_POS = case_when(
        CIGAR_1_LEN == 1 & SEQ_RCed_1 == TRUE ~ POS_1,
        CIGAR_1_LEN == 1 &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - 1),
        CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "M" & SEQ_RCed_1 == TRUE ~ POS_1,
        CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "M" &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - 1),
        
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == TRUE &
          CIGAR_1_L_first_element == "M" & CIGAR_1_L_last_element == "S" ~ POS_1,
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "S" &
          CIGAR_1_L_last_element == "M" ~ as.integer(POS_1 + CIGAR_1_N_last_element -1),
        
        SEQ_RCed_1 == TRUE &
          CIGAR_1_L_first_element == "S" &
          CIGAR_1_L_last_element == "M" & SATAG_1_1 != "SA" ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "S"  & SATAG_1_1 != "SA" ~ NF_NUMBER,
        
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == TRUE &
          CIGAR_1_L_last_element == "S" &
          CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_second_element == "M" ~ POS_1,
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "S" &
          CIGAR_1_L_last_element == "M" &
          CIGAR_1_L_second_element == "M" ~ as.integer(POS_1 + READ_SPAN_MINUS_I -
                                                            CIGAR_1_N_first_element - 1),
        
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "M" ~ as.integer(as.integer(SATAG_1_4) + SATAG_1_6_N_last_element -1),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "M" ~ as.integer(as.integer(SATAG_1_4) + SATAG_1_6_N_last_element -1),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        TRUE ~ ERROR_NUMBER
      )
    ) %>%
    
    #calculate the position of the B flank on its chromosome/plasmid. THis is the position closest to the junction. This is not the final start position, if on the focus contig, this may been adjusting.

    mutate(
      FLANK_B_START_POS = case_when(
        CIGAR_1_LEN == 1 & FLANK_A_ORIENT=="FW" ~ as.integer(FlankAUltEnd+1),
        CIGAR_1_LEN == 1 & FLANK_A_ORIENT=="RV" ~ as.integer(FlankAUltEnd-1),
        CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "M" &
          SEQ_RCed_1 == TRUE ~ as.integer(POS_1 + CIGAR_1_N_first_element -1),
        CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "M" &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - CIGAR_1_N_last_element),
        
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == TRUE &
          CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "S" ~ as.integer(POS_1 + CIGAR_1_N_first_element -
                                                   1),
        SEQ_RCed_1 == TRUE &
          CIGAR_1_L_first_element == "S" &
          CIGAR_1_L_last_element == "M" & SATAG_1_1 != "SA" ~ NF_NUMBER,
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "S" & CIGAR_1_L_last_element == "M" ~ POS_1,
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "S" &
          CIGAR_1_L_last_element == "M" &
          CIGAR_1_L_second_element == "M" ~ POS_1 + CIGAR_1_N_second_element,
        SEQ_RCed_1 == FALSE &
          CIGAR_1_L_first_element == "M" &
          CIGAR_1_L_last_element == "S"  & SATAG_1_1 != "SA" ~ NF_NUMBER,
        
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "M" ~ as.integer(as.integer(SATAG_1_4) + SATAG_1_6_N_first_element -
                                                     1),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_first_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "+" &
          SATAG_1_6_L_last_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "M" ~ as.integer(as.integer(SATAG_1_4) + SATAG_1_6_N_first_element -
                                                     1),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          CIGAR_1_L_last_element == "S"  &
          SATAG_1_5 == "-" &
          SATAG_1_6_L_first_element == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        TRUE ~ ERROR_NUMBER
      )
    ) %>%

    rowwise()%>%
     mutate(FLANK_B_ORIENT = if (FLANK_B_START_POS != ERROR_NUMBER & FLANK_B_START_POS != NF_NUMBER & FLANK_B_END_POS != ERROR_NUMBER & FLANK_B_END_POS != NF_NUMBER){
      if (FLANK_B_START_POS < FLANK_B_END_POS){
        "FW"
      }else{
        "RV"
      }
    }else{
      "NOT_FOUND"
    }) %>%
    ungroup()

  
  function_time("Step 5 took ")
  
  ###############################################################################
  #Process data: step 6
  ###############################################################################
  
     data_improved4 = data_improved3 %>%
    #Determine the chromosome the B flank of the mate is mapped to
    mutate(
      MATE_FLANK_B_CHROM = (case_when(
        CIGAR_2_LEN == 1 ~ RNAME_2,
        (CIGAR_2_LEN > 1 &
           SEQ_RCed_2 == TRUE & CIGAR_2_L_first_element == "M") ~ RNAME_2,
        (CIGAR_2_LEN > 1 &
           SEQ_RCed_2 == FALSE & CIGAR_2_L_last_element == "M") ~ RNAME_2,
        (
          CIGAR_2_LEN > 1 &
            SEQ_RCed_2 == TRUE &
            SATAG_2_1 != "SA" & CIGAR_2_L_first_element == "S"
        ) ~ RNAME_2,
        (
          CIGAR_2_LEN > 1 &
            SEQ_RCed_2 == FALSE &
            SATAG_2_1 != "SA" & CIGAR_2_L_last_element == "S"
        ) ~ RNAME_2,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" & SATAG_2_5 == "+"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" & SATAG_2_5 == "-"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_first_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_last_element == "M" &
            SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" & SATAG_2_5 == "+"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" & SATAG_2_5 == "-"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_last_element == "M" &
            SATAG_2_5 == "+" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_first_element == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            CIGAR_2_L_last_element == "S" &
            SATAG_2_6_L_first_element == "M" &
            SATAG_2_5 == "-" &
            SATAG_2_12_L_last_element == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        TRUE ~ "ERROR"
      ))
    ) %>%
    #test whether B flanks from read and mate are mapped to same chromosome
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM,
                                              TRUE,
                                              FALSE))

  function_time("Step 6 took ")
  
  ###############################################################################
  #Process data: step 7
  ###############################################################################
  
  #remove rows where the chromosome of flank B as determined with read1 does not agree with read2
  if(DEBUG==FALSE){
    data_improved5 = data_improved4 %>% filter(MATE_FLANK_B_CHROM_AGREE == TRUE)
  }else{
    data_improved5 = data_improved4
  }

  data_improved5b = data_improved5 %>%
    
    #Search SEQ_1 for a sequence surrounding the DSB (meaning the cut has not been made or repaired perfectly)
    #search with allowing 1bp mismatch, but give different output whether the match is perfect or not
    rowwise() %>%
    mutate(DSB_AREA_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1))) %>%
    mutate(DSB_AREA_COUNT = length(DSB_AREA_CHECK@ranges)) %>%
    mutate(DSB_AREA_HIT = if(DSB_AREA_COUNT>0){
      if (SEQ_1_LEN >= (DSB_AREA_CHECK@ranges@start+DSB_AREA_CHECK@ranges@width-1)){
        as.character(DSB_AREA_CHECK[[1]])}else{""}
    }else{""} ) %>%
    ungroup()%>%
    mutate(DSB_AREA_INTACT = if_else(DSB_AREA_HIT == DSB_AREA_SEQ,
                                     "TRUE",
                                     "FALSE"))%>%
    mutate(DSB_AREA_1MM = if_else(
      DSB_AREA_INTACT == FALSE & DSB_AREA_COUNT>0,
      "TRUE",
      "FALSE")) %>%
    mutate(CASE_WT = if_else((
      DSB_AREA_INTACT==TRUE | DSB_AREA_1MM==TRUE),
      TRUE,
      FALSE)) %>%
  
    #select only columns that are used from now to save space
    select(
      PRIMER_SEQ,
      SEQ_1,
      DSB_AREA_INTACT,
      DSB_AREA_1MM,
      CASE_WT,
      MATE_FLANK_B_CHROM,
      FLANK_B_START_POS,
      FLANK_B_END_POS,
      FLANK_B_ORIENT,
      FLANK_B_CHROM,
      FILE_NAME,
      Max_BaseQual_1,
      MATE_FLANK_B_CHROM_AGREE,
      SEQ_2_first,
      QNAME_first,
      SEQ_1_LEN,
      FLANK_A_START_POS,
      hasPROBLEM,
      DEDUP_METHOD,
      TRIM_LEN,
      countEvents_init
    ) %>%
    
    
    mutate(FLANK_A_REF_LEN = as.integer(nchar(FLANK_A_REF))) %>%
    
    
    #Find how much SEQ_1 matches with FLANK_A_REF. Allow 1 bp mismatch somewhere, if the alignment after that continues for at least another 10 bp.
    rowwise() %>%
    mutate(FLANK_A_MATCH = matcher_skipper(FLANK_A_REF, SEQ_1)) %>%
    ungroup() %>%
    mutate(FLANK_A_LEN = nchar(FLANK_A_MATCH)) %>%
    mutate(FLANK_A_END_POS = case_when(FLANK_A_ORIENT == "FW" ~ as.integer(FLANK_A_START_POS + (FLANK_A_LEN -1)),
                                       FLANK_A_ORIENT == "RV" ~ as.integer(FLANK_A_START_POS - (FLANK_A_LEN -1)),
                                       TRUE ~ ERROR_NUMBER)) %>%
    
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        CASE_WT == TRUE ~ as.integer(0),
        CASE_WT != TRUE & FLANK_A_LEN != ERROR_NUMBER & FLANK_A_LEN != NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
        TRUE ~ ERROR_NUMBER
      )
    ) 
  
  function_time("Step 7 took ")
  
  ###############################################################################
  #Process data: step 8
  ###############################################################################
  
  data_improved6 = data_improved5b %>%
    mutate(TOTAL_REF = GLOBAL_TOTAL_REF) %>%
    #FLANK_B_REF. This ref includes homology.
    rowwise() %>%
    mutate(FLANK_B_REF =
             if (FLANK_B_CHROM == "NOT_FOUND" | FLANK_B_START_POS == ERROR_NUMBER | FLANK_B_END_POS == ERROR_NUMBER | CASE_WT == TRUE) {
               "NA"
             }else if (FLANK_B_CHROM == FOCUS_CONTIG & (abs(FlankAUltEnd - FLANK_B_END_POS) < MAX_DIST_FLANK_B_END) & FLANK_A_ORIENT == FLANK_B_ORIENT){
               if (FLANK_B_ORIENT == "FW" & FLANK_A_START_POS < FLANK_B_START_POS & FLANK_A_END_POS < FLANK_B_END_POS) {
                 substr(TOTAL_REF, start = 1, stop = (FLANK_B_END_POS-(FLANK_A_START_POS-1)))
               }else if (FLANK_B_ORIENT == "RV" & FLANK_A_START_POS > FLANK_B_START_POS & FLANK_A_END_POS > FLANK_B_END_POS){
                 substr(TOTAL_REF, start = 1, stop = (FLANK_A_START_POS-(FLANK_B_END_POS-1)))
               }else{
                 "NA"
               }
             }else{
               "NA"
             }) %>%
    mutate(FLANK_B_MATCH = if (FLANK_B_REF != "NA"){stri_reverse(matcher_skipper(stri_reverse(FLANK_B_REF), stri_reverse(SEQ_1)))
    }else{""}) %>%
    ungroup() %>%
    mutate(FLANK_B_MATCH_LEN = nchar(FLANK_B_MATCH)) %>% 

    mutate(NO_MATCH = case_when(FLANK_B_MATCH_LEN == 0 & FLANK_B_REF != "NA" ~ TRUE,
                                TRUE ~ FALSE)) %>% 
    

    mutate(
      hasPROBLEM = if_else(NO_MATCH == FALSE,
                           as.logical(hasPROBLEM),
                           TRUE)) %>%
    
    #flank b start position including MH
    mutate(FLANK_B_START_POS_MH = case_when(FLANK_B_REF != "NA" & FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_END_POS-(FLANK_B_MATCH_LEN-1)),
                                            FLANK_B_REF != "NA" & FLANK_B_ORIENT == "RV" ~ as.integer(FLANK_B_END_POS+(FLANK_B_MATCH_LEN-1)),
                                            TRUE ~ FLANK_B_START_POS)) %>%

    rowwise() %>%
    mutate(FLANK_B_LEN_MH = if_else(FLANK_B_ORIENT == "FW",
                                    FLANK_B_END_POS-(FLANK_B_START_POS_MH-1),
                                    FLANK_B_START_POS_MH-(FLANK_B_END_POS-1))) %>%
    
    #Extract the MH sequence
    mutate(MH = if (CASE_WT == FALSE & FLANK_B_CHROM != "NOT_FOUND" & FLANK_B_CHROM != "ERROR") {
      if (FLANK_B_CHROM == FOCUS_CONTIG & (abs(FlankAUltEnd - FLANK_B_END_POS) < MAX_DIST_FLANK_B_END) & FLANK_A_ORIENT == FLANK_B_ORIENT & (FLANK_A_LEN > (SEQ_1_LEN - FLANK_B_LEN_MH))) {
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
    } else{
      #if wt or flank b not identified
      ""
    }) %>%   
    #flank B start pos excluding MH for del calculation
    mutate(FLANK_B_START_POS_DEL = if (CASE_WT == TRUE){
      if (FLANK_A_ORIENT=="FW"){
        FLANK_A_END_POS+1
      }else{
        FLANK_A_END_POS-1
      }
    }else{
      if (FLANK_B_REF != "NA"){#if focus locus, right orientation
        if (nchar(MH) == 0){
          if (FLANK_B_ORIENT == "FW" & FLANK_A_END_POS >= FLANK_B_START_POS_MH){#case of insertion 
            FLANK_A_END_POS+1
          }else if (FLANK_B_ORIENT == "RV" & FLANK_A_END_POS <= FLANK_B_START_POS_MH){#case of insertion  
            FLANK_A_END_POS-1
          }else{#in case of deletion
            FLANK_B_START_POS_MH
          }
        }else{#if there is MH
          if (FLANK_B_ORIENT == "FW"){
            FLANK_B_START_POS_MH + nchar(MH)
          }else{#if RV
            FLANK_B_START_POS_MH - nchar(MH)
          }
        }
      }else{#in non focus locus cases or weird duplications
        FLANK_B_START_POS
      }})   %>%
    mutate(FLANK_B_LEN_DEL = if_else(FLANK_B_ORIENT == "FW",
                                     FLANK_B_END_POS-(FLANK_B_START_POS_DEL-1),
                                     FLANK_B_START_POS_DEL-(FLANK_B_END_POS-1))) %>%
    ungroup()
  
  function_time("Step 8 took ")
  
  ###############################################################################
  #Process data: step 9
  ###############################################################################
  
  data_improved7 = data_improved6 %>%
    mutate(FLANK_B_DEL = case_when(CASE_WT == TRUE ~ as.integer(0),
                                 CASE_WT != TRUE & FLANK_B_REF !="NA" & FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_START_POS_DEL - (FlankAUltEnd+1)),
                                   CASE_WT != TRUE & FLANK_B_REF !="NA" & FLANK_B_ORIENT == "RV" ~ as.integer((FlankAUltEnd-1) - FLANK_B_START_POS_DEL),
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

    #determine the filler sequence 
    mutate(
      FILLER = if_else(
        CASE_WT == FALSE & FLANK_B_CHROM != "NOT_FOUND" & FLANK_B_CHROM != "ERROR",					 
        substr(SEQ_1_WO_A, start = 1, stop = nchar(SEQ_1_WO_A) - FLANK_B_LEN_DEL),
        ""
      )) %>%
    mutate(insSize = nchar(FILLER)) %>%
    
    mutate(homologyLength = if_else(insSize == 0,
                                    nchar(MH),
                                    as.integer(-1))) %>%  
    
    #mark cases where the FLANK_B_LEN is too short 
    mutate(
      hasPROBLEM = if_else(FLANK_B_LEN_DEL >= FLANK_B_LEN_MIN,
                           as.logical(hasPROBLEM),
                           TRUE))  %>%
    
    mutate(SEQ_1_A = substr(SEQ_1, 1, FLANK_A_LEN)) %>%
    mutate(SEQ_1_B = substr(SEQ_1, SEQ_1_LEN - (FLANK_B_MATCH_LEN -1), SEQ_1_LEN)) %>%
    
    mutate(mismatch_found = case_when(SEQ_1_A != FLANK_A_MATCH ~ TRUE,
                                      SEQ_1_B != FLANK_B_MATCH ~ TRUE,
                                      TRUE ~ FALSE)) 
    
  function_time("Step 9 took ")
  
  ###############################################################################
  #Process data: step 10
  ###############################################################################
  
  data_improved8 = data_improved7 %>%
  
    #then apply a fix when the DSB is not set at 0, but when the total deletion length is 0. Also change the names for SIQPlotter. Also set the deletion length to 0 if the deletion length is negative, but the DSB is not allowed to occur downstream
    mutate(
      delRelativeStart = case_when(delSize == 0 & insSize == 0 ~ as.integer(0),
                                   (delSize !=0 | insSize != 0) & FLANK_A_ORIENT == "FW" & FLANK_A_DEL<0 ~ as.integer(0),
                                   TRUE ~ as.integer(FLANK_A_DEL * -1)),
      delRelativeEnd = if_else(
        delSize == 0 & insSize == 0,
        as.integer(0),
        as.integer(FLANK_B_DEL)
      )
    ) %>%
    
    #Add a check for insertion size = deletion size. That may indicate an incorrect call due to mismatches. Throw away those junctions resembling wt sequences.
    
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
    }) %>%
    ungroup() %>%
    mutate(CASE_WT = if_else(FAKE_DELIN_CHECK == TRUE,
                             TRUE,
                             CASE_WT))%>%


    #fix a bunch of things because of the updated CASE_WT
    mutate(
      delRelativeStart = if_else(CASE_WT == TRUE,
                                 as.integer(0),
                                 delRelativeStart),
      delRelativeEnd = if_else(CASE_WT == TRUE,
                               as.integer(0),
                               delRelativeEnd),
      FILLER = if_else(CASE_WT == TRUE,
                       "",
                       FILLER),
      MH = if_else(CASE_WT == TRUE,
                   "",
                   MH),
      insSize = if_else(CASE_WT == TRUE,
                        as.integer(0),
                        insSize),
      delSize = if_else(CASE_WT == TRUE,
                        as.integer(0),
                        delSize),
      homologyLength = if_else(CASE_WT == TRUE,
                               as.integer(0),
                               homologyLength),
      FLANK_B_ISFORWARD = if_else(FLANK_B_START_POS < FLANK_B_END_POS,
                                  TRUE,
                                  FALSE)) %>%
    rowwise() %>%
    #fix the FLANK_B_START_POS again if the event is wt
    mutate(FLANK_B_START_POS = if (delSize == 0 & insSize == 0){
      if (FLANK_A_ORIENT == "FW"){
        FlankAUltEnd + 1
      }else{
        FlankAUltEnd - 1
      }
    }else{
      FLANK_B_START_POS_DEL
    }) %>%
    ungroup() %>%
    #recalculate the FLANK B del
    mutate(FLANK_B_LEN_DEL = case_when(FLANK_B_ORIENT == "FW" ~ (FLANK_B_END_POS - FLANK_B_START_POS),
                                     TRUE ~ FLANK_B_START_POS - FLANK_B_END_POS)) %>%
    
    #Filter short flank Bs again, with the newly calculated flank b len
    filter(FLANK_B_LEN_DEL >= FLANK_B_LEN_MIN) %>%
    
    #flag erroneous events, and filter these away
    mutate(
      hasPROBLEM = if_else(
        FLANK_B_CHROM == "ERROR" |
          MATE_FLANK_B_CHROM == "NOT_FOUND" |
          MATE_FLANK_B_CHROM == "ERROR" |
          MATE_FLANK_B_CHROM_AGREE == FALSE |
          (CASE_WT == TRUE & FLANK_B_CHROM != FOCUS_CONTIG),
        TRUE,
        as.logical(hasPROBLEM)
      )
    ) 
    #remove rows with problems
    if(DEBUG==FALSE){
      data_improved8b = data_improved8 %>% filter(hasPROBLEM == FALSE)
    }else{
      data_improved8b = data_improved8
    }
    
    data_improved8c=data_improved8b %>%
    #sort the rows in order to have the ones with the highest base quality at the top
    arrange(desc(SEQ_1_LEN), desc(Max_BaseQual_1)) %>%
    
    #calculate the number of events
    group_by(
      FILE_NAME,
      PRIMER_SEQ,
      delRelativeStart,
      delRelativeEnd,
      FLANK_B_CHROM,
      FLANK_B_START_POS,
      MATE_FLANK_B_CHROM,
      MATE_FLANK_B_CHROM_AGREE,
      FLANK_B_ISFORWARD,
      FILLER,
      MH,
      insSize,
      delSize,
      homologyLength,
      mismatch_found,
      FAKE_DELIN_CHECK,
      DSB_AREA_INTACT,
      DSB_AREA_1MM,
      hasPROBLEM,
      DEDUP_METHOD,
      TRIM_LEN
    ) %>%
    summarize(
      countEvents = sum(countEvents_init),
      Max_BaseQual_1_max = max(Max_BaseQual_1),
      Name = first(QNAME_first),
      SEQ_1_first = first(SEQ_1),
      SEQ_2_first_first = first(SEQ_2_first)
    ) %>%
    
    #Add a subject and type column. Also add two extra columns for SIQplotter that are equal to the other del columns.
    mutate(
      Subject = FOCUS_LOCUS,
      Type = case_when(
        delSize == 0 & insSize == 0 & mismatch_found == FALSE & FAKE_DELIN_CHECK == FALSE & DSB_AREA_1MM==FALSE ~ "WT",
        (delSize == 0 & insSize == 0 & mismatch_found == TRUE) | FAKE_DELIN_CHECK==TRUE | DSB_AREA_1MM==TRUE | (DSB_AREA_INTACT == TRUE & FLANK_A_ORIENT=="FW" & FLANK_B_START_POS != FlankAUltEnd+1) | (DSB_AREA_INTACT == TRUE & FLANK_A_ORIENT=="RV" & FLANK_B_START_POS != FlankAUltEnd-1)  ~ "SNV",
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
      countEvents,
      FLANK_B_CHROM,
      FLANK_B_START_POS,
      FLANK_B_ISFORWARD,
      MATE_FLANK_B_CHROM,
      MATE_FLANK_B_CHROM_AGREE,
      mismatch_found,
      FAKE_DELIN_CHECK,
      FILLER,
      MH,
      insSize,
      homologyLength,
      delSize,
      Type,
      Subject,
      Max_BaseQual_1_max,
      SEQ_1_first,
      SEQ_2_first_first,
      DSB_AREA_INTACT,
      DSB_AREA_1MM,
      hasPROBLEM,
      DEDUP_METHOD,
      TRIM_LEN
    ) %>%
    
    #add required columns for SIQplotter
    mutate(getHomologyColor= "grey")
  
  function_time("Step 10 took ")
  
  ###############################################################################
  #Process data: step 11
  ###############################################################################
  
  #calculate the fraction of reads with the same outcomes
  data_improved9 = data_improved8c %>% group_by(FILE_NAME, PRIMER_SEQ) %>% summarize(SumCountEvents =
                                                                                  sum(countEvents))
  function_time("Step 11 took ")
  
  ###############################################################################
  #Process data: step 12
  ###############################################################################
  
  data_improved10 = left_join(data_improved8c, data_improved9, by = c("FILE_NAME", "PRIMER_SEQ")) %>%
    mutate(fraction = countEvents / SumCountEvents) %>%
    mutate(SumCountEvents = NULL,
           FlankAUltEnd = FlankAUltEnd,
           FLANK_A_ORIENT = FLANK_A_ORIENT,
           FOCUS_CONTIG = FOCUS_CONTIG,
           Genotype = Genotype,
           DNASample = DNASample,
           Ecotype = Ecotype,
           RunID = RunID,
           program_version = hash,
           Alias = paste0(Library, "_", RunID))
  
  #filter out duplicate position artefacts
  data_improved11 = data_improved10 %>% group_by(FLANK_B_CHROM, FLANK_B_START_POS) %>% summarize(countEventsPosSum = sum(countEvents))
  data_improved12 = data_improved10 %>% group_by(FLANK_B_CHROM, FLANK_B_START_POS) %>% summarize(countEventsPosMax = max(countEvents))
  data_improved13 = left_join(data_improved11, data_improved12, by = c("FLANK_B_CHROM", "FLANK_B_START_POS")) %>% mutate(MajorFractionAtPos = countEventsPosMax / countEventsPosSum) 
  data_improved14 = left_join(data_improved10, data_improved13, by = c("FLANK_B_CHROM", "FLANK_B_START_POS")) %>%mutate(FractionAtPos=countEvents/countEventsPosSum) 
  
  
  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, data_improved10)
  saveWorkbook(work_book, file = paste0(output_dir, Sample, "_", RunID, "_CISGUIDE_V_", hash_little, ".xlsx"), overwrite = TRUE)
  
  function_time("Step 12 took ")
}

sample_list = list.files(path=output_dir, pattern = "\\.xlsx")
wb = tibble()

for (i in sample_list){
  wb=rbind(wb, read.xlsx(paste0(output_dir, i)))
}
work_book2 <- createWorkbook()
addWorksheet(work_book2, "rawData")
writeData(work_book2, sheet = 1, wb)
saveWorkbook(work_book2, file = paste0(output_dir, "CISGUIDE_V_", hash_little, ".xlsx"), overwrite = TRUE)

#Write an additional sheet with read number info
read_numbers_info = read.csv(paste0(input_dir, "read_numbers.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
wb_numbers = read_numbers_info %>% 
  mutate(Alias = paste0(Sample, "_", RunID))%>%
  mutate(Sample = NULL,
         RunID = NULL)

addWorksheet(work_book2, "Information")
writeData(work_book2, sheet = 2, wb_numbers)
saveWorkbook(work_book2, file = paste0(output_dir, "CISGUIDE_V_", hash_little, ".xlsx"), overwrite = TRUE)

