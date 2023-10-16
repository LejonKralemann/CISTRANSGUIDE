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
MINBASEQUAL = 30 #minimum base quality (phred)
MAX_DIST_FLANK_B_END = 10000 #distance from end of flank B to DSB, determines max deletion size and also affects maximum insertion size
FLANK_B_LEN_MIN = 15 #minimum length of flank B. Also affects size of DSB_AREA_SEQ (2x FLANK_B_LEN_MIN)
LOCUS_WINDOW = 1000 #size of the window centered on the DSB, RB nick, or LB nick to determine locus info
DEBUG=FALSE #if true, reads will not be discarded when a problem has been detected, but flagged.
CISSWITCH=FALSE #if true, one step will be skipped to allow the analysis of different events at the same genomic location.

###############################################################################
#set parameters - non-adjustable
###############################################################################
NF_NUMBER = as.integer(-99999999) #don't change
ERROR_NUMBER = as.integer(99999999) #don't change
hash=system("git rev-parse HEAD", intern=TRUE)
hash_little=substr(hash, 1, 8)
sample_info = read.csv(paste0(input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
TIME_START=as.numeric(Sys.time())*1000

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
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
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
  REF = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ref))
  FlankAUltEnd = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FlankAUltEnd))
  FlankBUltStart = as.integer(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(FlankBUltStart))
  DNASample = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(DNA))
  Ecotype = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Ecotype))
  Library = as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Sample))
  
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
  Primer_seq = str_replace_all(as.character(sample_info %>% filter(row.names(sample_info) %in% i) %>% select(Primer)), "TCAGACGTGTGCTCTTCCGATCT", "")
  Primer_seq_len = nchar(Primer_seq)
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
  if (Primer_seq != "NA"){
  FASTA_MODE = FALSE
  }else{
  FASTA_MODE = TRUE
  }
  
  if (FASTA_MODE == FALSE){
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
  
  #set the minimum length of a read
  if (FASTA_MODE == FALSE){
  MINLEN = PRIMER_TO_DSB_GLOBAL+FLANK_B_LEN_MIN
  }else{
  MINLEN = 60
  message("MINLEN set to 60 because no primer seq available")
  }
  
  #the following things need to be obtained for fasta mode: PRIMER_TO_DSB (step7), FLANK_A_REF (step7), GLOBAL_TOTAL_REF (step8)
  
  function_time("Step 1 took ")

  ###############################################################################
  #Process data: step 2
  ###############################################################################
  
  data_improved_a  = data %>%
    
    #Count number of Ns and remove any reads with Ns
    mutate(NrN = str_count(SEQ_1, pattern = "N"),
           SEQ_1_LEN = nchar(SEQ_1)) %>%
    filter(NrN < 1) %>%
    
    #filter away reads that are too short
    filter(!(SEQ_1_LEN < MINLEN)) %>%
    #calculate the average base quality
    rowwise() %>%
    mutate(AvgBaseQual_1 = mean(utf8ToInt(QUAL_1)-33)) %>%
    mutate(AvgBaseQual_2 = mean(utf8ToInt(QUAL_2)-33)) %>%
    ungroup()
  
  data_improved_b = data_improved_a %>%
    
    
    #calculate the end position of the B flank in the mate. This is the end the furthest away from the junction.
    mutate(MATE_B_END_POS = case_when(MATE_B_ORIENT=="FW" ~ as.integer(MATE_B_POS),
                                      MATE_B_ORIENT=="RV" ~ as.integer(MATE_B_POS + 29),
                                      TRUE ~ ERROR_NUMBER))  %>%
    
    #test whether B flanks from read and mate are mapped to same chromosome, and agree on orientation.
    #note that the two reads are in two different orientations, so the B orientation should not agree
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM & FLANK_B_ORIENT != MATE_B_ORIENT,
                                              TRUE,
                                              FALSE)) 
    
    
    #remove rows where the chromosome of flank B as determined with read1 does not agree with read2
    if(DEBUG==FALSE){
      data_improvedc = data_improved_b %>% filter(MATE_FLANK_B_CHROM_AGREE == TRUE)
    }else{
      data_improvedc = data_improved_b
    }
    
  data_improved1 = data_improvedc %>%
  
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
      SEQ_1_LEN) %>%
    arrange(desc(AvgBaseQual_1), desc(AvgBaseQual_2)) %>%
    summarize(QNAME_first = dplyr::first(QNAME),
              SEQ_2_first = dplyr::first(SEQ_2),
              MATE_B_END_POS_list = toString(MATE_B_END_POS)
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
    #Search SEQ_1 for a sequence surrounding the DSB (meaning the cut has not been made or has been repaired perfectly)
    #search with allowing 1bp mismatch, but give different output whether the match is perfect or not
    #note it only checks the first hit
    rowwise() %>%
    mutate(DSB_AREA_CHECK = list(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1))) %>%
    mutate(DSB_AREA_COUNT = length(DSB_AREA_CHECK@ranges)) %>%
    mutate(DSB_AREA_HIT = if(DSB_AREA_COUNT>0){
      if (SEQ_1_LEN >= (DSB_AREA_CHECK@ranges@start[1]+DSB_AREA_CHECK@ranges@width[1]-1)){
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
    mutate(DSB_HIT_MULTI = if_else(DSB_AREA_COUNT>1,
                                   "TRUE",
                                   "FALSE")) %>%
  
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
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        CASE_WT == TRUE ~ as.integer(0),
        CASE_WT != TRUE & FLANK_A_LEN != ERROR_NUMBER & FLANK_A_LEN != NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
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
             if (CASE_WT == TRUE) {
               "NA"
             }else if (FLANK_B_CHROM == FOCUS_CONTIG & 
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
    mutate(FLANK_B_MATCH = if (FLANK_B_REF != "NA"){stri_reverse(matcher_skipper(stri_reverse(FLANK_B_REF), stri_reverse(SEQ_1)))
    }else{""}) %>%
    ungroup() %>%
    mutate(FLANK_B_MATCH_LEN = nchar(FLANK_B_MATCH)) %>% 

    mutate(NO_MATCH = case_when(FLANK_B_MATCH_LEN == 0 & FLANK_B_REF != "NA" ~ TRUE,
                                TRUE ~ FALSE))  %>%
    mutate(hasPROBLEM = if_else(NO_MATCH == FALSE,
                       FALSE,
                       TRUE)) 

  
    data_improved4 = data_improved3b %>%
      
    #flank b start position including MH
    mutate(FLANK_B_START_POS_MH = case_when(FLANK_B_REF != "NA" & FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_END_POS-(FLANK_B_MATCH_LEN-1)),
                                            FLANK_B_REF != "NA" & FLANK_B_ORIENT == "RV" ~ as.integer(FLANK_B_END_POS+(FLANK_B_MATCH_LEN-1)),
                                            CASE_WT == TRUE & FLANK_A_ORIENT == "FW" ~ as.integer(FlankAUltEnd+1), 
                                            CASE_WT == TRUE & FLANK_A_ORIENT == "RV" ~ as.integer(FlankAUltEnd-1), 
                                            TRUE ~ ERROR_NUMBER)) %>%

    rowwise() %>%
    mutate(FLANK_B_LEN_MH = if_else(FLANK_B_ORIENT == "FW",
                                    FLANK_B_END_POS-(FLANK_B_START_POS_MH-1),
                                    FLANK_B_START_POS_MH-(FLANK_B_END_POS-1))) %>%
    
    #Extract the MH sequence
    mutate(MH = if (CASE_WT == FALSE & FLANK_B_CHROM != "NOT_FOUND" & FLANK_B_CHROM != "ERROR") {
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
    } else{
      #if wt or flank b not identified
      ""
    }) %>%   
    #flank B start pos excluding MH (for del calculation)
      mutate(FLANK_B_START_POS_DEL = if (CASE_WT == TRUE){
        if (FLANK_A_ORIENT=="FW"){
          FLANK_A_END_POS+1
        }else{
          FLANK_A_END_POS-1
        }
      }else{
        if (FLANK_B_ORIENT == "FW"){
          FLANK_B_START_POS_MH + nchar(MH)
        }else{#if RV
          FLANK_B_START_POS_MH - nchar(MH)
        }})%>%
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
    mutate(FLANK_B_DEL = case_when(CASE_WT == TRUE ~ as.integer(0),
                                   CASE_WT == FALSE & 
                                     FLANK_B_CHROM == FOCUS_CONTIG & 
                                     FLANK_B_CLOSE_BY==TRUE & 
                                     FLANK_A_ORIENT == FLANK_B_ORIENT & 
                                     FLANK_B_ORIENT == "FW" ~ as.integer(FLANK_B_START_POS_DEL - (FlankAUltEnd+1)),
                                   CASE_WT == FALSE & 
                                     FLANK_B_CHROM == FOCUS_CONTIG & 
                                     FLANK_B_CLOSE_BY==TRUE & 
                                     FLANK_A_ORIENT == FLANK_B_ORIENT & 
                                     FLANK_B_ORIENT == "RV" ~ as.integer((FlankAUltEnd-1) - FLANK_B_START_POS_DEL),
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
    
  function_time("Step 6 took ")
  
  ###############################################################################
  #Process data: step 7
  ###############################################################################
  
  data_improved6 = data_improved5 %>%
  
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
      FLANK_B_ISFORWARD = if_else(FLANK_B_ORIENT=="FW",
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
      FLANK_B_START_POS_MH
    }) %>%
    
      #flag erroneous events, and filter these away
    mutate(
      hasPROBLEM = if_else(
      CASE_WT == TRUE & FLANK_B_CHROM != FOCUS_CONTIG,
        TRUE,
        as.logical(hasPROBLEM)
      )
    )%>%
    
    #calculate the minimum length of the read to capture the entire outcome
    #then trim the reads to that minimum.
    #any that have less than this minimum get removed
    mutate(read_minimum_length = FLANK_A_LEN + insSize + 30)%>%
    filter(SEQ_1_LEN >= read_minimum_length)%>%
    mutate(SEQ_1_trimmed = substr(SEQ_1, 1, read_minimum_length))
  
  
  
    #remove rows with problems
    if(DEBUG==FALSE){
      data_improved7 = data_improved6 %>% filter(hasPROBLEM == FALSE)
    }else{
      data_improved7 = data_improved6
    }
    
  #here grouping will occur to determine the number of anchors.
  #for TRANSGUIDE a consensus outcome will be determined
         if (CISSWITCH == TRUE){
      
        data_improved8pre2 =data_improved7 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
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
        MATE_B_END_POS_max = max(as.integer(MATE_B_END_POS_list)),
        MATE_B_END_POS_min = min(as.integer(MATE_B_END_POS_list))
      ) 
      }else{
        
        data_improved8pre2 =data_improved7 %>%
          ungroup()%>%
          separate_longer_delim(cols="MATE_B_END_POS_list", delim = ",") %>%
          group_by(
            FILE_NAME,
            PRIMER_SEQ,
            FLANK_B_CHROM,
            FLANK_B_START_POS,
            FLANK_B_ISFORWARD,
            FLANK_B_ORIENT,
            TRIM_LEN
          )%>%
          summarize(
            ReadCount = n(),
            Name = dplyr::first(QNAME_first),
            Count_consensus = max(table(SEQ_1_trimmed)),
            AnchorCount = n_distinct(MATE_B_END_POS_list),
            SEQ_1_con = names(which.max(table(SEQ_1_trimmed))),
            SEQ_2_con = names(which.max(table(SEQ_2_first))),
            MATE_B_END_POS_max = max(as.integer(MATE_B_END_POS_list)),
            MATE_B_END_POS_min = min(as.integer(MATE_B_END_POS_list)),
            delRelativeStart_con = names(which.max(table(delRelativeStart))),
            delRelativeEnd_con = names(which.max(table(delRelativeEnd))),
            FILLER_con = names(which.max(table(FILLER))),
            MH_con =names(which.max(table(MH))),
            insSize_con = names(which.max(table(insSize))),
            delSize_con = names(which.max(table(delSize))),
            homologyLength_con = names(which.max(table(homologyLength))),
            FAKE_DELIN_CHECK_con = names(which.max(table(FAKE_DELIN_CHECK))),
            DSB_AREA_INTACT_con = names(which.max(table(DSB_AREA_INTACT))),
            DSB_AREA_1MM_con = names(which.max(table(DSB_AREA_1MM))),
            DSB_HIT_MULTI_con = names(which.max(table(DSB_HIT_MULTI)))
          )%>%
          mutate(Consensus_freq = Count_consensus/ReadCount)%>%
          #rename columns to that code below is compatible with both CISSWITCH==TRUE and CISSWITCH==FALSE
          rename(delRelativeStart = delRelativeStart_con,
                 delRelativeEnd = delRelativeEnd_con,
                 FILLER = FILLER_con,
                 MH = MH_con,
                 insSize = insSize_con,
                 delSize = delSize_con,
                 homologyLength = homologyLength_con,
                 FAKE_DELIN_CHECK = FAKE_DELIN_CHECK_con,
                 DSB_AREA_INTACT = DSB_AREA_INTACT_con,
                 DSB_AREA_1MM = DSB_AREA_1MM_con,
                 DSB_HIT_MULTI = DSB_HIT_MULTI_con
                 )%>%
          #flag junctions when the dominant SEQ_1 sequence is not at least 3/4 of total SEQ_1's at this junction
          mutate(hasProblem = if_else(Consensus_freq <0.75,
                                      TRUE,
                                      FALSE))
      }
        
    data_improved8 =data_improved8pre2 %>%
      
      #then calculate the difference between the end of the primary read and the end of the mate.
      mutate(ANCHOR_DIST = case_when(FLANK_B_ORIENT=="FW" ~ as.integer(MATE_B_END_POS_max - FLANK_B_START_POS),
                                     FLANK_B_ORIENT=="RV" ~ as.integer(FLANK_B_START_POS - MATE_B_END_POS_min),
                                     TRUE ~ ERROR_NUMBER)) %>%
      #flag as problematic
      #when anchors are not far away enough
      #when anchor count is too low
      mutate(hasProblem = case_when(AnchorCount < 3 ~ TRUE,
                                    ANCHOR_DIST < 150 ~ TRUE,
                                    TRUE ~ as.logical(hasProblem)))%>%
      
    #Add a subject and type column. Also add two extra columns for SIQplotter that are equal to the other del columns.
    mutate(
      Subject = FOCUS_LOCUS,
      Type = case_when(
        delSize == 0 & insSize == 0 & FAKE_DELIN_CHECK == FALSE & DSB_AREA_1MM==FALSE ~ "WT",
        FAKE_DELIN_CHECK==TRUE | DSB_AREA_1MM==TRUE | (DSB_AREA_INTACT == TRUE & FLANK_A_ORIENT=="FW" & FLANK_B_START_POS != FlankAUltEnd+1) | (DSB_AREA_INTACT == TRUE & FLANK_A_ORIENT=="RV" & FLANK_B_START_POS != FlankAUltEnd-1)  ~ "WT",
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
      FAKE_DELIN_CHECK,
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
      hasProblem
    ) %>%
    

    #add required columns for SIQplotter
    mutate(getHomologyColor= "grey")
  
  function_time("Step 7 took ")
  
  ###############################################################################
  #Process data: step 8
  ###############################################################################
  
  #calculate the fraction of reads with a certain outcome within a library
  data_improved9 = data_improved8 %>% group_by(FILE_NAME) %>% summarize(ReadCountTotal =
                                                                                  sum(ReadCount))
  function_time("Step 8 took ")
  
  ###############################################################################
  #Process data: step 9
  ###############################################################################
  
  data_improved10 = left_join(data_improved8, data_improved9, by = "FILE_NAME") %>%
    mutate(fraction = ReadCount / ReadCountTotal) %>%
    mutate(countReadsTotal = NULL,
           FlankAUltEnd = FlankAUltEnd,
           FLANK_A_ORIENT = FLANK_A_ORIENT,
           FOCUS_CONTIG = FOCUS_CONTIG,
           Genotype = Genotype,
           DNASample = DNASample,
           Ecotype = Ecotype,
           RunID = RunID,
           program_version = hash,
           Plasmid = PLASMID,
           Plasmid_alt = PLASMID_ALT,
           Alias = paste0(Library, "_", RunID))
  

  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, data_improved10)
  saveWorkbook(work_book, file = paste0(output_dir, Sample, "_", RunID, "_CISGUIDE_V_", hash_little, "_", as.integer(Sys.time()), ".xlsx"), overwrite = TRUE)
  
  function_time("Step 9 took ")
}

sample_list = list.files(path=output_dir, pattern = "\\.xlsx")
wb = tibble()

for (i in sample_list){
  wb=rbind(wb, read.xlsx(paste0(output_dir, i)))
}
work_book2 <- createWorkbook()
addWorksheet(work_book2, "rawData")
writeData(work_book2, sheet = 1, wb)

#Write an additional sheet with read number info
read_numbers_info = read.csv(paste0(input_dir, "read_numbers.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
wb_numbers = read_numbers_info %>% 
  mutate(Alias = paste0(Sample, "_", RunID))%>%
  mutate(Sample = NULL,
         RunID = NULL)

addWorksheet(work_book2, "Information")
writeData(work_book2, sheet = 2, wb_numbers)
saveWorkbook(work_book2, file = paste0(output_dir, "Data_combined_CISGUIDE_V_", hash_little, "_", as.integer(Sys.time()), ".xlsx"), overwrite = TRUE)

