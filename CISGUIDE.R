library(BSgenome)
library(Biostrings)
library(stringi)
library(stringdist)
library(tidyverse)
library(openxlsx)

#The data should be "pre-treated" by running the script "Reverse TRANSGUIDE primary script" in bash. This performs some filtering like removing supplemental reads and requiring reads to start with the primer sequence. 
#change the chromosome names of the fasta file to add Chr, and remove all the other crap, also change - to _
input_dir= "C:/Users/lejon/Documents/Scripts/CISGUIDE/input/"
output_dir= "C:/Users/lejon/Documents/Scripts/CISGUIDE/output/"
version_no="6.15"
NF_NUMBER = as.integer(-99999999) #don't change
ERROR_NUMBER = as.integer(99999999) #don't change
MINLEN = as.integer(90) #minimum length of a read
MINBASEQUAL = 0.75 #minimum base quality
MAX_DIST_FLANK_B_END = 10000 #distance from end of flank B to DSB, determines max deletion size and also affects maximum insertion size
FLANK_B_LEN_MIN = 15 #minimum length of flank B
LOCUS_WINDOW = 1000 #size of the window centered on the DSB, RB nick, or LB nick to determine locus info

#Reads information you provide about the samples
sample_info = read.csv(paste0(input_dir, "Sample_information.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)

#note that below the flank A is the flank that starts with the primer until the break, with optional deletion.
#flank B starts with the break, and ends wherever the read 1 sequence ends.

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


for (i in sample_info$Sample){
  
  if (file.exists(paste0(input_dir, i, "_A.txt"))==FALSE){
    message("Primary processed file not found")
    next
  }
  
  data = read.csv(paste0(input_dir, i, "_A.txt"), sep = "\t", header=T, stringsAsFactors = FALSE)
  FOCUS_CONTIG = str_replace_all(as.character(sample_info %>% filter(Sample==i) %>% select(DSB_chrom)), "-", "_")
  FOCUS_LOCUS = as.character(sample_info %>% filter(Sample==i) %>% select(Locus_name))
  Genotype = as.character(sample_info %>% filter(Sample==i) %>% select(Genotype))
  PLASMID = str_replace_all(as.character(sample_info %>% filter(Sample==i) %>% select(Plasmid)), "-", "_")
  FlankAUltEnd = as.integer(sample_info %>% filter(Sample==i) %>% select(FlankAUltEnd))
  FlankBUltStart = as.integer(sample_info %>% filter(Sample==i) %>% select(FlankBUltStart))
  genomeseq = readDNAStringSet(paste0(input_dir, str_replace_all(PLASMID,"_", "-"),".fa") , format="fasta")
  contig_seq = as.character(eval(parse(text = paste0("genomeseq$", FOCUS_CONTIG))))
  Primer_seq = str_replace_all(as.character(sample_info %>% filter(Sample==i) %>% select(Primer)), "TCAGACGTGTGCTCTTCCGATCT", "")
  Primer_match = as.data.frame(matchPattern(pattern = Primer_seq, subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  Primer_RC_match = as.data.frame(matchPattern(pattern = as.character(reverseComplement(DNAString(Primer_seq))), subject = DNAString(contig_seq), max.mismatch = 0, fixed=TRUE))
  FLANK_A_ORIENT = as.character(sample_info %>% filter(Sample==i) %>% select(FLANK_A_ORIENT))
  PrimerType = as.character(sample_info %>% filter(Sample==i) %>% select(PrimerType))

  
  DSB_AREA_SEQ = (if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= FlankAUltEnd - 14, stop= FlankAUltEnd + 15)
  }else if (FLANK_A_ORIENT == "RV"){
    as.character(reverseComplement(DNAString(substr(contig_seq, start= FlankAUltEnd - 15, stop= FlankAUltEnd + 14))))
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
  
  PRIMER_TO_DSB = if (FLANK_A_ORIENT == "FW"){
    FlankAUltEnd - (Primer_pos -1)
  }else if (FLANK_A_ORIENT=="RV"){
    Primer_pos - (FlankAUltEnd -1)
  }else{
    ERROR_NUMBER
  }
  
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
  TOTAL_REF = if (FLANK_A_ORIENT == "FW"){
    substr(contig_seq, start= Primer_pos, stop= Primer_pos+MAX_DIST_FLANK_B_END+PRIMER_TO_DSB)
  }else if (FLANK_A_ORIENT=="RV"){
    as.character(reverseComplement(DNAString(
      substr(contig_seq, start= Primer_pos-(MAX_DIST_FLANK_B_END+PRIMER_TO_DSB), stop= Primer_pos))
    ))
  }else{
    ""
  } 
  
  data_improved  = data %>%
    
    #filter(QNAME == "A00379:436:H3CHWDMXY:2:1110:10131:22138") %>%
    #remove reads with Ns
    rowwise() %>%
    mutate(NrN = str_count(SEQ_1, pattern = "N"),
           SEQ_1_LEN = nchar(SEQ_1)) %>%
    ungroup() %>%
    filter(NrN < 1) %>%
    
    #filter away reads that are too short
    filter(!(SEQ_1_LEN < MINLEN)) %>%
    
    #calculate the average base quality
    rowwise() %>%
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
    mutate(QNAME = str_replace_all(QNAME, "-", "_"))%>%
    ungroup()
  
  
  data_improved1 = data_improved %>%
    
    #select informative columns
    select(
      QNAME,
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
      SEQ_2,
      FILE_NAME,
      PRIMER_SEQ,
      AvgBaseQual_1,
      SEQ_1_LEN
    ) %>%
    rowwise() %>%
    mutate(RNAME_1 = as.character(RNAME_1))%>%
    mutate(RNAME_2 = as.character(RNAME_2))%>%
    ungroup() %>%
    
    #collapse optical duplicates (UMI collapsing was already done before mapping)
    
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
      SEQ_2,
      FILE_NAME,
      PRIMER_SEQ,
      SEQ_1_LEN
    ) %>%
    summarize(
      QNAME_1st = first(QNAME),
      NrOpticalDuplicates = n(),
      AvgBaseQual_1_Max = max(AvgBaseQual_1)
    ) %>%
    
    #filter away all reads under a certain quality threshold
    filter(AvgBaseQual_1_Max > MINBASEQUAL) %>%
    
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
    )
  
  
  data_improved2 = data_improved1 %>%
    
    #temp filters
    filter(SATAG_1_1 != "XA") %>%
    filter(SATAG_2_1 != "XA") %>%
    
    #calculate read span length. This is basically equal to read length, plus deletionlength within an alignment.
    rowwise() %>%
    mutate(READ_SPAN = sum(as.integer(CIGAR_1_N))) %>%
    
    #determine which chromosome the B flank is aligning to
    
    mutate(
      FLANK_B_CHROM = str_replace_all((case_when(
        (CIGAR_1_LEN == 1) ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == TRUE & head(CIGAR_1_L, 1) == "M") ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == TRUE  & head(CIGAR_1_L, 1) == "S" & SATAG_1_1 != "SA") ~ "NOT_FOUND",
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == FALSE & tail(CIGAR_1_L, 1) == "M") ~ RNAME_1,
        (CIGAR_1_LEN > 1 & SEQ_RCed_1 == FALSE  & tail(CIGAR_1_L, 1) == "S" & SATAG_1_1 != "SA") ~ "NOT_FOUND",
        (SEQ_RCed_1 == TRUE & SATAG_1_1 == "SA" & head(CIGAR_1_L, 1) == "S" & tail(SATAG_1_6_L, n = 1) == "M" & SATAG_1_5 == "+") ~ SATAG_1_3,
        (SEQ_RCed_1 == TRUE & SATAG_1_1 == "SA" & head(CIGAR_1_L, 1) == "S" & tail(SATAG_1_6_L, n = 1) == "S" & SATAG_1_5 == "+" & is.na(SATAG_1_12)) ~ "NOT_FOUND",
        (
          SEQ_RCed_1 == TRUE &
            SATAG_1_1 == "SA" &
            head(CIGAR_1_L, 1) == "S" &
            head(SATAG_1_6_L, n = 1) == "M" & SATAG_1_5 == "-"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == TRUE &
            SATAG_1_1 == "SA" &
            head(CIGAR_1_L, 1) == "S" &
            head(SATAG_1_6_L, n = 1) == "S" &
            SATAG_1_5 == "-" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            tail(CIGAR_1_L, 1) == "S" &
            tail(SATAG_1_6_L, n = 1) == "M" & SATAG_1_5 == "+"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            tail(CIGAR_1_L, 1) == "S" &
            tail(SATAG_1_6_L, n = 1) == "S" &
            SATAG_1_5 == "+" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            tail(CIGAR_1_L, 1) == "S" &
            head(SATAG_1_6_L, n = 1) == "M" & SATAG_1_5 == "-"
        ) ~ SATAG_1_3,
        (
          SEQ_RCed_1 == FALSE &
            SATAG_1_1 == "SA" &
            tail(CIGAR_1_L, 1) == "S" &
            head(SATAG_1_6_L, n = 1) == "S" &
            SATAG_1_5 == "-" & is.na(SATAG_1_12)
        ) ~ "NOT_FOUND",
        TRUE ~ "ERROR"
      )), "-", "_")
    ) %>%
    mutate(FLANK_B_CHROM =if_else(
      FLANK_B_CHROM != "Mt" & FLANK_B_CHROM != "Pt" & FLANK_B_CHROM != PLASMID & FLANK_B_CHROM != "NOT_FOUND",
      paste0("Chr", FLANK_B_CHROM),
      FLANK_B_CHROM)) %>%
    ungroup() %>%
    
    #filter away cases where flank B chrom is not found
    filter(FLANK_B_CHROM != "NOT_FOUND") %>%
    
    #write down the various potential insertions
    rowwise() %>%
    mutate(
      POS2_I = case_when(
        CIGAR_1_LEN > 2  &
          tail(head(CIGAR_1_L, 2), 1) == "I" ~ as.integer(tail(head(CIGAR_1_N, 2), 1)),
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
    mutate(READ_SPAN_MINUS_I = as.integer(READ_SPAN - (POS2_I + POS3_I + POS4_I + POS5_I + POS6_I + POS7_I)))
  
  data_improved3 = data_improved2 %>%
    
    #calculate the end position of the B flank, in order to get the ref seq. This is the end the furthest away from the junction.
    
    mutate(
      FLANK_B_END_POS = case_when(
        CIGAR_1_LEN == 1 & SEQ_RCed_1 == TRUE ~ POS_1,
        CIGAR_1_LEN == 1 &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - 1),
        head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "M" & SEQ_RCed_1 == TRUE ~ POS_1,
        head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "M" &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - 1),
        
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == TRUE &
          head(CIGAR_1_L, 1) == "M" & tail(CIGAR_1_L, 1) == "S" ~ POS_1,
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "S" &
          tail(CIGAR_1_L, 1) == "M" ~ as.integer(POS_1 + as.integer(tail(CIGAR_1_N, 1)) -1),
        
        SEQ_RCed_1 == TRUE &
          head(CIGAR_1_L, 1) == "S" &
          tail(CIGAR_1_L, 1) == "M" & SATAG_1_1 != "SA" ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "S"  & SATAG_1_1 != "SA" ~ NF_NUMBER,
        
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == TRUE &
          tail(CIGAR_1_L, 1) == "S" &
          head(CIGAR_1_L, 1) == "M" &
          tail(head(CIGAR_1_L, 2), 1) == "M" ~ POS_1,
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "S" &
          tail(CIGAR_1_L, 1) == "M" &
          tail(head(CIGAR_1_L, 2), 1) == "M" ~ as.integer(POS_1 + READ_SPAN_MINUS_I -
                                                            as.integer(head(CIGAR_1_N, 1)) - 1),
        
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "M" ~ as.integer(as.integer(SATAG_1_4) + as.integer(tail(SATAG_1_6_N, 1)) -1),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "M" ~ as.integer(as.integer(SATAG_1_4) + as.integer(tail(SATAG_1_6_N, 1)) -1),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        TRUE ~ ERROR_NUMBER
      )
    ) %>%
    
    #calculate the position of the B flank on its chromosome/plasmid. THis is the position closest to the junction. This is not the final start position, if on the focus contig, this may been adjusting.
    
    mutate(
      FLANK_B_START_POS = case_when(
        CIGAR_1_LEN == 1 & FLANK_A_ORIENT=="FW" ~ as.integer(FlankAUltEnd+1),
        CIGAR_1_LEN == 1 & FLANK_A_ORIENT=="RV" ~ as.integer(FlankAUltEnd-1),
        head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "M" &
          SEQ_RCed_1 == TRUE ~ as.integer(POS_1 + as.integer(head(CIGAR_1_N, 1)) -1),
        head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "M" &
          SEQ_RCed_1 == FALSE ~ as.integer(POS_1 + READ_SPAN_MINUS_I - as.integer(tail(CIGAR_1_N, 1))),
        
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == TRUE &
          head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "S" ~ as.integer(POS_1 + as.integer(head(CIGAR_1_N, 1)) -
                                                   1),
        SEQ_RCed_1 == TRUE &
          head(CIGAR_1_L, 1) == "S" &
          tail(CIGAR_1_L, 1) == "M" & SATAG_1_1 != "SA" ~ NF_NUMBER,
        CIGAR_1_LEN == 2 &
          SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "S" & tail(CIGAR_1_L, 1) == "M" ~ POS_1,
        CIGAR_1_LEN > 2 &
          SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "S" &
          tail(CIGAR_1_L, 1) == "M" &
          tail(head(CIGAR_1_L, 2), 1) == "M" ~ POS_1 + as.integer(tail(head(CIGAR_1_N, 2), 1)),
        SEQ_RCed_1 == FALSE &
          head(CIGAR_1_L, 1) == "M" &
          tail(CIGAR_1_L, 1) == "S"  & SATAG_1_1 != "SA" ~ NF_NUMBER,
        
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "M" ~ as.integer(as.integer(SATAG_1_4) + as.integer(head(SATAG_1_6_N, 1)) -
                                                     1),
        SEQ_RCed_1 == TRUE &
          SATAG_1_1 == "SA" &
          head(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "M" ~ as.integer(SATAG_1_4),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "+" &
          tail(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "M" ~ as.integer(as.integer(SATAG_1_4) + as.integer(head(SATAG_1_6_N, 1)) -
                                                     1),
        SEQ_RCed_1 == FALSE &
          SATAG_1_1 == "SA" &
          tail(CIGAR_1_L, 1) == "S"  &
          SATAG_1_5 == "-" &
          head(SATAG_1_6_L, 1) == "S" & is.na(SATAG_1_12) ~ NF_NUMBER,
        TRUE ~ ERROR_NUMBER
      )
    ) %>%
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
  
  data_improved4 = data_improved3 %>%
    
    #search the SEQ_1 sequence for an intact sequence surrounding the DSB. Then do the same but allow 1 mismatch
    rowwise() %>%
    mutate(DSB_AREA_INTACT = if_else(
      length(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 0)) > 0,
      "TRUE",
      "FALSE"))%>%
    mutate(DSB_AREA_1MM = if_else(
      length(matchPattern(DNAString(DSB_AREA_SEQ), DNAString(SEQ_1), max.mismatch = 1)) > 0 & DSB_AREA_INTACT==FALSE,
      "TRUE",
      "FALSE")) %>%
    mutate(CASE_WT = if_else((
      DSB_AREA_INTACT==TRUE | DSB_AREA_1MM==TRUE),
      TRUE,
      FALSE))%>%
    
    #Determine the chromosome the B flank of the mate is mapped to
    mutate(
      MATE_FLANK_B_CHROM = str_replace_all((case_when(
        CIGAR_2_LEN == 1 ~ RNAME_2,
        (CIGAR_2_LEN > 1 &
           SEQ_RCed_2 == TRUE & head(CIGAR_2_L, 1) == "M") ~ RNAME_2,
        (CIGAR_2_LEN > 1 &
           SEQ_RCed_2 == FALSE & tail(CIGAR_2_L, 1) == "M") ~ RNAME_2,
        (
          CIGAR_2_LEN > 1 &
            SEQ_RCed_2 == TRUE &
            SATAG_2_1 != "SA" & head(CIGAR_2_L, 1) == "S"
        ) ~ "NOT_FOUND",
        (
          CIGAR_2_LEN > 1 &
            SEQ_RCed_2 == FALSE &
            SATAG_2_1 != "SA" & tail(CIGAR_2_L, 1) == "S"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" & SATAG_2_5 == "+"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" & SATAG_2_5 == "-"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == TRUE &
            SATAG_2_1 == "SA" &
            head(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            tail(SATAG_2_12_L, n = 1) == "M" &
            SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" & SATAG_2_5 == "+"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" & SATAG_2_5 == "-"
        ) ~ SATAG_2_3,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" & is.na(SATAG_2_12)
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            tail(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "+" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ "NOT_FOUND",
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "-"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            head(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ SATAG_2_9,
        (
          SEQ_RCed_2 == FALSE &
            SATAG_2_1 == "SA" &
            tail(CIGAR_2_L, 1) == "S" &
            head(SATAG_2_6_L, n = 1) == "M" &
            SATAG_2_5 == "-" &
            tail(SATAG_2_12_L, n = 1) == "M" & SATAG_2_11 == "+"
        ) ~ "NOT_FOUND",
        TRUE ~ "ERROR"
      )), "-", "_")
    ) %>%
    mutate(MATE_FLANK_B_CHROM =if_else(
      MATE_FLANK_B_CHROM != "Mt" & MATE_FLANK_B_CHROM != "Pt" & MATE_FLANK_B_CHROM != PLASMID & MATE_FLANK_B_CHROM != "NOT_FOUND",
      paste0("Chr", MATE_FLANK_B_CHROM),
      MATE_FLANK_B_CHROM)) %>%
    ungroup()
  
  data_improved5 = data_improved4 %>%
    
    #test whether B flanks from read and mate are mapped to same chromosome
    mutate(MATE_FLANK_B_CHROM_AGREE = if_else(FLANK_B_CHROM == MATE_FLANK_B_CHROM,
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
      AvgBaseQual_1_Max,
      MATE_FLANK_B_CHROM_AGREE,
      NrOpticalDuplicates,
      SEQ_2,
      QNAME_1st,
      SEQ_1_LEN,
      FLANK_A_START_POS
    ) %>%
    
    rowwise() %>%
    mutate(FLANK_A_REF_LEN = as.integer(nchar(FLANK_A_REF))) %>%
    
    #Find how much SEQ_1 matches with FLANK_A_REF. Allow 1 bp mismatch somewhere, if the alignment after that continues for at least another 10 bp.
    mutate(FLANK_A_MATCH = matcher_skipper(FLANK_A_REF, SEQ_1)) %>%
    mutate(FLANK_A_LEN = nchar(FLANK_A_MATCH)) %>%
    mutate(FLANK_A_END_POS = if_else(FLANK_A_ORIENT == "FW",
                                     FLANK_A_START_POS + (FLANK_A_LEN -1),
                                     FLANK_A_START_POS - (FLANK_A_LEN -1))) %>%
    mutate(SEQ_1_WO_A = substr(SEQ_1, start = FLANK_A_LEN + 1, stop = SEQ_1_LEN)) %>%
    
    #calculate FLANK A DEL length
    mutate(
      FLANK_A_DEL = case_when(
        CASE_WT == TRUE ~ as.integer(0),
        CASE_WT != TRUE & FLANK_A_LEN != ERROR_NUMBER & FLANK_A_LEN != NF_NUMBER ~ as.integer(PRIMER_TO_DSB - FLANK_A_LEN),
        TRUE ~ ERROR_NUMBER
      )
    ) %>% ungroup()
  
  data_improved6 = data_improved5 %>%
    rowwise() %>%
    #FLANK_B_REF. This ref includes homology.
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
    
    mutate(FLANK_B_MATCH_LEN = nchar(FLANK_B_MATCH)) %>% 
    
    mutate(NO_MATCH = if_else(FLANK_B_MATCH_LEN == 0 & FLANK_B_REF != "NA",
                              TRUE,
                              FALSE)) %>%
    
    filter(NO_MATCH == FALSE) %>%
    
    #flank b start position including MH
    mutate(FLANK_B_START_POS_MH = if (FLANK_B_REF != "NA"){
      if (FLANK_B_ORIENT == "FW"){
        FLANK_B_END_POS-(FLANK_B_MATCH_LEN-1)
      }else{
        FLANK_B_END_POS+(FLANK_B_MATCH_LEN-1)
      }
    }else{
      FLANK_B_START_POS
    }) %>%
    
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
  
  data_improved7 = data_improved6 %>%
    rowwise() %>%
    
    
    
    #calculate FLANK B DEL length.
    mutate(FLANK_B_DEL =
             if (CASE_WT == TRUE) {
               0
             }else if (FLANK_B_REF!="NA"){
               if (FLANK_B_ORIENT == "FW"){
                 as.integer(FLANK_B_START_POS_DEL - (FlankAUltEnd+1))
               }else {
                 as.integer((FlankAUltEnd-1) - FLANK_B_START_POS_DEL)
               }
             }else{
               ERROR_NUMBER
             }) %>%
    
    #calculate total deletion length
    mutate(delSize = if_else(
      FLANK_B_DEL != ERROR_NUMBER,
      as.integer(FLANK_A_DEL + FLANK_B_DEL),
      ERROR_NUMBER
    )) %>%
    #correct MH if delsize =0
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
    
    #filter away cases where the FLANK_B_LEN is too short 
    filter(FLANK_B_LEN_DEL >= FLANK_B_LEN_MIN) %>%
    
    
    #count the number of mismatches that were ignored
    mutate(mismatch_found = if (FLANK_B_MATCH != ""){
      if ( (nrow(as.data.frame(matchPattern(pattern= FLANK_A_MATCH, subject=SEQ_1, max.mismatch=0))))>0 & (nrow(as.data.frame(matchPattern(pattern= FLANK_B_MATCH, subject=SEQ_1, max.mismatch=0))))>0){
        FALSE
      }else{
        TRUE
      }}else{
        FALSE
      })
  
  
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
    rowwise() %>%
    mutate(
      POT_FAKE_INS = if_else(insSize == delSize & delSize > 9,
                             as.character(substr(TOTAL_REF, start = FLANK_A_LEN, stop = FLANK_A_LEN + insSize -1)),
                             "")) %>%
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
    #recalculate the FLANK B del
    mutate(FLANK_B_LEN_DEL = if_else(FLANK_B_ORIENT == "FW",
                                     FLANK_B_END_POS - FLANK_B_START_POS,
                                     FLANK_B_START_POS - FLANK_B_END_POS)) %>%
    
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
        FALSE
      )
    ) %>%
    filter(hasPROBLEM == FALSE) %>%
    
    #sort the rows in order to have the ones with the highest base quality at the top
    arrange(desc(SEQ_1_LEN), desc(AvgBaseQual_1_Max)) %>%
    
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
      DSB_AREA_1MM
    ) %>%
    summarize(
      countEvents = sum(NrOpticalDuplicates),
      Max_BaseQual_1 = max(AvgBaseQual_1_Max),
      QNAME_first = first(QNAME_1st),
      SEQ_1_first = first(SEQ_1),
      SEQ_2_first = first(SEQ_2)
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
      QNAME_first,
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
      Max_BaseQual_1,
      SEQ_1_first,
      SEQ_2_first,
      DSB_AREA_INTACT,
      DSB_AREA_1MM
    ) %>%
    
    #rename column for SIQplotter
    rename(Alias = FILE_NAME)
  
  
  #calucate the fraction of duplicates compared to total
  data_improved9 = data_improved8 %>% group_by(Alias, PRIMER_SEQ) %>% summarize(SumCountEvents =
                                                                                  sum(countEvents))
  
  data_improved10 = left_join(data_improved8, data_improved9, by = c("Alias", "PRIMER_SEQ")) %>%
    mutate(fraction = countEvents / SumCountEvents) %>%
    mutate(SumCountEvents = NULL,
           FlankAUltEnd = FlankAUltEnd,
           FLANK_A_ORIENT = FLANK_A_ORIENT,
           FOCUS_CONTIG = FOCUS_CONTIG,
           Genotype = Genotype,
           program_version = version_no)
  
  
  #write an excel sheet as output
  work_book <- createWorkbook()
  addWorksheet(work_book, "rawData")
  writeData(work_book, sheet = 1, data_improved10)
  saveWorkbook(work_book, file = paste0(output_dir, i, "_CISGUIDE_V", version_no, ".xlsx"), overwrite = TRUE)
  
}

sample_list = list.files(path=output_dir, pattern = "\\.xlsx")
wb = tibble()

for (i in sample_list){
  wb=rbind(wb, read.xlsx(paste0(output_dir, i)))
}
work_book2 <- createWorkbook()
addWorksheet(work_book2, "rawData")
writeData(work_book2, sheet = 1, wb)
saveWorkbook(work_book2, file = paste0(output_dir, "CISGUIDE_V", version_no, ".xlsx"), overwrite = TRUE)
