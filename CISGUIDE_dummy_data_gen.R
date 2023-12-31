###############################################################################
#install and load packages
###############################################################################
if (require(BSgenome)==FALSE){install.packages("BSgenome")}
if (require(Biostrings)==FALSE){install.packages("Biostrings")}
if (require(stringi)==FALSE){install.packages("stringi")}
if (require(stringdist)==FALSE){install.packages("stringdist")}
if (require(tidyverse)==FALSE){install.packages("tidyverse")}
if (require(openxlsx)==FALSE){install.packages("openxlsx")}

#producing dummy data

WT_SEQ=  "GATCTTTTATCAAAACATTTATAGCAGGGTAGACGAAGATATGACTCACCTAGAGAAGTAATGTATGTCAGTCCCATTTTTGTGTATTTGTATGCGATCCCCATCAAATTTGCATTCAGCAACTACATCTTTCCCATGCAGCTGTGTCGAGGAAGAGACTCATGCATTGGTTAATGAAATAAAACACCAATACAACAACATTAACAAACAAACAAAATACTGCACATGGAAATGGAGTAAGCTGACAATAGAGTTGATGGAATTTACCTTCTTCCAAGCAGCATTTACGTCACCAATACGCATAGCTAACTGGGGTCGTACGGCTTTTCCAACCTCTATATCCTACAAGAGAAGCACAGTTTTTATGCTTATCAGTTTAAGATGATTTCAATAAATCGATTGTATCTCGAGGAAATGAAACAGTAGTTTGTAAGTTGTCATTCAAGAATTTCTAAATACAAAACCTGGCGCTTGTGCCGTTGGTGTCGATCCCTCAACTTTTCACAGACTAGTTTCAAATCACATGTGACATTAAACAAGTCCTCAGCATCTGGATGGAACTCCTGAAAAATGCTCTTCTCACTCATTCCCAGTTTCAAATCTGTAAACCATTTCCCAACAAAATCAGATTGATCTCAATAACTTTAAAGAAAAAAGATGGAAACTTTATGTGAAGAAACTGCAAAGTTCCTACCTTTTAGAATAATCCTAATGACCCACTTCATCTCCTGAGCATTTGTCTTCTGAATCAATGTAGAAAGAACCAAAGTTTTCTCCGCTCTAGATCCCAAAAGAAAACAAGATAAGAGAGTAAACAGAAACAAAATTGATTCTTTATGATTAGAGAAGCACCTGTTCTCACTTGAAGCCAAACGATCAAGCAAATCATTCAATTCCTTAATAGTCAAACCGCCAGAAGCCATTCCTTGTCTACGTTGCAATACCTATGAGCACAAAACTGTTCAAAAACCCTAAATCCCCCAAAAACCAAAACCCAGAAAGAAAGTCAAACAAATCCTAAAACCTCAGCAGCAATTAAGGAGAAGTTTCCAGCATTGGCTCCAGCTTTAGCAGTTCCTCCTTTACGCCAATTGAGAAGACGAACAGCATCAGGAGCATCACGCGAGATACCAAGTGCGTCGATCAGACACGTGGCGAGCACTGACTCTTTAAGGCCGTAGCTACCTCGTTCTCGATCAAGCGACGGAATGATTAAACGGACGGCGACGAAGTAGTCGGAGGGTTTGCAGTAAGTGTCGAGAAACTTGCGGAATTTGGATCGTTTCTGAGATGAGGTTTTGCTTTTTTGAATCCAGTTAAAGAGAGATACCAGTACGCTGAATTTGATCTCCTCCGTCATCTTACTTTGCCGGAATTATAGATTCTTCCGGAGAAGATGGTTAATCAGGGAAATAGTGACGGAGAATCACTGAGCAGCGGCGTGAATTGCTTGCCGGAGAAGGAAGTCGCGTTTGAAAAGTTGCACGGCCGGTTTTG"
LIG4_SEQ="GATCTTTTATCAAAACATTTATAGCAGGGTAGACGAAGATATGACTCACCTAGAGAAGTAATGTATGTCAGTCCCATTTTTGTGTATTTGTGGTGTAAACAAATTGACGCTTAGACAACTTAATAACACATTGCGGACGTTTTTAATGTACTGGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGCCCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCAAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCCCGATCTAGTAACATAGATGACACCGCGCGCGATAATTTATCCTAGTTTGCGCGCTATATTTTGTTTTCTATCGCGTATTAAATGTATAATTGCGGGACTCTAATCATAAAAACCCATCTCATAAATAACGTCATGCATTACATGTTAATTATTACATGCTTAACGTAATTCAACAGAAATTATATGATAATCATCGCAAGACCGGCAACAGGATTCAATCTTAAGAAACTTTATTGCCAAATGTTTGAACGATCGGGGAAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCCCCCGTGTTCTCTCCAAATGAAATGAACTTCCTTATATAGAGGAAGGGTCTTGCGAAGGATAGTGGGATTGTGCGTCATCCCTTACGTCAGTGGAGATATCACATCAATCCACTTGCTTTGAAGACGTGGTTGGAACGTCTTCTTTTTCCACGATGCTCCTCGTGGGTGGGGGTCCATCTTTGGGACCACTGTCGGCAGAGGCATCTTCAACGATGGCCTTTCCTTTATCGCAATGATGGCATTTGTAGGAGCCACCTTCCTTTTCCACTATCTTCACAATAAAGTGACAGATAGCTGGGCAATGGAATCCGAGGAGGTTTCCGGATATTACCCTTTGTTGAAAAGTCTCAATTGCCCTTTGGTCTTCTGAGACTGTATCTTTGATATTTTTGGAGTAGACAAGTGTGTCGTGCTCCACCATGTTGACGAAGATTTTCTTCTTGTCATTGAGTCGTAAGAGA"

#R2 reads start with the normal seq, R1 is on the other side and RCed


#first a column with fragment lengths and another column with genotype
#then based on these columns, generate the fragment seq
#then based on the fragment seq, generate the read seqs
#generate quality scores (just all Gs for the same length)
#generate read names ("@M02948:174:000000000-JBDYN:1:1101:15770:"plus rownumber)
#then order the data in a way that is suitabf

n_list=c(0.75, 0.90, 0.95, 0.99)

for (i in n_list) {
n_high=i
n_high_unders=str_replace(formatC(n_high, format="f", digits=2), "\\.", "_")

output.file1 <- file(paste0("dummy/dummy_data_", n_high_unders, "_R1.fastq"), "wb")  
output.file2 <- file(paste0("dummy/dummy_data_", n_high_unders, "_R2.fastq"), "wb")  


draw_box = c("wt", "lig4")
frag_len = rep(150:600, each=50)
total_reads=length(frag_len)
read_num=1:total_reads
geno=sample(draw_box, size=total_reads, replace=TRUE, prob=c(n_high, 1-n_high))
QUAL=rep(paste0(rep("G", each=150), collapse=""), total_reads)



data1=tibble(read_num, frag_len, geno, QUAL) %>% 
  rowwise()%>%
  mutate(frag_seq= case_when(geno=="wt" ~ substr(WT_SEQ, 1, frag_len),
                             TRUE ~ substr(LIG4_SEQ, 1, frag_len)))%>%
  mutate(R2_seq = substr(frag_seq, 1, 150))%>%
  mutate(R1_seq = as.character(reverseComplement(DNAString(substr(frag_seq, nchar(frag_seq)- 149, nchar(frag_seq)))))) %>%
  mutate(read_name = paste0("@M02948:174:000000000-JBDYN:1:1101:15770:", formatC(read_num, width=4, format="d", flag="0"))) %>%
  mutate(plus = "+")

data1_R1=data1 %>%
  select(read_name, R1_seq, plus, QUAL) %>%
  pivot_longer(everything()) %>%
  select(!name)

data1_R2=data1 %>%
  select(read_name, R2_seq, plus, QUAL) %>%
  pivot_longer(everything()) %>%
  select(!name)



write.table(data1_R1, file=output.file1, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(data1_R2, file=output.file2, quote=FALSE, row.names=FALSE, col.names=FALSE)

close(output.file1)
close(output.file2)
}

