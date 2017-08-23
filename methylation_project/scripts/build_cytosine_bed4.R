# --------------------------------------------------------------
# Date: 2017-08-22
# Nadège Pulgar-Vidal
# 
# Create bed4 files from the bsmap methylation call files for use
# with bedtools map (to find mean meth.ratio across feature)
# --------------------------------------------------------------

# ==========================
# --- LIBRARIES
# ==========================

library(tidyverse)
library(magrittr)


# ==========================
# --- SETUP
# ==========================

glob_meth_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/GLOB_meth.txt'
mg_meth_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/MG_meth.txt'

glob_Cs_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/GLOB_meth_Cs.bed'
mg_Cs_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/MG_meth_Cs.bed'

  

# ==========================
# --- PROCESS
# ==========================

infile <- glob_meth_file
outfile <- glob_Cs_file

read_delim(infile, delim ="\t", col_names = T) %>%
  select(chr, pos, context, ratio, eff_CT_count) %>%
  filter(eff_CT_count >= 10) %>%
  mutate(start = pos, end = pos) %>%
  select(chr, start, end, ratio, eff_CT_count ->
  dat #for checking


dat %>%
  filter(sample == "GLOB") %>%
  select(chr, start, end, ratio) %>%
  write_delim(outfile, delim="\t", col_names = F)




















