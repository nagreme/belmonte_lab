# --------------------------------------------------------------
# Date: 2017-08-22
# Nadège Pulgar-Vidal
#
# Create bed4 files from the bsmap methylation call files for use
# with bedtools map (to find mean meth.ratio across feature)
#
# --------------------------------------------------------------
#
# EDIT 2019-07-25
#
# This is a slightly modified version that will be used to build
# bed4 files of Deirdre's bsmap methylation data so that we can
# look at the methylation levels in Dylan's data (siRNAs).
#
# The change is essentially the threshold used for filtering: the
# original script uses an eff_CT_count >= 10 but this one will used
# a threshold of 4 (?)
#
# --------------------------------------------------------------

# ==========================
# --- LIBRARIES
# ==========================

library(tidyverse)
library(magrittr)


# ==========================
# --- SETUP
# ==========================

glob_meth_file <- '/home/belmontes/Desktop/nadege/methylation_data/methylation_calls_bsmap/GLOB_meth.txt'
mg_meth_file <- '/home/belmontes/Desktop/nadege/methylation_data/methylation_calls_bsmap/MG_meth.txt'

glob_Cs_file <- '/home/belmontes/Desktop/nadege/belmonte_lab/siRNA_project/methylation_and_siRNAs/GLOB_meth_Cs_'
mg_Cs_file <- '/home/belmontes/Desktop/nadege/belmonte_lab/siRNA_project/methylation_and_siRNAs/MG_meth_Cs_'

FILTER_THRESHOLD <- 4



# ==========================
# --- PROCESS
# ==========================

# Selection
# infile <- glob_meth_file
# outfile <- glob_Cs_file

infile <- mg_meth_file
outfile <- mg_Cs_file

# Read in the file, filetring by context and create the needed start and end bed columns
read_delim(infile, delim ="\t", col_names = T) %>%
  select(chr, pos, context, ratio, eff_CT_count) %>%
  filter(eff_CT_count >= FILTER_THRESHOLD) %>%
  mutate(start = pos, end = pos) %>%
  select(chr, start, end, ratio, context) ->
  dat

# Write out each context to separate files
dat %>%
  filter(context == "CG") %>%
  select(chr, start, end, ratio) %>%
  write_delim(paste0(outfile, "CpG.bed"), delim="\t", col_names = F)

dat %>%
  filter(context == "CHG") %>%
  select(chr, start, end, ratio) %>%
  write_delim(paste0(outfile, "CHG.bed"), delim="\t", col_names = F)

dat %>%
  filter(context == "CHH") %>%
  select(chr, start, end, ratio) %>%
  write_delim(paste0(outfile, "CHH.bed"), delim="\t", col_names = F)
