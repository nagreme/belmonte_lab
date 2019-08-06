# ------------------------------------------------------------------
# Date: 2019-08-02
# Nadège Pulgar-Vidal
#
# R script to compute the average methylation ratio per siRNA locus
# from Shortstack output based on the output of intersectBed using
# the bed files of both data sets I prepared.
#
# intersectBed -wao -a 2of3_MG.bed -b MG_meth_CpG.bed MG_meth_CHG.bed MG_meth_CHH.bed -names CpG CHG CHH > MG_intersect_results.txt
#
# intersectBed -wao -a 2of3_GLOB.bed -b GLOB_meth_CpG.bed GLOB_meth_CHG.bed GLOB_meth_CHH.bed -names CpG CHG CHH > GLOB_intersect_results.txt
#
# ------------------------------------------------------------------

library(tidyverse)
library(magrittr)

infile <- "~/Desktop/nadege/belmonte_lab/siRNA_project/intersection_output/GLOB_intersect_results.tsv"
outfile <- "~/Desktop/nadege/belmonte_lab/siRNA_project/intersection_output/GLOB_intersect_avg_meth_ratios.tsv"
# infile <- "~/Desktop/nadege/belmonte_lab/siRNA_project/intersection_output/MG_intersect_results.tsv"
# outfile <- "~/Desktop/nadege/belmonte_lab/siRNA_project/intersection_output/MG_intersect_avg_meth_ratios.tsv"

A_subgenome_chrs <- c('N1','N2','N3','N4','N5','N6','N7','N8','N9','N10')
C_subgenome_chrs <- c('N11','N12','N13','N14','N15','N16','N17','N18','N19')

# read file and do initial clean up
read_tsv(infile, col_names=c('A_chr','A_start','A_end','A_locus','context','B_chr','B_start','B_end','B_meth_ratio','ignore')) %>%
  select(-ignore) %>% # removes ignore column
  mutate(B_meth_ratio=replace(B_meth_ratio, B_meth_ratio == '.', -1)) -> # '.' means intersect so say -1 meth_ratio to reflect that
dat

# make it behave
dat$B_meth_ratio <- as.numeric(dat$B_meth_ratio)

# can't group by only locus name because then all the nohits loci get lumped together
dat %>%
  group_by(A_chr,A_start,A_end,A_locus,context) %>% # each unique combination of these will be a row
  summarize(avg_meth_ratio=mean(B_meth_ratio)) %>% # for each row, the average of the group meth ratios is calculated
  mutate(subgenome = case_when(A_chr %in% A_subgenome_chrs ~ "A", A_chr %in% C_subgenome_chrs ~ "C")) -> # label subgenomes
intersection_meth_ratios

intersection_meth_ratios %>%
  # I can't omit A_start and A_end apparently because they're grouping values
  select(subgenome, A_chr, A_start, A_end, A_locus, context, avg_meth_ratio) %>%
  write_tsv(outfile)
