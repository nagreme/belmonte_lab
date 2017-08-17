# --------------------------------------------------------------
# Date: 2017-08-16
# Nadège Pulgar-Vidal
# 
# Build a BED5 file for gene expression (differential expression)
# from a bunch of files containing different parts of the necessary 
# info (Frankenstein).
# And visualize gene expression levels in bins across chromosomes
# --------------------------------------------------------------

# ==========================
# --- LIBRARIES
# ==========================

library(readr)
library(dplyr)
library(magrittr)
library(stringr)


# ==========================
# --- SETUP
# ==========================

# Input

# contains a 0/1 matrix for whether a gene is differentially expressed (only look at relevant column)
diff_exp_gene_matrix_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/gene_expression/differential_exp_sig.txt'

# contains GLOB and MG fpkm (integers) by gene
gene_fpkm_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/gene_expression/GLOB_MG_FPKM.txt'

# contains chr, start, end, by gene name (includes scaffolds)
gene_annotation_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/annotation/DH12075_annotation_genes.gff3'


# Output

# BED5 file output
bed5_outfile <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/gene_expression/GLOB_vs_MG_genes_with_scaffolds.bed'

# ==========================
# --- READ/PREP
# ==========================

# Read the relevant column in the gene matrixand pull out 
# a list of differentially expressed genes
read_delim(diff_exp_gene_matrix_file, delim = "\t") %>% 
  select(gene_id, GLOBvsMG) %>% 
  filter(GLOBvsMG == 1) %>% 
  select(gene_id) %>%
  unlist() ->
  diff_exp_genes

# Read in the fpkms and filter for differentially expressed genes (^^^)
read_delim(gene_fpkm_file, delim = "\t") %>% 
  filter(gene_id %in% diff_exp_genes) ->
  diff_exp_genes_fpkms

# Read in the relevant columns of the annotation file
read_delim(gene_annotation_file, delim = "\t",
           col_names = c("chr", "src", "feat_type", "start" ,"end", 
                         "score", "strand", "myst", "id_name")) %>% 
  select(chr, start, end, id_name) ->
  # mutate(gene_id = str_match(id_name, "ID=(Bna.*);Name=.*")[[,2]]) %>% # separate the gene id for later joins
  gene_annotation_full

# Haven't figured out how (if possible) to do this in a mutate call, so...
gene_annotation_full$gene_id <- unlist(lapply(gene_annotation_full$id_name, function(x) str_match(x, "ID=(Bna.*);Name=.*")[[2]]))


# ==========================
# --- CERATE BED5 files
# ==========================

gene_annotation_full %>%
  select(chr, start, end, gene_id) %>% 
  inner_join(diff_exp_genes_fpkms, .) %>%
  mutate(score = paste0(GLOB_fpkm,";",MG_fpkm)) %>%
  select(chr, start, end, gene_id, score) %>%
  write_delim(bed5_outfile, delim = "\t")

# Note: Don't forget to remove the top row (column names) before using with bedtools





