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

# Bin size (absolute positions)
bin_size <- 100000

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

# Name and location of output plot file
outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/visualizations/output_plot.pdf"

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

all_genes_fpkms <- read_delim(gene_fpkm_file, delim = "\t")

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
# --- CREATE BED5 FILE
# ==========================

gene_annotation_full %>%
  select(chr, start, end, gene_id) %>% 
  left_join(diff_exp_genes_fpkms, .) %>%
  mutate(score = paste0(GLOB_fpkm,";",MG_fpkm)) %>%
  select(chr, start, end, gene_id, score) %>%
  write_delim(bed5_outfile, delim = "\t")

# Note: Don't forget to remove the top row (column names) before using with bedtools


# ==========================
# --- PROCESS 
# ==========================

# Note: For some reason there are slightly more genes in our annotation than in our expression data
# To get these: anti_join(gene_annotation_full, all_genes_fpkms)
# There are 346 of them

left_join(all_genes_fpkms, gene_annotation_full) %>%
  select(gene_id, chr, start, end, GLOB, MG) %>%
  gather(sample, fpkm, -gene_id, -chr, -start, -end) ->
  genes_full

# Put the data into bins and calculate avg expression in the bin
genes_full %>%
  filter(fpkm <= 1000)
  mutate(bin = as.integer(start/bin_size)) %>%
  group_by(sample, chr, bin) %>%
  summarise(avg_fpkm = log(mean(fpkm)+1, 10)) ->
  genes_full_binned

# For reordering the chr names
genes_full_binned$chr_order <- factor(genes_full_binned$chr, 
                                            levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                                       "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

genes_full_binned <- na.omit(genes_full_binned)

# Scale bin counts so scale is consistent when splitting A and C genome
# Use logs for fpkms
# genes_full_binned %>%

# genes_full_binned$avg_fpkm <- unlist(lapply(genes_full_binned$avg_fpkm, function(x) x/max(genes_full_binned$avg_fpkm)))


# Separate data into 2 sets by chromosome (N1-N9, N10-N19)
genes_full_n1_n10 <- filter(genes_full_binned, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"))
genes_full_n11_n19 <- filter(genes_full_binned, chr %in% c("N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))



# ==========================
# --- VISUALIZE 
# ==========================

plot_data <- genes_full_binned
plot_data <- genes_full_n1_n10
plot_data <- genes_full_n11_n19


# Bar Chart
graph <- ggplot(plot_data, aes(bin, avg_fpkm)) +
  facet_grid(sample ~ chr_order, scales = "free_x") +
  geom_bar(stat = "identity", aes(fill = avg_fpkm, colour = avg_fpkm)) +
  guides(fill=FALSE) + 
  labs(title = "Gene Expression",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Average bin expression") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_color_gradient2(low="#021637", 
                        midpoint = 1,
                        mid="#00517C",
                        high="#C5F6D5") #white low, dark purple high


ggsave(outfile_path, height = 3, width = 20) 

