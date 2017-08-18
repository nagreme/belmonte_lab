# --------------------------------------------------------------
# Date: 2017-08-16
# Nadège Pulgar-Vidal
# 
# Visualize correlation between methylation and gene expression
# --------------------------------------------------------------

# ==========================
# --- LIBRARIES
# ==========================

library(magrittr)
library(tidyverse)
library(stringr)

# ==========================
# --- SETUP
# ==========================

# Input files

gene_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CpG_hypo.txt.orig'
gene_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CpG_hyper.txt.orig'
gene_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CHG_hypo.txt.orig'
gene_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CHG_hyper.txt.orig'
gene_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CHH_hypo.txt.orig'
gene_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_genes_CHH_hyper.txt.orig'
flank_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CpG_hypo.txt.orig'
flank_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CpG_hyper.txt.orig'
flank_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CHG_hypo.txt.orig'
flank_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CHG_hyper.txt.orig'
flank_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CHH_hypo.txt.orig'
flank_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/orig_lists/GLOB_vs_MG_flanks_CHH_hyper.txt.orig'


files <- list(gene_CpG_hypo_file, gene_CpG_hyper_file,
                  gene_CHG_hypo_file, gene_CHG_hyper_file,
                  gene_CHH_hypo_file, gene_CHH_hyper_file,
                  flank_CpG_hypo_file, flank_CpG_hyper_file,
                  flank_CHG_hypo_file, flank_CHG_hyper_file,
                  flank_CHH_hypo_file, flank_CHH_hyper_file)

# Output File

outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/visualizations/scatterplots/"


# Processing Parameters

tags <- c("CpG;gene;hypo","CpG;gene;hyper",
          "CHG;gene;hypo","CHG;gene;hyper",
          "CHH;gene;hypo","CHH;gene;hyper",
          "CpG;flank;hypo","CpG;flank;hyper",
          "CHG;flank;hypo","CHG;flank;hyper",
          "CHH;flank;hypo","CHH;flank;hyper")

context_index <- 1
feat_type_index <- 2
meth_type_index <- 3 


# Plotting Parameters

samples <- c("GLOB", "MG")
contexts <- c("CpG", "CHG", "CHH")
feat_types <- c("gene", "flank")
# meth_types <- c("hypo", "hyper")
chrs <- c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
          "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19")


# ==========================
# --- READ
# ==========================

# --- Expected input format ---
# chr.1	name	meth.diff	fpkms
# N1	BnaN01g09510	-33.64905284147557	1;3
# N1	BnaN01g17220	-60	4;47
# N1	BnaN01g37390	-42.13653583465869	7;67
# N1	BnaN01g38070	-31.47363465160075	20;7
# N1	BnaN01g38350	-43.68139653585639	49;77
# N1	BnaN01g38350	-38.43762145355616	49;77
# N10	BnaN10g12760	-32.64486517498567	10;29
# N10	BnaN10g20050	-39.32350718065004	94;239

data_list <- list()

# TODO: Find a better way to do this? Or more R-like?
tag_index <- 1 

# For each file
for (file in files)
{
  # read in as tibble/data.frame
  dat <- read_delim(file, delim = "\t")
  
  curr_tags <- str_split(tags[tag_index], ";")[[1]]
  
  # extract the columns we want and add a sample tag
  dat %<>%
    mutate(context = curr_tags[context_index]) %>%
    mutate(feat_type = curr_tags[feat_type_index]) %>%
    mutate(meth_type = curr_tags[meth_type_index]) 
  
  # add it to our growing list
  data_list[[tag_index]] <- dat
  
  # move to the next tag (and next position in data_list)
  tag_index <- tag_index + 1
}

# Collect data into one big table
data_full <- bind_rows(data_list)


# ==========================
# --- PREP
# ==========================

data_full$GLOB <- unlist(lapply(data_full$fpkms, function(x) as.integer(str_split(x, ";")[[1]][1])))
data_full$MG <- unlist(lapply(data_full$fpkms, function(x) as.integer(str_split(x, ";")[[1]][2])))

data_full %<>%
  select(-name, -fpkms) %>%
  gather(sample, fpkm, -chr.1, -context, -feat_type, -meth_type, -meth.diff) 


# --- Expected format after processing --- 
# chr.1 meth.diff context feat_type meth_type sample  fpkm
# <chr>     <dbl>   <chr>     <chr>     <chr>  <chr> <int>
#   1    N1 -30.51615     CpG      gene      hypo   GLOB     4
# 2    N1 -38.77150     CpG      gene      hypo   GLOB     1
# 3    N1 -34.11888     CpG      gene      hypo   GLOB     6
# 4    N1 -26.00619     CpG      gene      hypo   GLOB     4
# 5    N1 -32.05882     CpG      gene      hypo   GLOB     7
# 6    N1 -40.35088     CpG      gene      hypo   GLOB     5
# 7    N1 -34.44444     CpG      gene      hypo   GLOB     8

# To ensure chromosomes are in the right order
data_full$chr_order <- factor(data_full$chr.1, 
                              levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                         "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))



# ==========================
# --- VISUALIZATION
# ==========================


plot_scatter <- function(all_data, chr, contxt, ft_type, sampl)
{
  plot_data <- filter(all_data, chr.1 == chr, context == contxt,
                      feat_type == ft_type, sample == sampl)
  
  graph <- ggplot(plot_data, aes(fpkm, meth.diff)) +
    geom_point(aes(col=meth_type)) +
    # facet_grid(context + sample ~ chr_order) +
    theme(panel.background = element_rect(fill = "white")) +
    labs(title = "Methylation Ratio and Gene Expression Correlation",
         subtitle = paste0(contxt, " methylation on ", 
                           ft_type,"s on chromosome ", chr, " in ", sampl),
         x = paste0("fpkm"),
         y = "meth.diff") 
  
  ggsave(paste0(outfile_path, contxt, "_", sampl, "_", chr,"_", ft_type,"_scatter.pdf"), width= 7, height = 7) 
}


# Plot all posibilities

for (chrom in chrs)
{
  for (con in contexts)
  {
    for (ftt in feat_types)
    {
      for (s in samples)
      {
        plot_scatter(data_full, chrom, con, ftt, s)
      }
    }
  }
}




