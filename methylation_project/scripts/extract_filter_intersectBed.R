# --------------------------------------------------------------
# Date: 2017-08-15
# Nadège Pulgar-Vidal
# 
# Read in the results of our bedtools intersect and extract gene 
# names/id/accession separating hyper/hypo methylated genes.
# 
# Note: In our data GLOB was treated as control and MG was treatment
#
# Ran sed -i.orig -r "s/.*ID=(Bna.*);.*$/\1/" on all output 
# files to get simple gene list
# sed -i.orig -r "s/.*(Bna.*)\t.*\t.*$/\1/" * (for the second set)
#
# This was used to remove the first line of each file: (run in output dir)
# for FILE in *; do tail -n +2 "$FILE" > "$FILE.tmp" && mv "$FILE.tmp" "$FILE"; done
# --------------------------------------------------------------

# ==========================
# --- LIBRARIES
# ==========================

library(tidyverse)
library(magrittr)
library(stringr)


# ==========================
# --- SETUP
# ==========================

# --- Input files 
gene_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/genes_CpG_intersect.txt'
gene_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/genes_CHG_intersect.txt'
gene_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/genes_CHH_intersect.txt'
flank_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/flanks_CpG_intersect.txt'
flank_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/flanks_CHG_intersect.txt'
flank_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/flanks_CHH_intersect.txt'

# gene_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_genes_CpG_intersect.txt'
# gene_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_genes_CHG_intersect.txt'
# gene_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_genes_CHH_intersect.txt'
# flank_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_flanks_CpG_intersect.txt'
# flank_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_flanks_CHG_intersect.txt'
# flank_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/GLOB_vs_MG_flanks_CHH_intersect.txt'


in_files <- list(gene_CpG_file, gene_CHG_file, gene_CHH_file, 
                 flank_CpG_file, flank_CHG_file, flank_CHH_file)


# --- Output files
gene_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CpG_hypo.txt'
gene_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CpG_hyper.txt'
gene_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CHG_hypo.txt'
gene_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CHG_hyper.txt'
gene_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CHH_hypo.txt'
gene_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/genes_CHH_hyper.txt'
flank_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CpG_hypo.txt'
flank_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CpG_hyper.txt'
flank_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CHG_hypo.txt'
flank_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CHG_hyper.txt'
flank_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CHH_hypo.txt'
flank_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/annotation/methylated_lists/flanks_CHH_hyper.txt'

# gene_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CpG_hypo.txt'
# gene_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CpG_hyper.txt'
# gene_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CHG_hypo.txt'
# gene_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CHG_hyper.txt'
# gene_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CHH_hypo.txt'
# gene_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_genes_CHH_hyper.txt'
# flank_CpG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CpG_hypo.txt'
# flank_CpG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CpG_hyper.txt'
# flank_CHG_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CHG_hypo.txt'
# flank_CHG_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CHG_hyper.txt'
# flank_CHH_hypo_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CHH_hypo.txt'
# flank_CHH_hyper_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/bedtools_intersects/differential_expression/methylated_lists/GLOB_vs_MG_flanks_CHH_hyper.txt'


out_files <- list(gene_CpG_hypo_file, gene_CpG_hyper_file,
                  gene_CHG_hypo_file, gene_CHG_hyper_file,
                  gene_CHH_hypo_file, gene_CHH_hyper_file,
                  flank_CpG_hypo_file, flank_CpG_hyper_file,
                  flank_CHG_hypo_file, flank_CHG_hyper_file,
                  flank_CHH_hypo_file, flank_CHH_hyper_file)


# ==========================
# --- PROCESS 
# ==========================

index <- 1

for (file in in_files)
{
   # Read one file
  read.table(file, header=F, sep="\t", 
             col.names = c("chr.1", "DMR_start", "DMR_end", "meth.diff",
                           "chr.2", "src", "feature_type", "feat_start",
                           "feat_end", "score", "strand", "myst", "name",
                           "num_overlap_bp")) %>%
    select(chr.1, name, meth.diff, num_overlap_bp) ->
             # col.names = c("chr.1", "DMR_start", "DMR_end", "meth.diff",
             #               "chr.2", "start", "end", "name", "fpkms",
             #               "num_overlap_bp")) %>%
    # select(chr.1, name, meth.diff, fpkms) -> # pull out the pieces we need
    dat
  
  #Run this if the name column is already simple other than "_flank" tags to remove those
  dat$name <- unlist(lapply(dat$name, function(x) str_replace(x, "_flank", "")))
  
  
  # Filter hypomethylated genes and write those out
  dat %>%
    filter(meth.diff < 0) %>%
    write_delim(out_files[[index]], delim="\t")
  
  index <- index + 1
  
  # Filter hypermethylated genes and write those out
  dat %>%
    filter(meth.diff > 0) %>%
    write_delim(out_files[[index]], delim="\t")
  
  # Note: The index+1 is because the hypo methylated file appears before the hyper one
  # This assumes that files are in order
  
  index <- index + 1
}



