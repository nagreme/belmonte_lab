# --------------------------------------------------------------
# Date: 2017-08-23
# Nadège Pulgar-Vidal
# 
# Visualize correlation between mean methylation on a gene/flank 
# and gene expression
#
# Filtering cuffdiff's genes.fpkm_tracking (prep step)
# tracking_id    class_code      nearest_ref_id  gene_id#        gene_short_name tss_id  locus   length  coverage        GLOB_FPKM*      GLOB_conf_lo    GLOB_conf_hi    GLOB_status     MG_FPKM*        MG_conf_lo      MG_conf_hi      MG_status
# Bna00158bs0010 -       -       Bna00158bs0010  Bna00158bs0010  -       Scaffold00158b:34951-36345      -       -       1.47416 0.490069        2.45825 OK      0.113054        0       0.43742 OK
# 
# for FILE in genes.fpkm_tracking; do awk '{print $4"\t"$10"\t"$14}' "$FILE" > "$FILE.less"; done
# 
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

# Mean Methylation Ratio Input Files

glob_gene_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_genes_CpG_mean_meth_ratio.txt.less'
glob_gene_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_genes_CHG_mean_meth_ratio.txt.less'
glob_gene_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_genes_CHH_mean_meth_ratio.txt.less'
glob_flank_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_flanks_CpG_mean_meth_ratio.txt.less'
glob_flank_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_flanks_CHG_mean_meth_ratio.txt.less'
glob_flank_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/GLOB_flanks_CHH_mean_meth_ratio.txt.less'

mg_gene_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_genes_CpG_mean_meth_ratio.txt.less'
mg_gene_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_genes_CHG_mean_meth_ratio.txt.less'
mg_gene_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_genes_CHH_mean_meth_ratio.txt.less'
mg_flank_CpG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_flanks_CpG_mean_meth_ratio.txt.less'
mg_flank_CHG_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_flanks_CHG_mean_meth_ratio.txt.less'
mg_flank_CHH_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/mean_meth_ratios/MG_flanks_CHH_mean_meth_ratio.txt.less'

# --- Mean Methylation Ratio Files Expected Format ---
# BnaN01g00010	0.02929113924
# BnaN01g00020	0.01047058824
# BnaN01g00030	0.006074626866

input_files <- c(glob_gene_CpG_file, glob_gene_CHG_file, glob_gene_CHH_file,
                 glob_flank_CpG_file, glob_flank_CHG_file, glob_flank_CHH_file,
                 mg_gene_CpG_file, mg_gene_CHG_file, mg_gene_CHH_file,
                 mg_flank_CpG_file, mg_flank_CHG_file, mg_flank_CHH_file)


# --- FPKM File
fpkm_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/2017_08_18_cuffdiff/genes.fpkm_tracking.less'

# --- FPKM File Expected Format ---
# gene_id	GLOB_FPKM	MG_FPKM
# Bna00158bs0010	1.47416	0.113054
# Bna00158bs0020	0	0
# Bna00158bs0030	0	0


# --- Output parameters

tags <- c("GLOB_genes_CpG", "GLOB_genes_CHG", "GLOB_genes_CHH",
          "GLOB_flanks_CpG", "GLOB_flanks_CHG", "GLOB_flanks_CHH",
          "MG_genes_CpG", "MG_genes_CHG", "MG_genes_CHH",
          "MG_flanks_CpG", "MG_flanks_CHG", "MG_flanks_CHH")

samples <- c(rep("GLOB", 6), rep("MG", 6))


outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/visualizations/scatterplots/"


# ==========================
# --- PROCESS
# ==========================

# Ready the fpkm file
read_delim(fpkm_file, delim = "\t", col_names = T) %>%
  gather(key = sample, value = FPKM, -gene_id) ->
  fpkm_data

# --- Format after gather 
# gene_id sample        FPKM
# <chr>  <chr>       <dbl>
#   1 Bna00158bs0010   GLOB    1.474160
# 2 Bna00158bs0020   GLOB    0.000000
# 3 Bna00158bs0030   GLOB    0.000000

index <- 1

# Read a file, join fpkms onto, and generate a scatterplot
for (file in input_files)
{
  read_delim(file, delim = "\t", col_names = c("gene_id", "mean_meth_ratio")) %>%
    filter(mean_meth_ratio != ".") %>% #just exclude the lines the features that don't have any data (.)
    mutate(mean_meth_ratio = as.double(mean_meth_ratio)) ->
    mean_meth_ratio_data
  
  fpkm_data %>%
    filter(sample == samples[[index]]) %>%
    inner_join(mean_meth_ratio_data, .) ->
    plot_data
  
  plot_data$FPKM <- log(plot_data$FPKM +1, 10)
    
    graph <- ggplot(plot_data, aes(mean_meth_ratio, FPKM)) +
      geom_point() +
      theme(panel.background = element_rect(fill = "white")) +
      labs(title = "Mean Methylation Ratio and Gene Expression Correlation",
           subtitle = tags[[index]],
           x = "Mean methylation ratio",
           y = "FPKM") 
    
    ggsave(paste0(outfile_path, tags[[index]], "_log_scatter.png")) #files are huge as pdfs so do pngs for now

    index <- index + 1
}







