# --------------------------------------------------------------
# Date: 2017-08-02
# Nadège Pulgar-Vidal
# 
# Methylation Heatmap and Bar Graph
# --------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)


# List of files
files <- c("/home/mark/Documents/Nadege/methylation_project/data/GLOB_meth.txt",
          "/home/mark/Documents/Nadege/methylation_project/data/MG_meth.txt")

# List of samples (matching files)
sample_tags <- c("GLOB",
                "MG")

data_list <- list()

# For each file
for (file in files)
{
  # read in as tibble/data.frame
  # clusts <- read.table(file=clusts_file, row.names = 1, header=T, sep="\t")
  dat <- read.table(file, header=T, sep="\t")
  
}









