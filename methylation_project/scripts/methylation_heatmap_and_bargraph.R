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
files <- c("/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/GLOB_test.txt",
          "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/MG_test.txt")

# List of samples (matching files)
sample_tags <- c("GLOB",
                "MG")

data_list <- list()

#Find a better way to do this plz...
tag_index <- 1 

# For each file
for (file in files)
{
  # read in as tibble/data.frame
  dat <- read.table(file, header=T, sep="\t")
  
  # extract the columns we want and add a sample tag
  dat %>%
    select(chr, pos, context, ratio) %>%
    mutate(sample = sample_tags[tag_index]) ->
    dat
  
  # add it to our growing list
  data_list[[tag_index]] <- dat
  
  # move to the next tag (and next position in data_list)
  tag_index <- tag_index + 1
}

data_full <- bind_rows(data_list)







