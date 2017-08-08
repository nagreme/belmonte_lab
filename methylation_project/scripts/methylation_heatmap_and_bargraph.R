# --------------------------------------------------------------
# Date: 2017-08-02
# Nadège Pulgar-Vidal
# 
# Methylation Heatmap (and Bar Graph)(wip)
# --------------------------------------------------------------

# ============================
# --- LIBRARIES
# ============================

library(magrittr) # for piping (%>%)
library(ggplot2) # for nicer plotting functionality
library(RColorBrewer) # for ready made colour palettes
library(dplyr) # for nice table manipulation
library(gplots) # for heatmap


# ============================
# --- SETUP/INPUT
# ============================

# List of files
files <- c("/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/GLOB_meth.txt",
          "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/MG_meth.txt")

# List of samples (matching files)
sample_tags <- c("GLOB",
                "MG")

# Name and location of output heatmap
outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/output_plot.pdf"

# Bin size (absolute positions)
bin_size <- 1000


# ============================
# --- READ/FORMAT DATA
# ============================

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

# Collect data into one big table
data_full <- bind_rows(data_list)

# Put cytosines into bins and calc bin mean ratio
data_full %>%
  #mutate(sample_f = factor())
  mutate(bin = as.integer(pos/bin_size)) %>%
  group_by(sample, chr, context, bin) %>%
  summarize(avg_ratio = mean(ratio)) ->
  data_binned

# Scale ratios relative to max
data_binned$avg_ratio <- unlist(lapply(data_binned$avg_ratio, function(x) x/max(data_binned$avg_ratio)))

# Make the sample tags factors so the order of the heatmap matches theb bar chart
data_binned$sample_f = factor(data_binned$sample, levels = c("MG", "GLOB"))


# ============================
# --- VISUALIZE
# ============================

# Bar Chart of average bin ratios
graph <- ggplot(data_binned, aes(bin, avg_ratio)) +
  facet_grid(context + sample ~ chr) +
  geom_bar(stat = "identity") +
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Average bin ratio")

# Heatmap of average bin methylation ratio
graph <- ggplot(data_binned, aes(bin, sample_f)) +
  facet_grid(context ~ chr, switch = "y") +
  geom_tile(aes(fill = avg_ratio, colour = avg_ratio)) + 
  # the colour parameter above colours the tile outlines same as fill
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "")
       
ggsave(outfile_path)





