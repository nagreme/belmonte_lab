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
files <- c("/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/GLOB_test.txt",
          "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/MG_test.txt")

# List of samples (matching files)
sample_tags <- c("GLOB",
                "MG")

# Name and location of output heatmap
outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/test_heatmap.pdf"

# Bin size (absolute positions)
bin_size <- 10


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
  mutate(bin = as.integer(pos/bin_size)) %>%
  group_by(sample, chr, context, bin) %>%
  summarize(avg_ratio = mean(ratio)) ->
  data_binned

#Scale ratios relative to max
data_binned$avg_ratio <- unlist(lapply(data_binned$avg_ratio, function(x) x/max(data_binned$avg_ratio)))


# ============================
# --- VISUALIZE
# ============================

# heatmap_plot <- ggplot(data_full, aes(pos, 0)) +
#   geom_tile(aes(fill = ratio)) +
#   scale_fill_brewer(palette = "PuBu") +
  # facet_grid(context + sample ~ chr, switch = "y")

# heatmap_plot <- ggplot(data_full, aes(pos,sample, fill=ratio)) + 
#   scale_fill_brewer(palette = "PuBu") +
#   facet_grid(context + sample ~ chr, switch = "y") +
#   geom_tile()


# Bar Chart of average bin ratios
graph <- ggplot(data_binned, aes(bin, avg_ratio)) +
  facet_grid(context + sample ~ chr) +
  geom_bar(stat = "identity") +
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Average bin ratio")

# Heatmap attempt (3)

ggsave(outfile_path)





