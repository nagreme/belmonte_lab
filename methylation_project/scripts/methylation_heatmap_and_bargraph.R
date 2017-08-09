# --------------------------------------------------------------
# Date: 2017-08-02
# Nadège Pulgar-Vidal
# 
# Methylation Heatmap and Bar Graph
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
outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/visualizations/output_plot.pdf"

# Bin size (absolute positions)
bin_size <- 100000


# ============================
# --- READ/FORMAT DATA
# ============================

data_list <- list()

# TODO: Find a better way to do this? Or more R-like?
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
  data_binned_full

# Scale ratios relative to max
data_binned_full$avg_ratio <- unlist(lapply(data_binned_full$avg_ratio, function(x) x/max(data_binned_full$avg_ratio)))

# Make the sample tags factors so the order of the heatmap matches theb bar chart
data_binned_full$sample_f = factor(data_binned_full$sample, levels = c("MG", "GLOB"))

# TODO: Can I write ^this^ as a mutate() call?


# --- Separate data into 2 sets by chromosome (N1-N9, N10-N19)
data_n1_9 <- filter(data_binned_full, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9"))
data_n10_19 <- filter(data_binned_full, chr %in% c("N10", "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# --- Separate CHH context so we can see it better
data_CHH <- filter(data_binned_full, context %in% c("CHH"))

# ============================
# --- VISUALIZE
# ============================

# Choose which set of data to plot
plot_data <- data_binned_full
plot_data <- data_n1_9
plot_data <- data_n10_19
plot_data <- data_CHH


# Bar Chart of average bin ratios
graph <- ggplot(plot_data, aes(bin, avg_ratio)) +
  facet_grid(context + sample ~ chr, scales = "free_x") +
  geom_bar(stat = "identity") +
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Average bin ratio")

# Heatmap of average bin methylation ratio
graph <- ggplot(plot_data, aes(x = bin, y = sample_f, fill = avg_ratio)) +
  facet_grid(context ~ chr, switch = "y", scales = "free_x") +
  geom_tile(aes(colour = avg_ratio)) + 
  # the colour parameter above colours the tile outlines same as fill
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "")  +
  #scale_fill_distiller(palette = "Spectral")

# TODO: Pick colour palette ^^^

# TODO: Should I use geom_raster instead of geom_tile? Can I? Apparently it's more efficient
       
ggsave(outfile_path, height = 7, width = 12)





