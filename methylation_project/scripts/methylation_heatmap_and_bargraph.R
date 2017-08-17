# --------------------------------------------------------------
# Date: 2017-08-02
# Nadège Pulgar-Vidal
# 
# Methylation Heatmap and Bar Graph
# And transposable element density!
# And gene density!
# And CpG island density!
#
# Note: You can't just run the entire content of the script: it
# won't work for everything. Go through section by section and 
# make sure to adjust the input/setup data and select the options 
# you want along the way
# --------------------------------------------------------------

# ============================
# --- LIBRARIES
# ============================

library(magrittr) # for piping (%>%)
library(ggplot2) # for nicer plotting functionality
library(RColorBrewer) # for ready made colour palettes
library(dplyr) # for nice table manipulation
library(readr) #for better reading of tabulated input


# ============================
# --- SETUP
# ============================

# Bin size (absolute positions)
bin_size <- 100000

# Name and location of output file
outfile_path <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/visualizations/output_plot.pdf"


# ***** Methlyation Data *****

# List of input files
files <- c("/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/GLOB_meth.txt",
          "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/MG_meth.txt")

# List of sample names (matching files)
sample_tags <- c("GLOB",
                "MG")


# ***** Transposable Elements *****

# Transposable elements file
trans_elem_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/DH12075_trans_elem.gff"
# I had to remove the comment lines with the gff version and genome name at the top of 
# the file for the reading to work


# ***** Genes *****

# Gene file
genes_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/annotation/DH12075_annotation_genes_no_scaffold.gff3'


# ***** CpG Islands *****
CpG_island_file <- '/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/DH12075_CpG_islands.txt'


# ============================
# --- READ/FORMAT DATA
# ============================

# ***** Methlyation Data *****

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

# For reordering the chr names
data_binned_full$chr_order <- factor(data_binned_full$chr, levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                                                      "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# Separate data into 2 sets by chromosome (N1-N9, N10-N19)
data_n1_10 <- filter(data_binned_full, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"))
data_n11_19 <- filter(data_binned_full, chr %in% c("N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# Separate CHH context so we can see it better
data_CHH <- filter(data_binned_full, context %in% c("CHH"))



# ***** Transposable Elements *****

read.table(trans_elem_file, header=F, sep="\t",
           col.names = c("chr", "src", "just", "start" ,"end", 
                         "score", "strand", "myst", "info")) %>%
  select(chr, start, end) -> #we only need chr and start for plotting purposes but I'm keeping end just in case for now 
  data_trans_elem

# Put the data into bins
data_trans_elem %>%
  mutate(bin = as.integer(start/bin_size)) %>%
  group_by(chr, bin) %>%
  summarise(bin_count = n()) ->
  data_trans_elem_binned

# Scale bin counts so scale is consistent when splitting A and C genome
data_trans_elem_binned$bin_count <- unlist(lapply(data_trans_elem_binned$bin_count, function(x) x/max(data_trans_elem_binned$bin_count)))


# For reordering the chr names
data_trans_elem_binned$chr_order <- factor(data_trans_elem_binned$chr, 
                                           levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                                      "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# Separate data into 2 sets by chromosome (N1-N9, N10-N19)
data_trans_elem_n1_10 <- filter(data_trans_elem_binned, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"))
data_trans_elem_n11_19 <- filter(data_trans_elem_binned, chr %in% c("N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))



# ***** Genes  *****

read.table(genes_file, header=F, sep="\t",
           col.names = c("chr", "src", "feat_type", "start" ,"end", 
                         "score", "strand", "myst", "info")) %>%
  select(chr, start, end) -> #we only need chr and start for plotting purposes but I'm keeping end just in case for now 
  data_genes

# Put the data into bins
data_genes %>%
  mutate(bin = as.integer(start/bin_size)) %>%
  group_by(chr, bin) %>%
  summarise(bin_count = n()) ->
  data_genes_binned

# Scale bin counts so scale is consistent when splitting A and C genome
data_genes_binned$bin_count <- unlist(lapply(data_genes_binned$bin_count, function(x) x/max(data_genes_binned$bin_count)))


# For reordering the chr names
data_genes_binned$chr_order <- factor(data_genes_binned$chr, 
                                           levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                                      "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# Separate data into 2 sets by chromosome (N1-N9, N10-N19)
data_genes_n1_10 <- filter(data_genes_binned, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"))
data_genes_n11_19 <- filter(data_genes_binned, chr %in% c("N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))



# ***** CpG Islands *****

read.table(CpG_island_file, header=F, sep="\t",
           col.names = c("chr", "start", "end", "label")) %>%
  select(chr, start, end) -> #we only need chr and start for plotting purposes but I'm keeping end just in case for now 
  data_CpG_islands

# Put the data into bins
data_CpG_islands %>%
  mutate(bin = as.integer(start/bin_size)) %>%
  group_by(chr, bin) %>%
  summarise(bin_count = n()) ->
  data_CpG_islands_binned

# Scale bin counts so scale is consistent when splitting A and C genome
data_CpG_islands_binned$bin_count <- unlist(lapply(data_CpG_islands_binned$bin_count, function(x) x/max(data_CpG_islands_binned$bin_count)))


# For reordering the chr names
data_CpG_islands_binned$chr_order <- factor(data_CpG_islands_binned$chr, 
                                      levels = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                                                 "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))

# Separate data into 2 sets by chromosome (N1-N9, N10-N19)
data_CpG_islands_n1_10 <- filter(data_CpG_islands_binned, chr %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"))
data_CpG_islands_n11_19 <- filter(data_CpG_islands_binned, chr %in% c("N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19"))


# ============================
# --- VISUALIZE
# ============================

# *** Methylation Data ***

# Choose which set of data to plot
plot_data <- data_binned_full
plot_data <- data_n1_10
plot_data <- data_n11_19
plot_data <- data_CHH


# Bar Chart of average bin ratios
graph <- ggplot(plot_data, aes(bin, avg_ratio)) +
  facet_grid(context + sample ~ chr_order, scales = "free_x") +
  geom_bar(stat = "identity", aes(fill = avg_ratio, colour = avg_ratio)) +
  guides(fill=FALSE) + 
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Average bin ratio") +
  # scale_color_gradient2(midpoint = 0.4, # for full
  scale_color_gradient2(midpoint = 0.05,# for CHH
                        # low="#C5F6D5", 
                        # mid="#00517C", 
                        # high="#021637") # reverse sad snape
                        low="#021637",
                        mid="#00517C",
                        high="#C5F6D5") #sad snape
                                                      

ggsave(outfile_path, height = 3, width = 20) #use width 20ish for full, 12 for one


sad_snape_colours <- colorRampPalette(c("#021637","#00517C", "#C5F6D5"))(100)

# Heatmap of average bin methylation ratio
graph <- ggplot(plot_data, aes(x = bin, y = sample_f) +
  facet_grid(context ~ chr_order, switch = "y", scales = "free_x") +
  geom_tile(aes(fill = avg_ratio), colour = avg_ratio)) + 
  # the colour parameter above colours the tile outlines same as fill
  labs(title = "Methylation ratio",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "") #+
# scale_fill_manual(values = sad_snape_colours)
  # scale_fill_distiller(palette = "Spectral")

# TODO: Figure out how to change the heatmap colours (I can get another colour scale but the tiles don't change)

# TODO: Should I use geom_raster instead of geom_tile? Can I? Apparently it's more efficient
       
ggsave(outfile_path, height = 7, width = 20)



# ***** Transposable Elements *****

plot_data <- na.omit(data_trans_elem_binned) # to exclude scaffolds
plot_data <- data_trans_elem_n1_10
plot_data <- data_trans_elem_n11_19

# Note: The chr_order col contains NAs if there are scaffolds in the input gff

# Heatmap
graph <- ggplot(plot_data, aes(x = bin, y="bin_count")) +
  geom_tile(aes(fill=bin_count,  colour = bin_count)) +
  facet_grid(~ chr_order, scales = "free_x") +
  labs(title = "Transposable Element Density",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "")

ggsave(outfile_path, height = 2, width = 12)

# Man these heatmaps are not that great and I can't get them to behave like I want them to
# I just want a barchart but the count is represented by tile colour rather than bar height...


# Bar Chart
graph <- ggplot(plot_data, aes(bin, bin_count)) +
  facet_grid(~ chr_order, scales = "free_x") +
  geom_bar(stat = "identity", aes(fill = bin_count, colour = bin_count)) +
  guides(fill=FALSE) + 
  scale_y_continuous(limits = c(0,1)) +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Transposable Element Density",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Bin density") +
  scale_color_gradient2(low="#FFFFFF", 
                        midpoint = 0.5,
                        mid="#9B7BC0",
                        high="#473E75",
                        limits = c(0,1)) #white low, dark purple high


ggsave(outfile_path, height = 2, width = 20) 



# ***** Genes  *****

plot_data <- data_genes_binned
plot_data <- data_genes_n1_10
plot_data <- data_genes_n11_19

# Bar Chart
graph <- ggplot(plot_data, aes(bin, bin_count)) +
  facet_grid(~ chr_order, scales = "free_x") +
  geom_bar(stat = "identity", aes(fill = bin_count, colour = bin_count)) +
  guides(fill=FALSE) + 
  scale_y_continuous(limits = c(0,1)) +
  theme(panel.background = element_rect(fill = "white")) +
  labs(title = "Gene Density",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Bin density") +
  scale_color_gradient2(low="#FFFFFF", 
                       midpoint = 0.7,
                       mid="#9B7BC0",
                       high="#473E75",
                       limits = c(0,1)) #white low, dark purple high


ggsave(outfile_path, height = 2, width = 20) 



# ***** CpG Islands *****

plot_data <- na.omit(data_CpG_islands_binned) # to exclude scaffolds
plot_data <- data_CpG_islands_n1_10
plot_data <- data_CpG_islands_n11_19

# Bar Chart
graph <- ggplot(plot_data, aes(bin, bin_count)) +
  facet_grid(~ chr_order, scales = "free_x") +
  geom_bar(stat = "identity", aes(fill = bin_count, colour = bin_count)) +
  guides(fill=FALSE) + 
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "CpG Island Density",
       x = paste0("Bin number\n(",bin_size,"bp bins)"),
       y = "Bin density") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_color_gradient2(low="#FFFFFF", 
                        midpoint = 0.5,
                        mid="#9B7BC0",
                        high="#473E75",
                        limits = c(0,1)) #white low, dark purple high


ggsave(outfile_path, height = 2, width = 20) 



#save.image('/home/mark/Documents/Nadege/belmonte_lab/methylation_project/Rdata/.RData_2017_08_16_methylation_heatmap_barplot')
