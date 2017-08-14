# Date: 2017/08/01
# Nadège Pulgar-Vidal (nadpulgarvidal@gmail.com)

# ==========================
# --- LIBRARIES
# ==========================

library(methylKit)
library(genomation)
library(magrittr)
library(dplyr)
library(readr)


# ==========================
# --- SETUP INPUT
# ==========================

# For now just handle the files separately

glob_CpG_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/GLOB_CG.txt"
glob_CHG_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/GLOB_CHG.txt"
glob_CHH_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/GLOB_CHH.txt"
mg_CpG_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/MG_CG.txt"
mg_CHG_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/MG_CHG.txt"
mg_CHH_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/MG_CHH.txt"

files_CpG <- list(glob_CpG_file, mg_CpG_file)
files_CHG <- list(glob_CHG_file, mg_CHG_file)
files_CHH <- list(glob_CHH_file, mg_CHH_file)

# Genes only file (filtered with grep from the full gff file)
gff_annotation_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/DH12075_annotation_genes_exons.gff3"

# Annotation files
# "DH12075_annotation.gff3" (full annotation)
# "DH12075_annotation_genes.gff3" (gene features only)
# "DH12075_annotation_exons.gff3" (exon features only)
# "DH12075_annotation_genes_exons.gff3" (gene and exon features)


gff_genes_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/DH12075_annotation_genes.gff3"

flank_genes_file <- "/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/flank_genes.txt"


# ==========================
# --- READ FILES
# ==========================

#For individual objects
# glob_CpG_obj <- methRead(glob_CpG_file, sample.id="GLOB", assembly="DH12075",
#                         header=TRUE, context="CpG", resolution="base",
#                         pipeline=list(fraction=TRUE, chr.col=1, start.col=2, end.col=2,
#                                       coverage.col=7, strand.col=3, freqC.col=5 ))

# Read objects grouped by context

meth_CpG_obj <- methRead(files_CpG, 
                         sample.id=list("GLOB", "MG"), 
                         assembly="DH12075",
                         treatment = c(0,1),
                         header=TRUE, 
                         context="CpG", 
                         resolution="base",
                         pipeline=list(fraction=TRUE, chr.col=1, start.col=2, end.col=2,
                                     coverage.col=7, strand.col=3, freqC.col=5 ))

meth_CHG_obj <- methRead(files_CHG, 
                         sample.id=list("GLOB", "MG"), 
                         assembly="DH12075",
                         treatment = c(0,1),
                         header=TRUE, 
                         context="CHG", 
                         resolution="base",
                         pipeline=list(fraction=TRUE, chr.col=1, start.col=2, end.col=2,
                                       coverage.col=7, strand.col=3, freqC.col=5 ))

meth_CHH_obj <- methRead(files_CHH, 
                         sample.id=list("GLOB", "MG"), 
                         assembly="DH12075",
                         treatment = c(0,1),
                         header=TRUE, 
                         context="CHH", 
                         resolution="base",
                         pipeline=list(fraction=TRUE, chr.col=1, start.col=2, end.col=2,
                                       coverage.col=7, strand.col=3, freqC.col=5 ))



# ==========================
# --- DIFFERENTIAL ANALYSIS
# ==========================

# --- Filter based on coverage and tile the individual objects

# Note: the tileMethylCounts function does not seem to use multiple cores even when specified

meth_CpG_obj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 1000, step.size = 1000, cov.bases = 0) ->
  tiles_CpG

meth_CHG_obj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 1000, step.size = 1000, cov.bases = 0) ->
  tiles_CHG

meth_CHH_obj %>%
  filterByCoverage(lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9) %>%
  tileMethylCounts(win.size = 1000, step.size = 1000, cov.bases = 0) ->
  tiles_CHH


# --- Unite across samples (keeping contexts separate)

meth_CpG <- unite(tiles_CpG)
meth_CHG <- unite(tiles_CHG)
meth_CHH <- unite(tiles_CHH)


# --- Differential methylation
# Note: use q-value of 0.05 

# Note: vvv uses multiple cores for the first part of calculations (^ cores == ^ mem usage)
diff_CpG <- calculateDiffMeth(meth_CpG, mc.cores = 12) 
diff_CHG <- calculateDiffMeth(meth_CHG, mc.cores = 12)
diff_CHH <- calculateDiffMeth(meth_CHH, mc.cores = 12)



# Do we want these separated or together or both?


hyper_diff_CpG <- getMethylDiff(diff_CpG, difference = 25, qvalue = 0.05, type = "hyper")
hypo_diff_CpG <- getMethylDiff(diff_CpG, difference = 25, qvalue = 0.05, type = "hypo")
all_diff_CpG25 <- getMethylDiff(diff_CpG, difference = 25, qvalue = 0.05, type = "all")

hyper_diff_CHG25 <- getMethylDiff(diff_CHG, difference = 25, qvalue = 0.05, type = "hyper")
hypo_diff_CHG25 <- getMethylDiff(diff_CHG, difference = 25, qvalue = 0.05, type = "hypo")
all_diff_CHG25 <- getMethylDiff(diff_CHG, difference = 25, qvalue = 0.05, type = "all")

hyper_diff_CHH25 <- getMethylDiff(diff_CHH, difference = 25, qvalue = 0.05, type = "hyper")
hypo_diff_CHH25 <- getMethylDiff(diff_CHH, difference = 25, qvalue = 0.05, type = "hypo")
all_diff_CHH25 <- getMethylDiff(diff_CHH, difference = 25, qvalue = 0.05, type = "all")




#diff meth per chr
# diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(diff_CpG, plot=T, qvalue.cutoff = 0.05, meth.cutoff = 25)
diffMethPerChr(diff_CHG, plot=T, qvalue.cutoff = 0.05, meth.cutoff = 25)
diffMethPerChr(diff_CHH, plot=T, qvalue.cutoff = 0.05, meth.cutoff = 25)





# ==========================
# --- ANNOTATION
# ==========================

# Uses genomation package
# Can read gff/gtf files instead of bed files

# Not right now
# # --- Get feature ranges information
# gff_GRanges_obj <- gffToGRanges(gff_annotation_file)
# 
# 
# # --- Annotate with gene parts
# # example: annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
# annotateWithGeneParts(as(diff_CpG, "GRanges"), GRangesList(gff_GRanges_obj))
# #the second param needs to be a GRangesList instead of a GRanges object for this to work
# 


# --- Extract needed columns from gene gff to a file
read.table(gff_genes_file, header=F, sep="\t", 
           col.names = c("chr", "src", "feature_type", "start", "end", "myst", "strand", "myst2", "name")) %>%
  select(chr, start, end, name) %>%
  write_delim(flank_genes_file, delim="\t")
# Note: If you find a way to suppress writing a header, add it here, otherwise delete manually before next step

gene_names <- as.data.frame(dat$name)


# --- Read feature flanks
# cpg.obj=readFeatureFlank(CpGisland, feature.flank.name=c("CpGi","shores"))
feat_flanks <- readFeatureFlank(flank_genes_file, 
                                feature.flank.name = c("gene", "flank"), 
                                flank = 3000)

# Add the gene names manually because readFeatureFlank Ignores them
mcols(feat_flanks[[1]]) <- gene_names


# --- Annotate feature flanks
# diffCpGann=annotateWithFeatureFlank(as(myDiffCpG25,"GRanges"),
#                                     cpg.obj$CpGi, cpg.obj$shores,
#                                     feature.name="CpGi",flank.name="shores")

# Doesn't really work for us -> no gene names???
diff_CpG_flank <- annotateWithFeatureFlank(target = as(all_diff_CpG25, "GRanges"),
                                           feature = feat_flanks$gene, feature.name = "Genes",
                                           flank = feat_flanks$flank, flank.name = "Flanks")

diff_CHG_flank <- annotateWithFeatureFlank(target = as(all_diff_CHG25, "GRanges"),
                                           feature = feat_flanks$gene, feature.name = "Genes",
                                           flank = feat_flanks$flank, flank.name = "Flanks")

diff_CHH_flank <- annotateWithFeatureFlank(target = as(all_diff_CHH25, "GRanges"),
                                           feature = feat_flanks$gene, feature.name = "Genes",
                                           flank = feat_flanks$flank, flank.name = "Flanks")

# Choose which data to plot
diff_flank_dat <- diff_CpG_flank
diff_flank_dat <- diff_CHG_flank
diff_flank_dat <- diff_CHH_flank

# The pie chart they provide
genomation::plotTargetAnnotation(diff_flank_dat, col=c("black","gray","white"),
                                 main="differential methylation annotation")

getFeatsWithTargetsStats(diff_flank_dat, percentage = TRUE)






# selectByOverlap(all_diff_CpG25, ranges = feat_flanks[[1]])

GenomicRanges::findOverlaps(as(all_diff_CpG25, "GRanges"), feat_flanks[[1]])





# ==========================
# --- METHODS NOT IN USE
# ==========================

# --- Tradiotnal way of reading in (meant for files from AMP pipeline?)

# file.listCG <- list("/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/GLOB_CG.txt",
#                     "/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/MG_CG.txt")
# 
# file.listCHG <- list("/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/GLOB_CHG.txt",
#                     "/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/MG_CHG.txt")
# 
# file.listCHH <- list("/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/GLOB_CHH.txt",
#                     "/media/deirdre5tb/methylome/bsmap-2.90/bsmap_full_run/MG_CHH.txt")
# 
# 
# # read the files to a methylRawList object
# methReadObjCG <- methRead(file.list,
#                sample.id=list("GLOB","MG"),
#                assembly="DH12075",
#                treatment=c(1,0),
#                context="CG")
# 
# 
# # read the files to a methylRawList object
# methReadObjCHG <- methRead(file.list,
#                           sample.id=list("GLOB","MG"),
#                           assembly="DH12075",
#                           treatment=c(1,0),
#                           context="CHG")
# 
# 
# # read the files to a methylRawList object
# methReadObjCHH <- methRead(file.list,
#                           sample.id=list("GLOB","MG"),
#                           assembly="DH12075",
#                           treatment=c(1,0),
#                           context="CHH")



# (Print) Descriptive Statistics
# getMethylationStats(methReadObjCG[[1]],plot=FALSE,both.strands=FALSE)
# getMethylationStats(methReadObjCHG[[1]],plot=FALSE,both.strands=FALSE)
# getMethylationStats(methReadObjCHH[[1]],plot=FALSE,both.strands=FALSE)


# --- Way of reading in Bismark Aligner data (sorted SAM or BAM files)

# But our .bam files are >20G x6...

# my.methRaw=processBismarkAln( location = 
#                                 system.file("extdata",
#                                             "test.fastq_bismark.sorted.min.sam", 
#                                             package = "methylKit"),
#                               sample.id="test1", assembly="hg18", 
#                               read.context="CpG", save.folder=getwd())



# --- Way to read in BSMAP files

# methobj_CpG <- methRead("/home/mark/Documents/Nadege/belmonte_lab/methylation_project/data/methylation_calls_bsmap/separated_contexts/GLOB_CG.txt",
#                sample.id="GLOB",
#                assembly="DH12075",
#                header=TRUE,
#                context="CpG",
#                resolution="base",
#                pipeline=list(fraction=TRUE,
#                              chr.col=1,
#                              start.col=2,
#                              end.col=2,
#                              coverage.col=7,
#                              strand.col=3,
#                              freqC.col=5 ))
# Note: The columns chosen for the pipeline list parameters are the right ones for the current version of bsmap (May 2017)

# --- Descriptive stats 
# Used as a test to verify the reading
# getMethylationStats(methobj_CpG, plot=TRUE, both.strands=FALSE)
# getCoverageStats(methobj_CpG, plot = TRUE, both.strands = FALSE)



