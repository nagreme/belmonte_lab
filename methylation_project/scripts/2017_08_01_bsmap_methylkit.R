# Date: 2017/08/01


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




# ==========================
# --- LIBRARIES
# ==========================

library(methylKit)
library(genomation)
library(magrittr)


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
# --- DIFFERENTAIL ANALYSIS
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


# Correct for overdispersion
# sim.methylBase1<-dataSim(replicates=6,sites=1000,
#                          treatment=c(rep(1,3),rep(0,3)),
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
# )

#the parameters here are geared towards treatments (also some knowledge/idea about treatment outcome is needed?)
# sim_meth <- dataSim(replicates = 2, sites = 1000, 
#                     sample.ids = c("GLOB", "MG"), assembly = "DH12075", context = "CpG")

# 
# my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
#                                overdispersion="MN",test="Chisq",mc.cores=1)


# differential methylation
# Note: use q-value of 0.05 ************** (In calculateDiffMeth or in getMethylDiff or both?)
# It's sort of possible here, but maybe weird? (Look at documentation)git git status

# Note: vvv uses multiple cores for the first part of calculations (^ cores == ^ mem usage)
diff_CpG <- calculateDiffMeth(meth_CpG, mc.cores = 12) 
diff_CHG <- calculateDiffMeth(meth_CHG, mc.cores = 12)
diff_CHH <- calculateDiffMeth(meth_CHH, mc.cores = 12)



# Do we want these separated or together or both?
# myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
# #
# # get hypo methylated bases
# myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
# #
# #
# # get all differentially methylated bases
# myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


#diff meth per chr
# diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
diffMethPerChr(diff_CpG, plot=T, qvalue.cutoff = 0.05, meth.cutoff = 25)


# ==========================
# --- ANNOTATION
# ==========================

# Uses genomation package
# Reads bed files (might require me to convert bed formats)
               
