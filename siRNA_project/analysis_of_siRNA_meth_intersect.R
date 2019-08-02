# ------------------------------------------------------------------
# Date: 2019-08-02
# Nadège Pulgar-Vidal
#
# R script to compute the average methylation ratio per siRNA locus
# from Shortstack output based on the output of intersectBed using
# the bed files of both data sets I prepared.
#
# intersectBed -wao -a 2of3_MG.bed -b MG_meth_CpG.bed MG_meth_CHG.bed MG_meth_CHH.bed -names CpG CHG CHH > MG_intersect_results.txt
#
# intersectBed -wao -a 2of3_GLOB.bed -b GLOB_meth_CpG.bed GLOB_meth_CHG.bed GLOB_meth_CHH.bed -names CpG CHG CHH > GLOB_intersect_results.txt
#
# ------------------------------------------------------------------

library(tidyverse)
library(magrittr)
