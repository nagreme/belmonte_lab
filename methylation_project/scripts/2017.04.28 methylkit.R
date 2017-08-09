save.image(file="20170519methylkit.RData")
load("20170518methylkit.RData")
setwd("/Volumes/Deirdre2TB/methylKit")

setwd("/Volumes/Deirdre2TB/")
load("/media/Deirdre2TB/methylKit/extdata/20170508methylkit.RData")


library(methylKit)
#read in files from the whole cytosine report - do not need the system.file
file.list=list( "/Volumes/Deirdre2TB/GLOB1CX.txt",
                "/Volumes/Deirdre2TB/GLOB2CX.txt",
                "/Volumes/Deirdre2TB/GLOB3CX.txt",
                "/Volumes/Deirdre2TB/MG1CX.txt",
                "/Volumes/Deirdre2TB/MG2CX.txt",
                "/Volumes/Deirdre2TB/MG3CX.txt")

#FOR CpG METHYLATION

CpG.obj=methRead(file.list,
         sample.id=list("GLOB1","GLOB2","GLOB3","MG1","MG2","MG3"),
         assembly="DH12075_DK",
         dbtype = NA, 
         pipeline = "bismarkCytosineReport",
         header = FALSE, 
         skip = 0, 
         sep = "\t", 
         context = "CpG",
         resolution = "base",
         treatment = c(0,0,0,1,1,1),
         dbdir = getwd(), 
         mincov = 20)

#getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
filtered.myobjCpG=filterByCoverage(myobjCpG,lo.count=30,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
CpG.GLOB1=getMethylationStats(filtered.myobjCpG[[1]],plot=TRUE,both.strands=FALSE)
CpG.GLOB2=getMethylationStats(filtered.myobjCpG[[2]],plot=FALSE,both.strands=FALSE)
CpG.GLOB3=getMethylationStats(filtered.myobjCpG[[3]],plot=FALSE,both.strands=FALSE)
CpG.MG1=getMethylationStats(filtered.myobjCpG[[4]],plot=FALSE,both.strands=FALSE)
CpG.MG2=getMethylationStats(filtered.myobjCpG[[5]],plot=FALSE,both.strands=FALSE)
CpG.MG3=getMethylationStats(filtered.myobjCpG[[6]],plot=FALSE,both.strands=FALSE)

write.table(CpG.GLOB1,"CpG.GLOB1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpG.GLOB2,"CpG.GLOB2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpG.GLOB3,"CpG.GLOB3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpG.mg1,"CpG.mg1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpG.mg2,"CpG.mg2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpG.mg3,"CpG.mg3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

#samples must be united to one before proceeding with analysis
methCpG=unite(filtered.myobjCpG, destrand=FALSE)
# #get that pretty Pearson correlation plot
# getCorrelation(methCpG,plot=TRUE)
# #hierarchical clustering, tree
# clusterSamples(methCpG, dist="correlation", method="ward", plot=TRUE)
# #get a PCA, probably more useful for this kind of analysis
# PCASamples(methCpG)
# #NOTE - myDiff=calculateDiffMeth(meth,num.cores=2) more cores = faster yo
myDiffCpG=calculateDiffMeth(methCpG)
#get hypermethylated bases
myDiffCpGhyper25.hyper=getMethylDiff(myDiffCpG,difference=25,qvalue=0.01,type="hyper")
#get hypomethylated bases
myDiffCpGhypo25=getMethylDiff(myDiffCpG,difference=25,qvalue=0.01,type="hypo")
#get ALL differentially methylated bases
myDiffCpG25=getMethylDiff(myDiffCpG,difference=25,qvalue=0.01)
#also if you want a plot, per chromosome, of differential methylation
chrdiff=diffMethPerChr(myDiffCpG,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)


#Similarly, we can read the CpG island annotation and annotate our differentially methylated bases/regions with them.

# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
library(genomation)
CpGisland="/Volumes/Deirdre2TB/methylKit/extdata/DH12075cpg.txt"
cpg.obj=readFeatureFlank(CpGisland, feature.flank.name=c("CpGi","shores"))

# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiffCpG25,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="differential methylation annotation")
getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)



#now we annotate! the gff file has been converted to a BED file
library(genomation)
gene.obj=readTranscriptFeatures("/Volumes/Deirdre2TB/methylKit/DH12075DKbednew.txt")
annotateWithGeneParts(as(myDiffCpG25,"GRanges"),gene.obj)
promotersCpG=regionCounts(myobjCpG,gene.obj$promoters)
head(promoters[[1]])

#or does it? can we use the gff as an annotation file after all?
GFFfile = "/Volumes/Deirdre2TB/Bna_genes_v3_1.gff3"
gff = gffToGRanges(GFFfile)
head(gff)

diffCpGann=annotateWithGeneParts(myDiffCpG25,gff)


#After getting the annotation of differentially methylated regions, we can get the distance to TSS and nearest 
#gene name using the  getAssociationWithTSS function from genomation package.
diffAnnCpG=annotateWithGeneParts(as(myDiff25CpG25,"GRanges"),gff)
head(getAssociationWithTSS(diffAnnCpG))

#It is also desirable to get percentage/number of differentially methylated regions that overlap with 
#intron/exon/promoters
getTargetAnnotationStats(diffAnnCpG,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnnCpG,precedence=TRUE,
                     main="differential methylation annotation")

getFeatsWithTargetsStats(diffAnnCpG,percentage=TRUE)



#For some situations, it might be desirable to summarize methylation information over tiling windows rather 
#than doing base-pair resolution analysis. methylKit provides functionality to do such analysis. 
#The function below tiles the genome with windows 1000bp length and 1000bp step-size and summarizes the 
#methylation information on those tiles. In this case, it returns a  methylRawList object which can be fed into 
#unite and calculateDiffMeth functions consecutively to get differentially methylated regions. The tilling 
#function adds up C and T counts from each covered cytosine and returns a total C and T count for each tile.
#tiles=tileMethylCounts(myobj,win.size=50,step.size=50)
#head(tiles[[1]],3)


#this was for running in windows
#DO FOR CHG METHYLATION

myobjCHG=methRead(file.list,
                  sample.id=list("GLOB1","GLOB2","GLOB3","MG1","MG2","MG3"),
                  assembly="DH12075_DK",
                  dbtype = NA, 
                  pipeline = "bismarkCytosineReport",
                  header = FALSE, 
                  skip = 0, 
                  sep = "\t", 
                  context = "CHG",
                  resolution = "base",
                  treatment = c(0,0,0,1,1,1),
                  dbdir = getwd(), 
                  mincov = 20)
#getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
filtered.myobjCHG=filterByCoverage(myobjCHG,lo.count=30,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)

CHG.GLOB1=getMethylationStats(filtered.myobjCHG[[1]],plot=FALSE,both.strands=FALSE)
CHG.GLOB2=getMethylationStats(filtered.myobjCHG[[2]],plot=FALSE,both.strands=FALSE)
CHG.GLOB3=getMethylationStats(filtered.myobjCHG[[3]],plot=FALSE,both.strands=FALSE)
CHG.MG1=getMethylationStats(filtered.myobjCHG[[4]],plot=FALSE,both.strands=FALSE)
CHG.MG2=getMethylationStats(filtered.myobjCHG[[5]],plot=FALSE,both.strands=FALSE)
CHG.MG3=getMethylationStats(filtered.myobjCHG[[6]],plot=FALSE,both.strands=FALSE)

write.table(CHG.GLOB1,"CHG.GLOB1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHG.GLOB2,"CHG.GLOB2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHG.GLOB3,"CHG.GLOB3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHG.mg1,"CHG.mg1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHG.mg2,"CHG.mg2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHG.mg3,"CHG.mg3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)


#samples must be united to one before proceeding with analysis
methCHG=unite(filtered.myobjCHG, destrand=FALSE)
#get that pretty Pearson correlation plot
getCorrelation(methCHG,plot=TRUE)
#hierarchical clustering, tree
clusterSamples(methCHG, dist="correlation", method="ward", plot=TRUE)
#get a PCA, probably more useful for this kind of analysis
PCASamples(methCHG)

#NOTE - myDiff=calculateDiffMeth(meth,num.cores=2) more cores = faster yo
myDiffCHG=calculateDiffMeth(methCHG)
#get hypermethylated bases
myDiffCHGhyper25.hyper=getMethylDiff(myDiffCHG,difference=25,qvalue=0.01,type="hyper")
#get hypomethylated bases
myDiffCHGhypo25=getMethylDiff(myDiffCHG,difference=25,qvalue=0.01,type="hypo")
#get ALL differentially methylated bases
myDiffCHG25=getMethylDiff(myDiffCHG,difference=25,qvalue=0.01)
#also if you want a plot, per chromosome, of differential methylation
diffMethPerChr(myDiffCHG,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)

chrdiffCHG=diffMethPerChr(myDiffCpG,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)




#now we annotate! the gff file has been converted to a BED file
library(genomation)
gene.obj=readTranscriptFeatures("E:/methylKit/DH12075.bed")
annotateWithGeneParts(as(myDiffCHG25,"GRanges"),gene.obj)


promotersCHG=regionCounts(myobjCHG,gene.obj$promoters)
head(promoters[[1]])

#After getting the annotation of differentially methylated regions, we can get the distance to TSS and nearest 
#gene name using the  getAssociationWithTSS function from genomation package.
diffAnnCHG=annotateWithGeneParts(as(myDiff25CHG25,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnnCHG))

#It is also desirable to get percentage/number of differentially methylated regions that overlap with 
#intron/exon/promoters
getTargetAnnotationStats(diffAnnCHG,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnnCHG,precedence=TRUE,
                     main="differential methylation annotation")

getFeatsWithTargetsStats(diffAnnCHG,percentage=TRUE)




#DO FOR CHH METHYLATION

file.list=list( "/Volumes/Deirdre2TB/GLOB1CX.txt",
                "/Volumes/Deirdre2TB/GLOB2CX.txt",
                "/Volumes/Deirdre2TB/GLOB3CX.txt",
                "/Volumes/Deirdre2TB/MG1CX.txt",
                "/Volumes/Deirdre2TB/MG2CX.txt",
                "/Volumes/Deirdre2TB/MG3CX.txt")
myobjCHH=methRead(file.list,
                  sample.id=list("GLOB1","GLOB2","GLOB3","MG1","MG2","MG3"),
                  assembly="DH12075_DK",
                  dbtype = NA, 
                  pipeline = "bismarkCytosineReport",
                  header = FALSE, 
                  skip = 0, 
                  sep = "\t", 
                  context = "CHH",
                  resolution = "base",
                  treatment = c(0,0,0,1,1,1),
                  dbdir = getwd(), 
                  mincov = 20)
#getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
filtered.myobjCHH=filterByCoverage(myobjCHH,lo.count=30,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)

CHH.GLOB1=getMethylationStats(filtered.myobjCHH[[1]],plot=FALSE,both.strands=FALSE)
CHH.GLOB2=getMethylationStats(filtered.myobjCHH[[2]],plot=FALSE,both.strands=FALSE)
CHH.GLOB3=getMethylationStats(filtered.myobjCHH[[3]],plot=FALSE,both.strands=FALSE)
CHH.MG1=getMethylationStats(filtered.myobjCHH[[4]],plot=FALSE,both.strands=FALSE)
CHH.MG2=getMethylationStats(filtered.myobjCHH[[5]],plot=FALSE,both.strands=FALSE)
CHH.MG3=getMethylationStats(filtered.myobjCHH[[6]],plot=FALSE,both.strands=FALSE)

write.table(CHH.GLOB1,"CHH.GLOB1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHH.GLOB2,"CHH.GLOB2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHH.GLOB3,"CHH.GLOB3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHH.mg1,"CHH.mg1.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHH.mg2,"CHH.mg2.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CHH.mg3,"CHH.mg3.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

#samples must be united to one before proceeding with analysis
methCHH=unite(filtered.myobjCHH, destrand=FALSE)
#get that pretty Pearson correlation plot
getCorrelation(methCHH,plot=TRUE)
#hierarchical clustering, tree
clusterSamples(methCHH, dist="correlation", method="ward", plot=TRUE)
#get a PCA, probably more useful for this kind of analysis
PCASamples(methCHH)

#NOTE - myDiff=calculateDiffMeth(meth,num.cores=2) more cores = faster yo
myDiffCHH=calculateDiffMeth(methCHH,mc.cores=12)
#get hypermethylated bases
myDiffCHHhyper25.hyper=getMethylDiff(myDiffCHH,difference=25,qvalue=0.01,type="hyper")
#get hypomethylated bases
myDiffCHHhypo25=getMethylDiff(myDiffCHH,difference=25,qvalue=0.01,type="hypo")
#get ALL differentially methylated bases
myDiffCHH25=getMethylDiff(myDiffCHH,difference=25,qvalue=0.01)
#also if you want a plot, per chromosome, of differential methylation
diffMethPerChr(myDiffCHH,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)

chrdiffCHH=diffMethPerChr(myDiffCpG,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

#now we annotate! the gff file has been converted to a BED file
library(genomation)
gene.obj=readTranscriptFeatures("E:/methylKit/DH12075.bed")
annotateWithGeneParts(as(myDiffCHH25,"GRanges"),gene.obj)


promotersCHH=regionCounts(myobjCHH,gene.obj$promoters)
head(promoters[[1]])

#After getting the annotation of differentially methylated regions, we can get the distance to TSS and nearest 
#gene name using the  getAssociationWithTSS function from genomation package.
diffAnnCHH=annotateWithGeneParts(as(myDiff25CHH25,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnnCHH))

#It is also desirable to get percentage/number of differentially methylated regions that overlap with 
#intron/exon/promoters
getTargetAnnotationStats(diffAnnCHH,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnnCHH,precedence=TRUE,
                     main="differential methylation annotation")

getFeatsWithTargetsStats(diffAnnCHH,percentage=TRUE)




#returns a  methylRawList object which can be fed into unite and calculateDiffMeth functions consecutively to get 
#differentially methylated regions. The tilling function adds up C and T counts from each covered cytosine and returns
#a total C and T count for each tile.

#this can be done for any C context
CpGwindow=tileMethylCounts(CpG.obj,
                 #integer for size of tiling windows
                 win.size=1000,
                 #integer for step size of tiling windows
                 step.size=1000,
                 #number of bases to be covered in a given window
                 cov.bases=0,
                 mc.cores=1,
                 #logical determineing whether stored as a flat file database
                 save.db=FALSE)
#unite samples prior to further analysis into a single data entity
CpGwindowmeth=unite(CpGwindow, destrand=FALSE)
#this makes the pretty pearson correlation plot
getCorrelation(CpGwindowmeth,plot=TRUE)

PCASamples(CpGwindowmeth)
CpGwindowdiff=calculateDiffMeth(CpGwindowmeth)
CpGwindodiffHYPER=getMethylDiff(CpGwindowdiff,difference=25,qvalue=0.01,type="hyper")
CpGwindodiffHYPO=getMethylDiff(CpGwindowdiff,difference=25,qvalue=0.01,type="hypo")
CpGwindodiffALL=getMethylDiff(CpGwindowdiff,difference=25,qvalue=0.01)
diffMethPerChr(CpGwindowdiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

write.table(CpGwindodiffHYPER,"CpGwindowHYPER.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpGwindodiffHYPO,"CpGwindowHYPO.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(CpGwindodiffALL,"CpGwindowALL.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

head(CpGwindodiffALL)

#for CpG islands at least, these ones are 300bp long I believe
CpGisland="/Volumes/Deirdre2TB/methylKit/extdata/DH12075cpg.txt"
cpg.obj=readFeatureFlank(CpGisland, feature.flank.name=c("CpGi","shores"))

library(GenomicRanges)

# convert methylDiff object to GRanges and annotate
diffCpGwindowann=annotateWithFeatureFlank(as(CpGwindodiffALL,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
plotTargetAnnotation(diffCpGwindowann,col=c("green","gray","white"),
                     main="differential methylation annotation")


getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)


#Since we are working with highly heterogenous samples it is likely that there is heavy overdispersion.
#by default correction for overdispersion is off in the methylKit package
#NOTE - myDiff=calculateDiffMeth(meth,num.cores=2) more cores = faster yo
ODDiffCpG=calculateDiffMeth(methCpG,
                            overdispersion="MN",
                            num.cores=10)
#get hypermethylated bases
ODDiffCpGhyper25=getMethylDiff(ODDiffCpG,difference=25,qvalue=0.01,type="hyper")
#get hypomethylated bases
ODDiffCpGhypo25=getMethylDiff(ODDiffCpG,difference=25,qvalue=0.01,type="hypo")
#get ALL differentially methylated bases
ODDiffCpG25=getMethylDiff(ODDiffCpG,difference=25,qvalue=0.01)
#also if you want a plot, per chromosome, of differential methylation
chrODdiff=diffMethPerChr(ODDiffCpG,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)

write.table(ODDiffCpGhyper25,"CpGodHYPER.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(ODDiffCpGhypo25,"CpGodHYPO.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(ODDiffCpG25,"CpGodALL.txt",append=FALSE,
            quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)


