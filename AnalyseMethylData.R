###############################################################################
#  Loading the data
###############################################################################
# load packages required for analysis
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# get the 450k annotation data
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
head(ann850k)
# read in the sample sheet for the experiment               
#targets <- read.metharray.sheet("/Users/viktorian.miok/Documents/My Folder/UMFT/UMF-biochim/metilareIulie2019/", pattern="SampleSheet_forR.csv")
setwd("/Users/viktorian.miok/Documents/My Folder/UMFT/UMF-biochim/metilareIulie2019")
list.files(dataDirectory, recursive=TRUE)
dataDirectory <- "./data"
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet_forR.csv")
targets

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)

# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

###############################################################################
#  Quality control
###############################################################################

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylim = c(0,0.3),ylab="Mean detection p-values")
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group,
         pdf="qcReport-suplement.pdf")


# remove poor quality samples
keep <- colMeans(detP) < 0.2
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)
###############################################################################
## Normalisation
###############################################################################
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
###############################################################################
## Data exploration
###############################################################################
# MDS plots to look at largest sources of variation
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

###############################################################################
## Filtering
###############################################################################
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.15) == ncol(mSetSq)
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# Once the data has been filtered and normalised, it is often useful to re-examine the MDS plots
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

# The next step is to calculate M-values and beta values.
mVals <- getM(mSetSqFlt)
head(mVals[,1:7])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:7])

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
###############################################################################
## Probe-wise differential methylation analysis
###############################################################################

# this is the factor of interest
tretment <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for individual <- factor(targets$Sample_Source)
# use the above to create a design matrix
design <- model.matrix(~0+tretment, data=targets)
colnames(design) <- c(levels(tretment))#,levels(individual)[-1])
# fit the linear model
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons contMatrix <- makeContrasts(naive-rTreg,
contMatrix <- makeContrasts(Group_T-Group_N, levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix) 
fit2 <- eBayes(fit2)
# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]
DMPs <- topTable(fit2,  num=Inf, coef=1, genelist=ann850kSub)
head(DMPs)

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)


