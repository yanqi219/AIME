library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
# enableWGCNAThreads()

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
# setwd("/u/home/q/qyan/AIME/HILIC_WGCNA_Comp1")
load(file = "HILIC_control_expo_unexpo_residual_WGCNA.RData")

datTraits <- sampleID[,-1]
rownames(datTraits) <- sampleID$SampleID

# Preprocess for WGCNA
##transpose dataset from wide to long
wide_save_residual$met <- paste("met_", c(1:nrow(wide_save_residual)),sep = "")

datExpr <- as.data.frame(t(wide_save_residual[,-c(1,2,ncol(wide_save_residual))]))
colnames(datExpr) <- wide_save_residual$met

## choose power for soft-threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Run WGCNA

net = blockwiseModules(datExpr, power = 2,
                       corType = "bicor",
                       maxBlockSize =5000,
                       networkType = "signed",
                       minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Comp1TOM", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "HILIC_WGCNA_feature_module.RData")
