library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
# enableWGCNAThreads()

setwd("C:/Users/QiYan/Dropbox/AIME/Panda_HILICpos/HILIC_Controls_ExpoUnexpo/PANDA_input")
# setwd("/u/home/q/qyan/AIME/HILIC_WGCNA_Comp1")
load(file = "HILIC_control_expo_unexpo_residual_nonorm_WGCNA.RData")

datTraits <- sampleID[,-1]
rownames(datTraits) <- sampleID$SampleID

###################################
# Part I: Preprocess for WGCNA
###################################

##transpose dataset from wide to long
wide_save_residual$met <- paste("met_", c(1:nrow(wide_save_residual)),sep = "")

datExpr <- as.data.frame(t(wide_save_residual[,-c(1,2,ncol(wide_save_residual))]))
colnames(datExpr) <- wide_save_residual$met

## choose power for soft-threshold
## Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
## Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## Plot the results:
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 <- 0.9;
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###################################
# Part II: Run WGCNA
###################################

net <- blockwiseModules(datExpr, power = 12,
                       corType = "bicor",
                       maxBlockSize =5000,
                       networkType = "signed",
                       minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Comp1TOM", 
                       verbose = 3)

## open a graphics window
sizeGrWindow(12, 9)
## Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(net$colors)
## add bicor bar between metabolites and trait
cor_trait <- as.data.frame(datTraits$factorcase)
row.names(cor_trait) <- row.names(datTraits)                 ## Define the target trait (factorcase)
colnames(cor_trait) <- "factorcase"
cor_trait$factorcase[cor_trait$factorcase=="Exposed"] <- 1   ## recode to numeric
cor_trait$factorcase[cor_trait$factorcase=="Unexposed"] <- 2
GS.trait <- as.numeric(bicor(datExpr, cor_trait, robustY = FALSE, maxPOutliers = 1))
GS.traitColor <- as.vector(numbers2colors(GS.trait, signed = T))
datColors <- data.frame(mergedColors[net$blockGenes[[1]]], GS.traitColor[net$blockGenes[[1]]])
## Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], colors = datColors,
                    groupLabels = c("Module colors", "GS.trait"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "HILIC_WGCNA_feature_module.RData")

###################################
# Part III: Associate modules with traits
###################################

load(file = "HILIC_WGCNA_feature_module.RData")

## Define numbers of genes and samples
nGenes <- ncol(datExpr);
nSamples <- nrow(datExpr);
## Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
## correlation between eigengenes
signif(cor(MEs, use = "p"), 2)
## plot eigengene dendrogram
dissimME <- (1-t(cor(MEs, method = "p")))/2
hclustdatME <- hclust(as.dist(dissimME), method = "average")
plot(hclustdatME, main = "Clustering tree based of the module eigengenes")
plotEigengeneNetworks(multiME = MEs, setLabels = "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2),
                      cex.lab = 0.8)
## merge close module
merge <- mergeCloseModules(exprData = datExpr, colors = moduleColors, cutHeight = 0.25)
moduleColors.merge <- merge$colors
MEs.merge <- merge$newMEs
## plot heatmap and eigengene
for(i in 1: length(unique(moduleColors.merge))){
  which.color <- unique(moduleColors.merge)[i]
  which.ME <- MEs.merge[, paste("ME", which.color, sep = "")]
  par(mfrow = c(2,1), mar = c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[, which(moduleColors.merge == which.color)])),
          nrgcols = 30, rlabels = F, rcols = which.color,
          main = which.color, cex.main = 2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(which.ME, col = which.color, main = "", cex.main = 2,
          ylab = "eigengene expression", xlab = "array sample")
}

cor_trait.combine <- data.frame(cor_trait, datTraits$maternal_age, datTraits$maternal_raceeth, datTraits$preterm)
colnames(cor_trait.combine) <- c("factor", "maternal_age", "maternal_race", "preterm")
moduleTraitCor <- bicor(MEs.merge, cor_trait.combine, maxPOutliers = 1, robustY = FALSE);
# moduleTraitCor = cor(MEs, cor_trait, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
## Will display correlations and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
## Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(cor_trait.combine),
               yLabels = names(MEs.merge),
               ySymbols = names(MEs.merge),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

save(moduleTraitCor, moduleTraitPvalue, file = "HILIC_WGCNA_ModuleTrait_corr.RData")

## which is cor_trait
## names (colors) of the modules
modNames <- substring(names(MEs.merge), 3)

geneModuleMembership <- as.data.frame(bicor(datExpr, MEs.merge, maxPOutliers = 1));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");

geneTraitSignificance <- as.data.frame(bicor(datExpr, cor_trait.combine$factor, maxPOutliers = 1, robustY = FALSE));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) <- paste("GS.", names(cor_trait.combine)[1], sep="");
names(GSPvalue) <- paste("p.GS.", names(cor_trait.combine)[1], sep="");

## Intramodular analysis: identifying genes with high GS and MM
module <- "purple"
column <- match(module, modNames);
moduleGenes <- moduleColors.merge==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
verboseScatterplot(abs(-log(MMPvalue[moduleGenes, column])),
                   abs(-log(GSPvalue[moduleGenes, 1])),
                   xlab = paste("Module Membership P value in", module, "module"),
                   ylab = "Gene significance P value",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

save(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue, file = "HILIC_WGCNA_MMGS.RData")

###################################
# Part IV: Conduct bircor between features and traits
###################################

y <- as.vector(cor_trait$factorcase)
x <- as.matrix(datExpr)
output <- standardScreeningBinaryTrait(x, y, 
                             corFnc = bicor, 
                             corOptions = list(use = 'p',robustY=F),
                             kruskalTest = FALSE, qValues = TRUE, var.equal=FALSE, 
                             na.action="na.exclude", getAreaUnderROC = TRUE)
output <- output[order(output$pvalueStudent.1.vs.2),]

