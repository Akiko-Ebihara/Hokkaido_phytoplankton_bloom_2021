###Library packages####
library(tidyverse)
library(magrittr)
library(vegan)
library(viridis)
library(cowplot)
library(ggsignif)
library(ggrepel)
library(scales)
library(ggbiplot)
library(ggforce)
library(corrplot)
library(psych)
library(Hmisc) #correlation calculate
library(spaa) #correlation test
library(compositions) #Centered log ratio (clr) transformation
library(zCompositions) #Replace zero
#library(RColorBrewer)

####Alpha-, beta-diversity####
###packages for Alpha-, beta-diversity
library(vegan)

####Co-occurrence network analysis####
library(WGCNA)
library(phyloseq)
library(ComplexHeatmap)

###plot
library(corrplot)
#####

###Co-occurrece analysis using WGCNA####
##Data to be used
#1.S16_Cruise.CLR.df
#2.ASV_16S.Phylum.RE
#3.S18_Cruise.CLR.df
#4.ASV_18S.Phylum.RE
#5.Env.data

###Make data frame
##Marge 16S and 18S clr-transformed data
dim(S16_Cruise.CLR.df) #90samples, 818 ASVs 99 834
dim(S18_Cruise.CLR.df) #60samples, 692 ASVs  66 680
#Total 1510 ASVs
S16_Marge <- S16_Cruise.CLR.df %>%
  dplyr::mutate(SampleID = rownames(.))
S18_Marge <- S18_Cruise.CLR.df %>%
  dplyr::mutate(SampleID = rownames(.))
Marge_CLR <- left_join(S16_Marge, S18_Marge, by = "SampleID") %>%
  dplyr::select(SampleID, everything())
rownames(Marge_CLR) <- Marge_CLR$SampleID
Marge_CLR %<>% dplyr::select(-SampleID)
dim(Marge_CLR) #81 samples, 1467 ASVs
#Marge_CLR:Marged 16S and 18S clr-transformed data. row is sampleID, colum is ASV. 


##Marge 16S and 18S Phylum data
Marge_Phylum <- rbind(ASV_16S.Phylum.RE, ASV_18S.Phylum.RE)
#Marge_Phylum:ASV and Taxa columns



####Input data
##Check if there are too many missing values using "goodSamplesGenes" function.
gsg = goodSamplesGenes(Marge_CLR, verbose = 1)
gsg$allOK #TRUE


##Check outlier samples using cluster analysis (Unweighted Pair Group Method with Arithmetic mean, UPGMA)
sampleTree = hclust(dist(Marge_CLR), method = "average")
#plot
par(cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off() #There is not an outlier sample.


##Read Environmental data which with match the sample order in "Marge_CLR".
names(Env.data)
TraitsData <- Env.data
rownames(TraitsData) <- TraitsData$SampleID
TraitsData %<>% dplyr::select(-SampleID)
dim(TraitsData)
names(TraitsData)

##Add Environmental data to sample cluster
sampleTree2 = hclust(dist(Marge_CLR), method = "average")
#plot
traitColors = numbers2colors(TraitsData, signed = F)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(TraitsData),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


###Make Network model
#Assuming the co-occurrence network is scale-free, set soft-Threshold using "pickSoftThreshold" function.
powers = c(c(1:10), seq(from = 11, to=20, by=1))

sft = pickSoftThreshold(Marge_CLR, powerVector = powers, verbose = 1, networkType = "signed hybrid")
#10

#plot
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed hybrid R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Soft-Threshold = 7


##Calculate adjacency using "signed hybrid" method; 
#the “signed hybrid” adjacency is exactly zero for all negative (or zero) correlations
#type = "signed hybrid"
softPower = 7
adjacency = adjacency(Marge_CLR, power = softPower, type = "signed hybrid") 
write.csv(adjacency, file = "adjacency.csv")


##Calculate TOM based on the adjacency (0 to 1) using "unsigned" method.
max(as.data.frame(adjacency))
min(as.data.frame(adjacency))
TOM = TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM = 1-TOM
#1-TOM; TOM dissimilarity


##Module identify
#Clustering based on the TOM dissimilarity
TaxaTree = hclust(as.dist(dissTOM), method = "average")
#plot
sizeGrWindow(12,9)
plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#Cut tree using 'cutreeDynamic' function. 
#The parameters in this function are deepSplit, pamStage, maxCoreScatter, minGap, and minClusterSize.
#Decide minModuleSize
minModuleSize = 30
#Identify module
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = 4, method = "hybrid", 
                            pamStage = T, pamRespectsDendro = T,
                            minClusterSize = minModuleSize)
table(dynamicMods) #The module named "0" is the unclassified ASVs. 

#Name the colors to modules
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) #The module named "grey" is the unclassified ASVs. 

#plot the result of module identifying to the tree.
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")
dev.off()


##Calculate module eigengenes
#module eigengenes; the first principal components of modules using the "moduleEigengenes" function
MEList = moduleEigengenes(Marge_CLR, colors = dynamicColors)
MEs = MEList$eigengenes
#MEs; module eigengenes


##Calculate Pearson's correlate coefficient between modules using 'pairwise.complete.obs' method.
#pairwise.complete.obs; calculate pairwise correlations ignore NA.
MEDiss = 1-cor(MEs, use = 'pairwise.complete.obs')
#Clustering modules based on the correlations.
METree = hclust(as.dist(MEDiss), method = "average")
#plot
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#Decide the similarity threshold; 0.20
MEDissThres = 0.20
abline(h=MEDissThres, col = "red")
dev.off()


##Merge similar modules
merge = mergeCloseModules(Marge_CLR, dynamicColors, cutHeight = MEDissThres, verbose = 1)
#Name the colors to modules
mergedColors = merge$colors
#Calculate module eigengenes
mergedMEs = merge$newMEs
#plot the result of merge modules
plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


##Summarize the results
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


##Add meta data of ASVs
annot = ASV_TaxID.RE
dim(annot)
names(annot)
probes = names(Marge_CLR)
probes2annot = match(probes, annot$ASV)
sum(is.na(probes2annot))
TaxaInfo0 <- data.frame(ASV = probes,
                        Taxonomy = annot$Taxonomy[probes2annot],
                        module = moduleColors)  %>%
  arrange(module, Taxonomy)
TaxaInfo1 <- left_join(TaxaInfo0, Marge_Phylum, by = "ASV")
TaxaInfo2 <- left_join(TaxaInfo1, module_names, by = "module")
TaxaInfo <- left_join(TaxaInfo2, SN_codes, by = "SNname")
head(TaxaInfo)
write.csv(TaxaInfo, file = "TaxaInfo.csv")


###Module-Traits relationship
nTaxa = ncol(Marge_CLR)
nSamples = nrow(Marge_CLR)
MEs0 = moduleEigengenes(Marge_CLR, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
write.csv(MEs, file = "MEs.csv")

rmgreyMEs <- MEs %>%
  dplyr::select(-MEgrey)
moduleTraitCor <- cor(rmgreyMEs, TraitsData, use = "pairwise.complete.obs", method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue_adj <- matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method="BH"),ncol=ncol(moduleTraitPvalue)) 
rownames(moduleTraitPvalue_adj) <- rownames(moduleTraitCor) 
colnames(moduleTraitPvalue_adj) <- colnames(moduleTraitCor) 
write.csv(moduleTraitCor, file = "moduleTraitCor.csv")
write.csv(moduleTraitPvalue_adj, file = "moduleTraitPvalue_adj.csv")

moduleTraitCor.sel <- as.data.frame(moduleTraitCor) %>%
  dplyr::select(May, CM, PC, BB, Sinking, Suspended, Free_living)# %>%
#dplyr::rename("Time" = "May")
moduleTraitPvalue_adj.sel <- as.data.frame(moduleTraitPvalue_adj) %>%
  dplyr::select(May, CM, PC, BB, Sinking, Suspended, Free_living)# %>%
#dplyr::rename("Time" = "May")

textMatrix = paste(signif(as.matrix(moduleTraitCor.sel), 2), "\n(",
                   signif(as.matrix(moduleTraitPvalue_adj.sel), 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor.sel)

labeledHeatmap(Matrix = moduleTraitCor.sel,
               xLabels = names(moduleTraitCor.sel),
               yLabels = names(rmgreyMEs),
               ySymbols = names(rmgreyMEs),
               xSymbols = names(moduleTraitCor.sel),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               yColorWidth = 0.7 * strwidth("M"),
               xColorWidth = 0.6 * strwidth("M"),
               cex.lab.x = 0.8,
               xLabelsAdj = 0.5,
               cex.lab.y = 0.8,
               main = paste("Module-trait relationships"))
dev.off()


##Cor plot (Module-Traits relationship)
rename_ModTraitsCor0 <- moduleTraitCor.sel %>%
  dplyr::select(Sinking, Suspended, Free_living, CM, PC, BB, May) %>%
  dplyr::rename("Free-living" = "Free_living") %>%
  dplyr::mutate(module = gsub(rownames(.), pattern="ME",replacement = ""))
rename_ModTraitsCor1 <- left_join(rename_ModTraitsCor0, module_names, by = "module") %>%
  transform(SNname= factor(SNname, levels = subnetwork_name)) %>%
  dplyr::arrange(SNname)
rownames(rename_ModTraitsCor1) <- rename_ModTraitsCor1$SNname
rename_ModTraitsCor <- rename_ModTraitsCor1 %>%
  dplyr::select(-module, -SNname)
#
rename_ModTraitsPval0 <- moduleTraitPvalue_adj.sel %>%
  dplyr::select(Sinking, Suspended, Free_living, CM, PC, BB, May) %>%
  dplyr::rename("Free-living" = "Free_living") %>%
  dplyr::mutate(module = gsub(rownames(.), pattern="ME",replacement = ""))
rename_ModTraitsPval1 <- left_join(rename_ModTraitsPval0, module_names, by = "module") %>%
  transform(SNname= factor(SNname, levels = subnetwork_name)) %>%
  dplyr::arrange(SNname)
rownames(rename_ModTraitsPval1) <- rename_ModTraitsPval1$SNname
rename_ModTraitsPval <- rename_ModTraitsPval1 %>%
  dplyr::select(-module, -SNname)

col <- colorRampPalette(c("deepskyblue", "white", "#FF3232"))(20)
graphics.off()
corrplot(as.matrix(rename_ModTraitsCor), method = "circle", col = col, tl.col="black", tl.srt=45, addCoef.col = "black",
         number.cex = 0.6, p.mat = as.matrix(rename_ModTraitsPval), sig.level = 0.05, insig = "blank")
graphics.off()
png("./figure/plot_ModTrait.png", height=18, width=11, units="cm", res=300)
corrplot(as.matrix(rename_ModTraitsCor), method = "circle", col = col, tl.col="black", tl.srt=45, addCoef.col = "black",
         number.cex = 0.6, p.mat = as.matrix(rename_ModTraitsPval), sig.level = 0.05, insig = "blank")
dev.off()
write.csv(rename_ModTraitsCor, file = "rename_ModTraitsCor.csv")
write.csv(rename_ModTraitsPval, file = "rename_ModTraitsPval.csv")
#####