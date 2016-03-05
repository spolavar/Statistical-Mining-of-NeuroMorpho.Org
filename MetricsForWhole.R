#loading the feature and metadata matrices
source("C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/plotutilfunctions.R")
source("C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/LoadMetricFiles.r")

#arbor type == all feature matrix
dim(metricDataWhole)
colnames(metricDataWhole)
###computing CV = SD/mean ###
LMmetricsWhole <- metricDataWhole
dim(LMmetricsWhole)

metricDataWhole[1:10,c("neuron_name","PathDistance_avg")]

chkNegVals(LMmetricsWhole[,2:ncol(LMmetricsWhole)])
t <- cntNRmvNAs(LMmetricsWhole,2)
dim(t)
NAneurons <- setdiff(LMmetricsWhole$neuron_name,t$neuron_name)
length(NAneurons)
LMmetricsWhole[which(LMmetricsWhole$neuron_name %in% NAneurons), c("neuron_name","Soma_Surface_total_sum","Depth_total_sum")]
#remove the -ve column = Helix_min and not useful column = Helix_avg
LMmetricsWhole <- LMmetricsWhole[,!(colnames(LMmetricsWhole) %in% c("Helix_min","Helix_avg"))]
#remove the neurons with NA vales in soma surface and depth columns
LMmetricsWhole <- cntNRmvNAs(LMmetricsWhole,2)
dim(LMmetricsWhole)

#if mean < 0, then remove that column from CV
chkMeanCols(LMmetricsWhole[,2:length(LMmetricsWhole)])
cvarr <- cvfun(LMmetricsWhole)
length(cvarr)

cvframe <- data.frame("metrics"=colnames(LMmetricsWhole[,2:length(LMmetricsWhole)]), "coefficient of variation" = cvarr)
cvframe <- cvframe[order(cvframe$coefficient.of.variation,decreasing=TRUE),]
dim(cvframe)
write.table(cvframe, file="cvvalues.txt",quote=FALSE, sep=",",row.names=FALSE)

whole_mu <- sapply(metricall[33:62],FUN = mean)
whole_sd <- sapply(metricall[33:62],FUN= sd)
whole_cv <- whole_sd/whole_mu
dim(metricall)
dim(metricDataWhole)

#reading feature matrix into new local variable
metricall <- metricDataWhole
colnames(metricall)
dim(metricall)
#remove rows that has all features as NA 
metricall <- cntNRmvNAs(metricall,114)

###############Data preprocessing for eliminating neurons that don't fit the metric analyses######################
#check for zero depth data or flat 2D neurons and remove them from analysis
flatwholeneuritespos <-  with(metricall, which(is.na(metricall$Depth_total_sum)))
flatwholeneurites <- metricall[flatwholeneuritespos, c("neuron_name","N_bifs_total_sum")]
print(paste("removing ",length(flatwholeneuritespos)," from", nrow(metricall), " neurons"))
metricall <- metricall[-flatwholeneuritespos,]

#metricall <- subset(metricall,!is.na(metricall$Depth_total_sum))
dim(metricall)
#remove CNG.swc from neuron_name
metricall <- clnData(metricall)
#32 metrics chosen for PCA analysis
pcawholeMetrics
length(pcawholeMetrics)
dim(metricall)
unique(usedMetaData$archive_name)
unique(usedMetaData$species)
unique(usedMetaData$region1)
unique(usedMetaData$region2)
unique(usedMetaData$region3)
unique(usedMetaData$cellclass1)
unique(usedMetaData$cellclass2)
unique(usedMetaData$cellclass3)



#merge metric data with meta data
metricall <- getPlotSubset(metricall,metrics=pcawholeMetrics)
colnames(metricall)
metricall[which(metricall$archive_name == "Rodger"),"neuron_name"]
#remove extra column that may have been inserted from merge process
metricall <- metricall[, !(colnames(metricall) %in% c("X"))]
#remove neurons with < 4 bifurcations
Nbranchoutlierpos <- with(metricall, which(metricall$N_bifs_total_sum < 4))
Nbranchoutlier <- metricall[Nbranchoutlierpos,c("neuron_name","N_bifs_total_sum","archive_name")]
str(Nbranchoutlier)
unique(Nbranchoutlier$archive_name)
print(paste("removing ",length(Nbranchoutlierpos)," from", nrow(metricall), " neurons"))
metricall <- metricall[-Nbranchoutlierpos,]
#remove Larkman archive because his neurons are dendrograms and hence aren't valid biological metrics.
larkmanpos <- with(metricall, which(metricall$archive_name=="Larkman"))
larkmanoutlier <- metricall[larkmanpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(larkmanoutlier)
print(paste("removing from larkman",length(larkmanpos)," from", nrow(metricall), " neurons"))
metricall <- metricall[-larkmanpos,]
#remove neurons with fractal==0
zerofractalpos <- with(metricall, which(metricall$Fractal_Dim_avg == 0))
zerofractaloutlier <- metricall[zerofractalpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(zerofractaloutlier)
print(paste("removing zero fractals",length(zerofractalpos)," from", nrow(metricall), " neurons"))
metricall <- metricall[-zerofractalpos,]
colnames(metricall)
#check for zero soma surface
zerosomasurface <- with(metricall, which(metricall$Soma_Surface_total_sum == 0))
zerosomaoutlier <-  metricall[zerosomasurface,c("neuron_name","archive_name")]
nrow(zerosomaoutlier)
unique(zerosomaoutlier$archive_name)
print(paste("removing zero soma surface",length(zerosomasurface)," from", nrow(metricall), " neurons"))
metricall <- metricall[-zerosomasurface,]

hist(metricall$Helix_max, breaks=30, xlab="Helix")
hist(metricall$Fractal_Dim_avg, breaks=30, xlab="fractal")

min(metricall$Helix_max)
helix_zeropos <- with(metricall, which(metricall$Helix_max == 0.00))
helixzero <- metricall[helix_zeropos,c("neuron_name","Fractal_Dim_avg","Contraction_avg","Branch_pathlength_max","archive_name")]
dim(helixzero)
dim(metricall)
unique(metricall$cellclass3)
#### log transform features that are NOT normal ####
#the normality is verified by the skewness of the distribution. 
plotDf <- plotTesting(metricall)
colnames(plotDf)
hist(plotDf[,62], breaks=30, xlab="Helix")
hist(plotDf[,64], breaks=30, xlab="fractal")
unique(plotDf$cellclass3)
colnames(plotDf)
dim(plotDf)
#standardize metric dataframe before computing PCA
for(i in 33:64){
  #get non-NA values for normalization
  colv <- plotDf[,i][!is.na(plotDf[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  plotDf[,i][!is.na(plotDf[,i])] <- colv/msd
  print(paste(names(plotDf)[i],msd))
}

#### check for metricall metadata content and cell count####
dim(metricall)
cf <- metricall[which(metricall$cellclass3 == "Climbing fiber"),]
dim(cf)
colnames(cf)
cf$Helix_max
#Assign Group labels to the feature matrix on 'mxdgroupnum' based on archive,celltype, brain region and species metadata
###########add cell types groups##############
aggregate(plotDf$expercond, by = list(plotDf$expercond), length)
#celltype hierarchy
primary <- list(expercond="Control")
x <- getSubset2(plotDf, primary, typeMetaData=usedMetaData)
level1 <- "cellclass1"
level2 <- "cellclass2"
level3 <- "cellclass3"
#make groups following hierarchical ordering of metadata
ctMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData, level1, level2,level3,minSize=100)
printGrps(ctMetricHierarchy)
#group neuron_name (primariy id) according the hierarchy lists
Grplst <- addgrps2MetricDF(plotDf,ctMetricHierarchy)
for(i in 1:length(Grplst)){
  print(Grplst[[i]]$GroupSpec)
}
#assigning a groupnum and groupttl to each neuron according to the hierarchy list
plotDf_ct <- metricDf_plotting(plotDf,Grplst)
dim(plotDf_ct)

###########add brain region groups############
primary <- list(expercond = "Control")
colnames(usedMetaData)
level1 <- "region1"
level2 <- "lobes"
level3 <- "region2"
level4 <- "layer"
unique(usedMetaData$lobes)
#a difference of 22 neurons was found between regions and lobes grouping on 'Hippocampus'. To fix this, find out the neuron names and add their missing lobes values.
ss <- getSubset2(plotDf, list(region1 = "Hippocampus"), typeMetaData=usedMetaData)
dim(ss)
ss1 <- getSubset2(plotDf, list(lobes = "Hippocampus"), typeMetaData=usedMetaData)
dim(ss1)
t <- setdiff(ss$neuron_name, ss1$neuron_name)
plotDf[which(plotDf$neuron_name %in% t),]['lobes'] <- "Not reported"
tt <- subset(ss, ss$neuron_name %in% t)
unique(subset(plotDf,plotDf$neuron_name %in% t)[,c("region1","lobes","region2","region3","layer")])
aggregate(plotDf[which(plotDf$lobes == "Hippocampus"),], by=list(plotDf[which(plotDf$lobes == "Hippocampus"),]$region2), FUN= length)
aggregate(plotDf, by=list(plotDf$region1), FUN= length)
#### removed the cutoff value to include all data in the cluster analysis. Due to this the groups will have "not reported" as part of the group label. For the crosshair plots the "not reported" labels are eliminated.####
brMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData,level2,level3,level4,minSize=300)

#brMetricHierarchy
printGrps(brMetricHierarchy)
#str(brMetricHierarchy)

Grplst <- addgrps2MetricDF(plotDf,brMetricHierarchy)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}
plotDf_br <- metricDf_plotting(plotDf,Grplst)
colnames(plotDf_br)
unique(plotDf_br$groupttl)
dim(plotDf_br)

###########add species groups#################
primary <- list(expercond = "Control")
level1 <- "order"
level2 <- "species"
level3 <- "strain"
spstrainMetricall <- makeHierarchyGroups(primary, plotDf,usedMetaData,level1,level2,minSize=55)#,level3,
printGrps(spstrainMetricall)
Grplst <- addgrps2MetricDF(plotDf,spstrainMetricall)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}
plotDf_sp <- metricDf_plotting(plotDf,Grplst)
dim(plotDf_sp)
colnames(plotDf_sp)
dim(plotDf_sp)

########prepare data matrix for PCA without NAs in groups and metrics##########

#add rownames for three dimensions
rownames(plotDf_ct) <- plotDf_ct$neuron_name
rownames(plotDf_br) <- plotDf_br$neuron_name
rownames(plotDf_sp) <- plotDf_sp$neuron_name
dim(plotDf_ct)
dim(plotDf_br)
dim(plotDf_sp)
colnames(plotDf_ct)

dim(subset(plotDf_ct, !is.na(plotDf_ct$groupnum)))

plotDf_ct <- subset(plotDf_ct, !is.na(plotDf_ct$groupnum))

dim(plotDf_ct)
colnames(plotDf_ct[,33:64])
#handles missing values by removing rows that have NAs - remove soma_surface from PCA analysis as not all data have soma_surface
pcametrics <- na.omit(plotDf_ct[,33:64])
dim(pcametrics)
pcametrics <- subset(plotDf_ct, plotDf_ct$neuron_name%in%rownames(pcametrics))
colnames(pcametrics)
names(pcametrics)[65] <- "ct_groupnum"
names(pcametrics)[66] <- "ct_groupttl"

#add brain region groups
nonpcaDf <- subset(plotDf_br, plotDf_br$neuron_name%in%rownames(pcametrics))
dim(nonpcaDf)
colnames(nonpcaDf)
pcametrics[,"br_groupnum"] <- nonpcaDf$groupnum
pcametrics[,"br_groupttl"] <- nonpcaDf$groupttl
#add species groups
#plotDf_sp <- subset(plotDf_sp,!is.na(plotDf_sp$groupnum))
dim(plotDf_sp)
colnames(plotDf_sp)
nonpcaDf <- subset(plotDf_sp, plotDf_sp$neuron_name%in%rownames(pcametrics))
dim(nonpcaDf)
colnames(nonpcaDf)
pcametrics[,"sp_groupnum"] <- nonpcaDf$groupnum
pcametrics[,"sp_groupttl"] <- nonpcaDf$groupttl

dim(pcametrics)
colnames(pcametrics)
#remove NA rows from br_groupnum and sp_groupnum
pcametrics <- subset(pcametrics, !is.na(pcametrics$br_groupnum))
pcametrics <- subset(pcametrics, !is.na(pcametrics$sp_groupnum))

#build data matrix with PCs, group tags

pcametrics <- editBrainRegionLabels(pcametrics,"br_groupttl")
unique(pcametrics$br_groupttl)
pcametrics <- editCelltypeLabels(pcametrics,"ct_groupttl")
unique(pcametrics$ct_groupttl)
pcametrics <- editspeciesLabels(pcametrics,"sp_groupttl")
unique(pcametrics$sp_groupttl)

###########add mixed groups####################
#make mixedgroups
colnames(pcametrics)
dim(pcametrics)
fourdims <- aggregate(pcametrics$neuron_name, as.list(pcametrics[,c("archive_name","ct_groupttl","br_groupttl","sp_groupttl")]), FUN = length)
fourdims
nrow(fourdims)
dim(fourdims)
min(fourdims$x)
max(fourdims$x)
#groups from 4D space combining archives/species/celltypes/brain regions
mixedgrpsDf <- subset(fourdims,fourdims$x>=1)#55
dim(mixedgrpsDf)
mixedgrpsDf
#prepare data matrix for computing pca and for cluster analysis
datamatrix <- addmixedgrpcol(mixedgrpsDf, pcametrics)
dim(datamatrix)
mxdgrpctlg <- aggregate(datamatrix$neuron_name, as.list(datamatrix[,c("archive_name","ct_groupttl","br_groupttl","mxdgroupnum")]), FUN = length)
write.table(mxdgrpctlg, file="mxdgrpTable.txt",quote=TRUE, sep=",")
colnames(datamatrix)
unique(datamatrix$mxdgroupnum)
#check for correlation before running PCA
prePCAcorr <- mostcorrelated(datamatrix[,33:64])
str(prePCAcorr)
colnames(prePCAcorr)
#eliminate one of the highly correlated pairs
removeFeatures <- eliminateCorrFeat(prePCAcorr)
#remove the selected correlated features from data matrix
datamatrix <- datamatrix[ , -which(names(datamatrix) %in% removeFeatures$feature)]
colnames(datamatrix)
###########pca for dimensionality reduction and determine #clusters for whole arbor######
colnames(datamatrix[,33:60])
colnames(pcametrics[,33:60])
dim(datamatrix)
pcamat <- prcomp(x=datamatrix[,33:60],center=T,scale=T)
names(str(pcamat))
pcasumm <- summary(pcamat)
unique(datamatrix$archive_name)
unique(pcametrics$archive_name)
nrow(pcametrics[which(pcametrics$archive_name == "Gulyas"),])


#make scree plot to choose most contributing PCs
plot(pcamat$sdev^2/sum(pcamat$sdev^2)*100, type = "b", ylim = c(0,100), xlim=c(0,32), main = "Scree plot for whole arbor",xlab="PCA Eigenvalues",ylab="% of variance")
#plot cumulative variance
lines(pcasumm$importance[3,]*100,main = "whole arbor metrics",type = "b", col = 2, xlab="cumulated PCA eigen values",ylab="% of variance")
legend(10,50,c("cumulative variance","variance"),lty=c(1,1),col=c("red","black"))

#six components were covering ~75% of variance
#choose eigenvalues >= 1
pcamat$sdev ^ 2
#writing the highest loadings on each component into a file
colnames(pcametrics)
loadall <- abs(pcamat$rotation)
dim(loadall)
PCloadingsarr <- matrix(nrow=28,ncol=3)
dim(PCloadingsarr)
for(c in 1:28){
  rowid <- which(loadall == max(loadall[,c]),arr.ind=TRUE)[1]
  colid <- which(loadall == max(loadall[,c]),arr.ind=TRUE)[2]
  metricname <- names(loadall[,colid])[rowid]
  PCname <- names(loadall[rowid,])[colid]
  PCloadingsarr[c,1] <- PCname
  PCloadingsarr[c,2] <- metricname
  PCloadingsarr[c,3] <- pcamat$rotation[rowid, colid]
}

PCloadingsarr
write.table(PCloadingsarr, file="PCloadingsarr.txt",quote=FALSE, sep=",")
#loadings of PC1 and PC2
pc1loadings <-sort(pcamat$rotation[,1],decreasing=F)
dotchart(pc1loadings,main="PC1 loadings",cex=0.7, xlab="variable loadings",col=2)
pc2loadings <-sort(pcamat$rotation[,2],decreasing=F)
dotchart(pc2loadings,main="PC2 loadings",cex=0.7, xlab="variable loadings",col=2)
pc3loadings <-sort(pcamat$rotation[,3],decreasing=F)
dotchart(pc3loadings,main="PC3 loadings",cex=0.7, xlab="variable loadings",col=2)
pc4loadings <-sort(pcamat$rotation[,4],decreasing=F)
dotchart(pc4loadings,main="PC4 loadings",cex=0.7, xlab="variable loadings",col=2)
pc5loadings <-sort(pcamat$rotation[,5],decreasing=F)
dotchart(pc5loadings,main="PC5 loadings",cex=0.7, xlab="variable loadings",col=2)
pc6loadings <-sort(pcamat$rotation[,6],decreasing=F)
dotchart(pc6loadings,main="PC6 loadings",cex=0.7, xlab="variable loadings",col=2)
pc7loadings <-sort(pcamat$rotation[,7],decreasing=F)
dotchart(pc7loadings,main="PC7 loadings",cex=0.7, xlab="variable loadings",col=2)
pc24loadings <-sort(pcamat$rotation[,24],decreasing=F)
dotchart(pc24loadings,main="PC24 loadings",cex=0.7, xlab="variable loadings",col=2)

pc8loadings <-sort(pcamat$rotation[,8],decreasing=F)
#pc8loadings <-sort(abs(pcamat$rotation[,8]),decreasing=F)
dotchart(pc8loadings,main="PC8 loadings",cex=0.7, xlab="variable loadings",col=2)

#create a data frame that has Principal components and metadata groups
dim(datamatrix)
colnames(datamatrix)
dim(pcamat$x)
wholeReduced <- data.frame(datamatrix[,1:32], pcamat$x, datamatrix[,61:67])
dim(wholeReduced)
unique(wholeReduced$ct_groupnum)
length(unique(wholeReduced$ct_groupnum))
unique(wholeReduced$br_groupnum)
length(unique(wholeReduced$br_groupnum))
unique(wholeReduced$sp_groupnum)
length(unique(wholeReduced$sp_groupnum))

############cluster the PC space using EM#####################
library(cluster)
library(fpc)
library("leaps")

colnames(wholeReduced)
dim(wholeReduced)
#cluster PCs of feature matrix
tmpEM <- Mclust(wholeReduced[,grep("[PC]", names(wholeReduced), value=TRUE)],G=1:15)

#emResult <- Mclust(clstrMatrix,G=1:15)
summary(tmpEM)
plot(tmpEM, xlab = "Number of clusters", data=wholeReduced[,grep("[PC]", names(wholeReduced), value=TRUE)], what='BIC')
plot(tmpEM, what="classification",dimens=1:3)
plot(tmpEM, what="classification",dimens=7:9)
plot(tmpEM, what="classification",dimens=10:12)
unique(tmpEM$classification)
#merge classification with PC feature matrix
wholeReduced$classification = tmpEM$classification
datamatrix$classification = wholeReduced$classification
#clustering on reduced dimensions accounting for 95% of variance in the data
dim(wholeReducedtmpEM_95 <- Mclust(wholeReduced[,33:50],G=1:15)
wholeReduced$classification = tmpEM_95$classification
plot(regsubsets(y=tmpEM_95$classification,x=wholeReduced[,33:50]))
write.table(wholeReduced, file="PCMetricFeatureDF.txt",quote=FALSE, sep=",",row.names=FALSE)
dim(tmpEM_95$classification)
summary(tmpEM_95)
unique(tmpEM_95$classification)
plot(tmpEM_95, data=wholeReduced[,33:50], what='BIC')
plot(tmpEM_95, what="classification",dimens=1:3)
plot(tmpEM, what="classification",dimens=7:9)
plot(tmpEM, what="classification",dimens=10:12)
#merge classification from 95% with PC feature matrix
wholeReduced$classification = tmpEM_95$classification
datamatrix$classification = wholeReduced$classification

#clustering on reduced dimensions accounting for 95% of variance in the data
tmpEM_7 <- Mclust(wholeReduced[,33:37],G=1:15)
summary(tmpEM_7)
unique(tmpEM_5$classification)
plot(tmpEM_5, data=wholeReduced[,33:50], what='BIC')
plot(tmpEM, what="classification",dimens=1:3)
plot(tmpEM, what="classification",dimens=7:9)
plot(tmpEM, what="classification",dimens=10:12)

#############Find significant group-cluster associations from binomial test##############

alpha <- 0.05
#takes input parameters: single metadata feature and the original data frame with PCs.
grpfeature <- 'mxdgroupnum'
clusfeature <- 'classification' #classification95'

#add cluster classification to featPCDf
groupmatrix <- aggregate(wholeReduced$neuron_name, by=list(groupnum=wholeReduced[,grpfeature]), FUN=length)
print(groupmatrix)
clustermatrix <- aggregate(wholeReduced$neuron_name, by=list(clusternum=wholeReduced[,clusfeature]), FUN=length)
print(clustermatrix)

grpclustermatrix <- aggregate(wholeReduced$neuron_name, by=list(groupnum=wholeReduced[,grpfeature],
                                                              clusternum=wholeReduced[,clusfeature]), FUN=length)
groupvals <- groupmatrix$groupnum #unique(featPCDf[,grpfeature])
clustervals <- clustermatrix$clusternum #unique(featPCDf$classification)
  
row_names <- groupvals
col_names <- clustervals
#determine the #group values in a single metadata feature
g_num <- length(groupvals)#length(unique(featPCDf[,grpfeature]))
#run EM based clustering on featPCDf only on PC columns
#determine the #clusters from that
c_num <- length(clustervals)#nrow(clustermatrix)
dataN <- nrow(wholeReduced)
#grpfeature aggregate on unique group values
print(grpclustermatrix)
#compute observed and expected matrices
tmplst <- getObservedAndExpectedMatrices(groupmatrix,clustermatrix,grpclustermatrix,dataN)
observedMat <- tmplst[[1]]
expectedMat <- tmplst[[2]]
expectedMat_integer <- tmplst[[3]]
print(dim(observedMat))
print(dim(expectedMat_integer))
print(paste(g_num,c_num))
###run chi sqaure test on #groups and #columns
chisqtest <- chisq.test(observedMat)
stdresMatrix_chi <- chisqtest$stdres        #standardized residuals 
observedMat_chi <- chisqtest$observed
expectedMat_chi <- floor(chisqtest$expected+0.5)

probMatrix_chi <- pnorm(-1 * abs(stdresMatrix_chi)) * 2  #compute p-vals by reversing to the -ve side and multiplying by 2
probMatrix_chi[which(probMatrix_chi==0,arr.ind=T)] <- pnorm(-37.5)*2 #minimum value for which p-val can be calculated

probMatrix_chi <- withinTableCorrectionBf(probMatrix_chi)
problog10Matrix_chi <- converttolog10(probMatrix_chi,observedMat_chi,expectedMat_chi/dataN,dataN)   #convert pvals to log10 values

###run binom test on #groups and #columns
#probMatrix <- runBinomtest(observedMat, expectedMat, dataN, g_num, c_num)
###convert pvals to log10(pvals)
#probMatrixlog10 <- converttolog10(probMatrix,observedMat,expectedMat,dataN)

###check for overshooting groups and undershooting on chisq table###
sigcond_over <- which(problog10Matrix_chi>0, arr.ind=TRUE)
#get groups and clusters indicies for all significant differences both over and under representation
sigcond_all <- which(!is.na(problog10Matrix_chi), arr.ind=TRUE)
nrow(sigcond_all)
str(sigcond_over)
sigcond_over[,'row']
sigcond_over[,'col']
nrow(sigcond_over)
sigMatrix <- getSignificantGroupsNClusters(sigcond_over,wholeReduced,clusfeature)
sigMatrix_all <- getSignificantGroupsNClusters(sigcond_all,wholeReduced, clusfeature)


###analyzing by groups: the ones that split in two groups###
library("car")
library("tourr")
sigMatrix <- calculateGCProportions(sigMatrix,observedMat_chi)

groupoccurences <- aggregate(sigMatrix$row, by=list(sigMatrix$row), FUN=length)
multiplegrps <- groupoccurences[with(groupoccurences, which(groupoccurences$x==2)),"Group.1"]
#select groups that are divided into 2 clusters

sigMatrix_s <- groupoccurences[groupoccurences$x==2,]
groupoccurences[groupoccurences$x==1,]
write.table(sigMatrix, file="OverRepresentedClusters.txt",quote=TRUE, sep=",",row.names=FALSE)
###print matrices to file####
#write the chisq test pval matrix to file
write.table(round(problog10Matrix_chi,2), file="chisq_log10pval.txt",quote=FALSE, sep=",",row.names=TRUE)
#write data matrices into output files for analysis
write.table(sigMatrix_s, file="sigMatrix_s.txt",quote=TRUE, sep=",",row.names=FALSE)
#The observed matrix
observedMatrix
write.table(observedMat_chi, file="ObservedMatrix.txt",quote=FALSE, sep=",",row.names=TRUE)
#the expected matrix converted to integer for ease of comparison
expectedMatrix
write.table(expectedMat_chi, file="expectedMatrix.txt",quote=FALSE, sep=",",row.names=TRUE)
#the delta matrix
write.table((observedMat_chi - expectedMat_chi), file="deltaMatrix.txt",quote=FALSE, sep=",",row.names=TRUE)
#the ratio matrix
dataN
ratioMatrix <- round(observedMat_chi/expectedMat_chi,2)
ratioMatrix
write.table(ratioMatrix, file="ratioMatrix.txt",quote=FALSE, sep=",",row.names=TRUE)
#probability matrix
problog10Matrix <- round(problog10Matrix,2)
write.table(problog10Matrix, file="problog10Matrix.txt",quote=FALSE, sep=",",row.names=FALSE)

#plot all 7 clusters for whole arbor data
dim(wholeReduced)
colarr <- colorcode(groupnum=wholeReduced$classification, length(unique(wholeReduced$classification)))

plot(wholeReduced[,33],wholeReduced[,34],col="white",cex=.55,pch=16,xlab=names(wholeReduced)[33],
     main=paste(dataN,"whole arbor data classified into 7 models"),ylab=names(wholeReduced)[34],cex.main=0.8)
text(wholeReduced[,33],wholeReduced[,34],labels=as.character(wholeReduced$classification),col=mypalette[colarr],cex=.7)
legend("topleft",text.width = strwidth("1"), legend = unique(wholeReduced$classification), pch=19, cex=0.7, col=unique(mypalette[colarr]))
clusnums <- unique(wholeReduced$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

colarr <- colorcode(groupnum=wholeReduced$classification95, length(unique(wholeReduced$classification95)))

plot(tmpEM,what="classification",dimens=1:2,cex =0.8,pch=1)
summary(tmpEM)
tmpEM$classification
clusnums <- unique(wholeReduced$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

#plotting clusters with 
colarr <- colorcode(sigMatrix_s$Group.1,length(unique(sigMatrix_s$Group.1)))
grpscluster <- subset(sigdf, sigdf$mxdgroupnum%in%unique(sigMatrix_s$Group.1))
dim(grpscluster)
unique(grpscluster$mxdgroupnum)
#plot PCs in cluster plotting
colarr <- colorcode(groupnum=grpscluster$mxdgroupnum, length(unique(grpscluster$mxdgroupnum)))
plot(grpscluster[,33],grpscluster[,34],col="white",cex=.55,pch=16,xlab=names(grpscluster)[33],
     main="11 out of 37 groups that are divided between 2 clusters",ylab=names(grpscluster)[34],cex.main=0.8)
text(grpscluster[,33],grpscluster[,34],labels=as.character(grpscluster$classification),col=mypalette[colarr],cex=.7)
legend("topleft",text.width = strwidth("10"), legend = unique(grpscluster$mxdgroupnum), pch=19, cex=0.7, col=unique(mypalette[colarr]))
clusnums <- unique(grpscluster$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

animate(grpscluster[,33:37], grand_tour(), display_xy(col=mypalette[colarr]))
#group labels
legend("bottomleft", 
       legend = unique(grpscluster$mxdgroupnum), 
       pch=19, cex=0.7, col=unique(mypalette[colarr]))

metagrpsList <- c("layer","cellclass3","age_class","gender","protocol", "slice_thickness","slicing_direction","stain", "magnification")
length(metagrpsList)
#unique(wholeReduced$protocol)

#add pmid column to wholeReduced 
pmidfile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/pmid.txt"
pmidcol <- read.csv(pmidfile,header = TRUE, sep="\t",na.strings=c("NA",""),quote="\"",stringsAsFactors=FALSE)
dim(pmidcol)
colnames(pmidcol)
wholeReduced <- merge(x=wholeReduced,y=pmidcol,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(wholeReduced)
colnames(wholeReduced)
#remove 15 cells from alexa fluor 594 data
svobodaS1cells <- sigdf[with(sigdf, which(sigdf$mxdgroupnum == 35 & sigdf$stain=='Alexa Fluor 594')), c('neuron_name')]
wholeReduced <- subset(wholeReduced, !(wholeReduced$neuron_name %in% alexa594data))
dim(wholeReduced)
grpsw2clusters <- compare2clusters(sigMatrix_s,problog10Matrix_chi,metagrpsList)
#extract the returned variables from groups that are divided into 2 clusters 
str(grpsw2clusters)
#extract the information for group#12 which is a 3rd element in the return list 
#dataframe - each element has 4 variables
ss <- grpsw2clusters[12]$subgrpDF
#pvalue of metatype- 1st variable
metatyp <- names(grpsw2clusters[9]$pvalue)
featurecol <- ss[,metatyp]
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
grp12MetricData <- subset(datamatrix, rownames(datamatrix) %in% rownames(ss))
#featurecol <- grp12MetricData[,"cellclass3"]
dim(grp12MetricData)
colnames(grp12MetricData)
grpN <- unique(grp12MetricData$mxdgroupnum)
grpttl <- paste(grp12MetricData$ct_groupttl[1], grp12MetricData$br_groupttl[1], grp12MetricData$sp_groupttl[1])
Xname <- 'ln(Length_total_sum)' #'Partition_asymmetry_avg' #'Bif_tilt_local_max' 
Yname <- 'ln(Branch_pathlength_avg)'#'Bif_ampl_remote_avg' #

#check correlation
test <- grp12MetricData[,c("Bif_tilt_local_max","Bif_ampl_remote_avg")]
test <- grp12MetricData[,c('Partition_asymmetry_avg','ln(Branch_pathlength_avg)')]
cor(test)
plot(prcomp(test)$x, col )

#rerun EM clustering on the subset of original metrics
#cluster PCs of feature matrix
reEM <- Mclust(grp12MetricData[,33:64],G=1:15)
#emResult <- Mclust(clstrMatrix,G=1:15)
summary(reEM)
plot(reEM, data=grp12MetricData[,33:64], what='BIC')
plot(reEM, what="classification",dimens=c(14,17,20,23),cex = 0.55, cex.sub=0.7,col=c("black", "darkorange"))
plot(reEM, what="classification",dimens=c(20,23),cex = 0.55, cex.sub=0.7,col=c("black", "darkorange"))
grp12MetricData$classification <- reEM$classification

#do EM on only the LY stain data
grp12LYOnly <- subset(grp12MetricData, grp12MetricData$stain == "Lucifer Yellow")
grp12GolgiOnly <- subset(grp12MetricData, grp12MetricData$stain == "Golgi")
dim(grp12LYOnly)
dim(grp12GolgiOnly)
reEM_LY <-  Mclust(grp12LYOnly[,33:64],G=1:15)
summary(reEM_LY)

reEM_Golgi <- Mclust(grp12GolgiOnly[,33:64],G=1:32)
summary(reEM_Golgi)

#plotting the original morphometrics
plot(grp12MetricData[,Xname],grp12MetricData[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(grp12MetricData[,Xname],grp12MetricData[,Yname],labels=as.character(grp12MetricData$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#overlay ellipsoids on the clusters
clusnums <- unique(grp12MetricData$classification)
parms <- reEM$parameters

for(i in 1:length(clusnums)){
  XYpcs <- grp12MetricData[with(grp12MetricData,which(grp12MetricData$classification == clusnums[i])),c(Xname,Yname)]
  XYmeans <- colMeans(XYpcs)
  
  #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
  #       radius=1, center.pch=as.character(clusnums[i]),col=1)
  ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

#####end of group 12#####
#extract the information for group#19 which is 4th element in the return list
#dataframe - each element has 4 variables
ss <- grpsw2clusters[16]$subgrpDF
#pvalue of metatype- 1st variable
metatyp <- names(grpsw2clusters[13]$pvalue)
featurecol <- ss[,metatyp]
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
grp19MetricData <- subset(datamatrix, rownames(datamatrix) %in% rownames(ss))
dim(grp19MetricData)
colnames(grp19MetricData)
grpN <- unique(grp19MetricData$mxdgroupnum)
grpttl <- paste(grp19MetricData$ct_groupttl[1], grp19MetricData$br_groupttl[1], grp19MetricData$sp_groupttl[1])
Xname <- 'ln(Bif_tilt_remote_max)' #'Partition_asymmetry_avg' #'Bif_tilt_local_max' #dominant feature on PC1
Yname <- 'ln(Branch_pathlength_avg)'#'Bif_ampl_remote_avg' #dominant feature on PC2

#check correlation
test <- grp19MetricData[,c("Bif_tilt_local_max","Bif_ampl_remote_avg")]
test <- grp19MetricData[,c('Partition_asymmetry_avg','ln(Branch_pathlength_avg)')]
cor(test)

#get apical & basal dendrites Only data
dim(Dendritedatamatrix)
dim(datamatrix)
grp19DendriteData <- subset(Dendritedatamatrix, Dendritedatamatrix$neuron_name %in% grp19MetricData$neuron_name)
dim(grp19DendriteData)
colnames(grp19DendriteData)
dim(grp19MetricData)
colnames(grp19MetricData)

#get axonOnlydata
dim(Axondatamatrix)
dim(datamatrix)
grp19AxonData <- subset(Axondatamatrix, Axondatamatrix$neuron_name %in% grp19MetricData$neuron_name)
dim(grp19AxonData)

#rerun EM clustering on the subset of original metrics
#cluster PCs of feature matrix
reEM <- Mclust(grp19MetricData[,33:64],G=1:15)
grp19MetricData$classification <- reEM$classification
reEM_dend <- Mclust(grp19DendriteData[,33:63],G=1:15)
grp19DendriteData[,"classification"] <- reEM_dend$classification
reEM_Ax <- Mclust(grp19AxonData[,33:63],G=1:15)
grp19AxonData[,"classification"] <- reEM_Ax$classification
summary(reEM_Ax)
#emResult <- Mclust(clstrMatrix,G=1:15)
summary(reEM)

#plotting whole arbor metrics
plot(grp19MetricData[,Xname],grp19MetricData[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(grp19MetricData[,Xname],grp19MetricData[,Yname],labels=as.character(grp19MetricData$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#plotting dendrite only morphometrics
plot(grp19DendriteData[,Xname],grp19DendriteData[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(grp19DendriteData[,Xname],grp19DendriteData[,Yname],labels=as.character(grp19DendriteData$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#plotting axon only data
plot(grp19AxonData[,Xname],grp19AxonData[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(grp19AxonData[,Xname],grp19AxonData[,Yname],labels=as.character(grp19AxonData$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#overlay ellipsoids on the clusters
clusnums <- unique(grp19DendriteData$classification)
parms <- reEM$parameters

for(i in 1:length(clusnums)){
  XYpcs <- grp12MetricData[with(grp12MetricData,which(grp12MetricData$classification == clusnums[i])),c(Xname,Yname)]
  XYmeans <- colMeans(XYpcs)
  
  #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
  #       radius=1, center.pch=as.character(clusnums[i]),col=1)
  ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

#distmat <- res[2]

###testing subgroups###
femaleInsula <- wholeReduced[with(wholeReduced, which(wholeReduced$mxdgroupnum==6 & wholeReduced$gender=="F")),c("PC1","PC2","ct_groupttl","br_groupttl","sp_groupttl")]
maleInsula <- wholeReduced[with(wholeReduced, which(wholeReduced$mxdgroupnum==6 & wholeReduced$gender=="M")),c("PC1","PC2","ct_groupttl","br_groupttl","sp_groupttl")]
dim(femaleInsula)
dim(maleInsula)
fxmin <- min(min(femaleInsula$PC1), min(maleInsula$PC2))
mxmin <- min(min(maleInsula$PC1),min(maleInsula$PC2))
fxmax <- max(max(femaleInsula$PC1),max(femaleInsula$PC2))
mxmax <- max(max(maleInsula$PC1),max(maleInsula$PC2))

hist(maleInsula$PC1, xlim=c(mxmin, mxmax), xlab="male PC1 & PC2", breaks=30, col="red",cex.main=0.8, main=paste(maleInsula$ct_groupttl[1],maleInsula$br_groupttl[1],maleInsula$sp_groupttl[1]))
hist(maleInsula$PC2, add = T, breaks=30, col=rgb(0, 1, 0, 0.5))
legend('topleft',c('PC1','PC2'),
      fill = c('red',rgb(0, 1, 0, 0.5)), bty = 'n',
       border = NA)

hist(femaleInsula$PC1, xlim=c(fxmin, fxmax), xlab="female PC1 & PC2", breaks=30, col="red",cex.main=0.8, main=paste(femaleInsula$ct_groupttl[1],femaleInsula$br_groupttl[1],femaleInsula$sp_groupttl[1]))
hist(femaleInsula$PC2, add = T, breaks=30, col=rgb(0, 1, 0, 0.5))
legend('topleft',c('PC1','PC2'),
       fill = c('red',rgb(0, 1, 0, 0.5)), bty = 'n',
       border = NA)

pc1xmin <- min(min(maleInsula$PC1), min(femaleInsula$PC1))
pc1xmax <- max(max(maleInsula$PC1), max(femaleInsula$PC1))
hist(maleInsula$PC1, xlim=c(pc1xmin, pc1xmax), xlab="PC1", breaks=30, col="red",cex.main=0.8, main=paste(maleInsula$ct_groupttl[1],maleInsula$br_groupttl[1],maleInsula$sp_groupttl[1]))
hist(femaleInsula$PC1, add = T, breaks=30, col=rgb(0, 1, 0, 0.5))
legend('topleft',c('male','female'),
       fill = c('red',rgb(0, 1, 0, 0.5)), bty = 'n',
       border = NA)

pc2xmin <- min(min(maleInsula$PC2), min(femaleInsula$PC2))
pc2xmax <- max(max(maleInsula$PC2), max(femaleInsula$PC2))
hist(maleInsula$PC2, xlim=c(pc2xmin, pc2xmax), xlab="PC2", breaks=30, col="red",cex.main=0.8, main=paste(maleInsula$ct_groupttl[1],maleInsula$br_groupttl[1],maleInsula$sp_groupttl[1]))
hist(femaleInsula$PC2, add = T, breaks=30, col=rgb(0, 1, 0, 0.5))
legend('topleft',c('male','female'),
       fill = c('red',rgb(0, 1, 0, 0.5)), bty = 'n',
       border = NA)

#FDR correction on given Matrix. A matrix 

#bonferroni correction for within table, the input parameter is a matrix
withinTableCorrectionBf <- function(pvaltable){
  dfcol = nrow(pvaltable) - 1
  dfrow = ncol(pvaltable) - 1
  df = dfcol * dfrow
  return(pvaltable*df)
}

#run binom test to compute probability matrix  
runBinomtest <- function(observedMat, expectedMat, dataN,rowsz, colsz ){
  #size of matrix
  bonferroni_correction <- rowsz*colsz
  #create probability matrix from binomial test
  probMatrix <- matrix(data = -9, nrow = rowsz, ncol = colsz, byrow=FALSE)
  #dimnames(probMatrix) <- list(row_names, col_names)
  #print(paste("probMatrix:",dim(probMatrix)))
  #row index
  for(i in 1:rowsz){
    #colindex refferring to #column
    for(j in 1:colsz){
      #print(paste(i,j))
      #print(paste(observedMat[i,j],expectedMat[i,j]))
      
      #populate probMatrix byrow with exact binomial test pval
      probMatrix[i,j] <- binom.test(observedMat[i,j],
                                    dataN,
                                    expectedMat[i,j],
                                    alternative="two.sided")$p.value*bonferroni_correction
      
    }
  }
  return (probMatrix)
}

#convert pvals to log10(pvals)
converttolog10 <- function(probMatrix,observed,expected,size){
  row_names <- row.names(probMatrix)
  col_names <- colnames(probMatrix)
  gind <- length(row_names)
  cind <- length(col_names)
  N <- size
  probMatrixlog10 <- log10(probMatrix)
  #print(probMatrixlog10)
  #the new cutoff is abs(probMatrixlog10) >= abs(log10(alpha)) == 1.3
  #set all values > log(1) == 0 to NA 
  probMatrixlog10[which(probMatrixlog10>0,arr.ind=T)] <- NA
  #set all values whose absolute value < 1.3 as NA (non-significant after bonferroni)
  probMatrixlog10[which(abs(probMatrixlog10)<abs(log10(alpha)),arr.ind=T)] <- NA
  #set -ve values for undershoot and +ve for overshoot as convention
  for(i in 1:gind){
    #colindex refferring to #column
    for(j in 1:cind){
      
      #if observed values < expected values, then add log10(pvalue)
      if(observed[i,j] < expected[i,j]*N)
        #use -ve convention for undershoot
        probMatrixlog10[i,j] <- probMatrixlog10[i,j]
      else
        #use +ve convention for overshoot
        probMatrixlog10[i,j] <- -1 * probMatrixlog10[i,j]
    }
  }
  dimnames(probMatrixlog10) <- list(row_names,col_names)
  return (probMatrixlog10)
}

getObservedAndExpectedMatrices <- function(groupdf,clusterdf, grpclusterdf,clusN){
  groupvals <- groupdf$groupnum
  clustervals <- clusterdf$clusternum
  g_num <- length(groupvals)
  c_num <- length(clustervals)
  dataN <- clusN
  row_names <- groupvals
  col_names <- clustervals
  #print(paste(g_num,c_num,dataN))
  #initialize the assignment matrix
  tmpMatrix <- matrix(-9, g_num, c_num)
  tmpMatrix_1 <- matrix(-9, g_num, c_num)
  #print(dim(tmpMatrix_1))
  tmpMatrix_2 <- matrix(-9, g_num, c_num)
  #create Observed and expected matrices
  for(i in 1:g_num){
    #computing gi which is percentage of Gi
    tmp_g <- groupdf[with(groupdf, which(groupdf$groupnum == groupvals[i])),"x"]/dataN
    #row_names[i] <- groupdf$groupnum[i]
    #print(tmp_g)
    #colindex refferring to #columns
    for(j in 1:c_num){
      #computing j which is percentage of Cj
      tmp_c <- clusterdf[with(clusterdf,which(clusterdf$clusternum == clustervals[j])),"x"]/dataN
      #col_names[j] <- clusterdf$clusternum[j]
      #print(paste("tmp_c",j,tmp_c))
      #probability of gi and cj
      tmp_prob <- tmp_g*tmp_c
      #print(paste("tmp_prob",j,tmp_prob))
      #count of occurence of Gi in Cj which is intersection
      tmp <- grpclusterdf[with(grpclusterdf,
                               which(grpclusterdf$groupnum==groupvals[i] & 
                                 as.numeric(grpclusterdf$clusternum)==clustervals[j])),"x"]
      #populate the prob. matrix byrow
      tmpMatrix_1[i,j] <- tmp_prob
      tmpMatrix_2[i,j] <- floor((tmp_prob*dataN)+0.5) #expected matrix is also converted to integer for comparison with observered
      #populate the #occurences matrix byrow
      if(length(tmp)>0)
        tmpMatrix[i,j] <- tmp
      else
        tmpMatrix[i,j] <- 0
    }
  }
  #print(row_names)
  #print(col_names)
  dimnames(tmpMatrix) <- list(row_names, col_names)
  dimnames(tmpMatrix_1) <- list(row_names, col_names)
  dimnames(tmpMatrix_2) <- list(row_names, col_names)
  print(tmpMatrix)
  return (list(tmpMatrix, tmpMatrix_1, tmpMatrix_2))
}

#function to create a group-cluster combination for significant associations
getSignificantGroupsNClusters <- function(sigcond, maindf, clustergrp){
  if(is.null(maindf$mxdgroupnum) | is.null(maindf[,clustergrp])){
    cat("Not a valid dataframe, mxdgroupnum and classification column were not found!")
    return
  }
  t <- NULL
  for(i in 1:nrow(sigcond)){
    t <- rbind(t, maindf[with(maindf,
                              which(maindf$mxdgroupnum==sigcond[i,1] & 
                                      maindf[,clustergrp]==sigcond[i,2])),
                         c("archive_name","ct_groupttl","br_groupttl","sp_groupttl","protocol","slice_thickness")][1,])
  }
  #print(dim(t))
  sigcond <- cbind(sigcond, t,row.names=NULL)
  #sigcond <- merge(x=sigcond,y=t,by.x=c('row','col'),by.y=c('mxdgroupnum','classification'))
  #rownames(sigcond) <- NULL
  #print(dim(sigcond))
  #print(colnames(sigcond))
  sigcond <- sigcond[order(sigcond$row,sigcond$col),]
  #print(sigcond)
  return (sigcond)
}

#function that calculates counts and percentages of grp-cluster associations
calculateGCProportions <- function(sigMatrix, observedMat){
  t1 <- data.frame()
  for(i in 1:nrow(sigMatrix)){
    #calculate the percentage of neurons that are from a group given the cluster.  
    ginc <- round(observedMat[sigMatrix$row[i],sigMatrix$col[i]]/sum(observedMat[sigMatrix$row[i],]),2)
    #calculate the percentage of neurons that are from a cluster given the group.
    cing <- round(observedMat[sigMatrix$row[i],sigMatrix$col[i]]/sum(observedMat[,sigMatrix$col[i]]),2)
    t1 <- rbind(t1,c(ginc,cing, observedMat[sigMatrix$row[i],sigMatrix$col[i]]))
  }
  
  length(t1)
  print(dim(t1))
  dim(sigMatrix)
  colnames(t1) <- c("grpperc", "clusperc","counts")
  sigMatrix <- cbind(sigMatrix, t1)
  dim(sigMatrix)
  colnames(sigMatrix)
  #re-organize the order of columns
  sigMatrix <- sigMatrix[,c(1:2,7:9,3:6)]
  sigMatrix
  rownames(sigMatrix) <- NULL
  return (sigMatrix)
}


#compute the distances between two clusters and plot cluster plots for visualization
compare2clusters <- function(subsetclus,problog10Matrix,metagrpsList){
  grpswith2clusters <- list()
 
  #initialize the distmatrix
  distMatrix <- data.frame(grpNum=character(0),metagrp=character(0),PC1=numeric(0),PC2=numeric(0),stringsAsFactors=FALSE)
  str(distMatrix)
  #loop for each group that is divided into two clusters
  for(i in 1:nrow(subsetclus)){
    #the two highly significant clusters
    c1 <- which(problog10Matrix[subsetclus$Group.1[i],]>0)[1]   #get the two cluster numbers that the group is divided into
    c2 <- which(problog10Matrix[subsetclus$Group.1[i],]>0)[2]
    print("***********************************")
    print(paste("checking clusters:",c1,c2))
    tmpclustrData <- subset(wholeReduced,wholeReduced$mxdgroupnum == subsetclus$Group.1[i] & wholeReduced$classification %in% c(c1,c2))
    #tmpMetricData <- subset(datamatrix,rownames(datamatrix) %in% rownames(tmpclustrData))
    no_of_tests <- 0
    #check for confounds in each group
    retlist <- checkforConfounds(tmpclustrData,metagrpsList)
    print(retlist)
    pvals <- c()
    ctables <- c()
    confounDtables <- c()
   
    if(length(retlist)>0){
      #store the pvals, contingency tables and confounds for each group
      for(j in 1:length(retlist)){
        metaT <- retlist[[j]]$metaD$metatype
        metaV <- retlist[[j]]$metaD$metavalues
        confDList <- retlist[[j]]$confounD
        no_of_tests <- no_of_tests + length(metaV)
        #perform chisq test on the non-confound metatypes and check if they are significant
        chitest <- calculateChisqpvals(metaT, tmpclustrData, "classification")
        
        if(chitest$pval <= alpha){
          #create an array of pvalues with names as metatypes
          pvals <- append(pvals,chitest$pval)
          names(pvals)[length(pvals)] <- metaT #add metatype as name to the last element of the array
          #create array of contingency tables
          ctables <- append(ctables, list(chitest$ctable))
          names(ctables)[length(ctables)] <- metaT
          #create array of confound meta types
          confounDtables <- append(confounDtables, list(retlist[[j]]$confounD))
          names(confounDtables)[length(confounDtables)] <- metaT
        }   
      }
      #order pvals list and do the multi-test correction on the corresponding contingency tables excluding the confounds
      if(length(pvals)>0){
        pvals <- pvals[order(pvals)]
        print(paste("metaT:",metaT))
        print("pvals:")
        print(pvals)
        #do FDR correction and bonferroni correction
        for(j in 1:length(pvals)){
          metatyp <- names(pvals)[j]
          ind <- match(metatyp,names(ctables))
          contintable <- ctables[[ind]]
          ind <- match(metatyp,names(confounDtables))
          conftable <- confounDtables[[ind]]
          #print(contintable)
          FDR <- no_of_tests/j
          Bf <- (nrow(contintable) - 1)*(ncol(contintable)-1)
          print(paste('Bf:', Bf, 'correction:',FDR * Bf))
          contintable <- round(contintable * FDR * Bf, 2)
          print("contingency table:")
          print(contintable)
          sigind <- which(contintable<=alpha,arr.ind=T) #check cells with significant values
          #print(unique(row.names(sigind)))
          #for groups that have atleast 2 significant meta-values plot cluster 
          if(length(unique(row.names(sigind)))>1){
            sigGrpsubset <- subset(tmpclustrData, tmpclustrData$classification %in% c(c1,c2) & 
              tmpclustrData[,metatyp] %in% row.names(sigind))
            #sigGrpMetricData <- subset(tmpMetricData, rownames(tmpMetricData) %in% rownames(sigGrpsubset))
            
            print(paste(metatyp, "is significant after correction"))
            if(length(conftable)>0){
              for(c in 1:length(conftable)){
                print(conftable[[c]]$conftypes)
                print(conftable[[c]]$confvals)
              }
            }
            #store the details of significant groups
            grpswith2clusters <- append(grpswith2clusters,list(pvalue = pvals, probMat = ctables, 
                                                               confoundTable = confounDtables,subgrpDF = sigGrpsubset))
          
            #plot 2D cluster plot
            #show only metavals that are significant
            featurecol <- sigGrpsubset[,metatyp]
            #print(unique(featurecol))
            colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
            grpN <- sigGrpsubset$mxdgroupnum[1]
            grpttl <- paste(sigGrpsubset$ct_groupttl[1], sigGrpsubset$br_groupttl[1],sigGrpsubset$sp_groupttl[1])
            Xname <- 'PC1'
            Yname <- 'PC2'
            #print(sigGrpsubset[,Xname])
            if(length(featurecol)>=2){
              PCOrderDf <- findAxisWLargestSeparation(sigGrpsubset, metatyp)            
              
              Xname <- as.character(PCOrderDf$PC[1]) #first highest separation
              Yname <- as.character(PCOrderDf$PC[2]) #second highest separation
            }
            
            plot(sigGrpsubset[,Xname],sigGrpsubset[,Yname],col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
            text(sigGrpsubset[,Xname],sigGrpsubset[,Yname],labels=as.character(sigGrpsubset$classification),col=mypalette[colarr],cex=.7)  
            legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))
            #grpEM <- Mclust(newPCDf[,33:64],G=1:7)
            #print(grpEM)
            #draw ellipses based on original classification
            clusnums <- unique(sigGrpsubset$classification)
            parms <- tmpEM_95$parameters
            for(i in 1:length(clusnums)){
              XYpcs <- sigGrpsubset[with(sigGrpsubset,which(sigGrpsubset$classification == clusnums[i])),c(Xname,Yname)]
              XYmeans <- colMeans(XYpcs)
              
              #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
               #       radius=1, center.pch=as.character(clusnums[i]),col=1)
              ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
                             radius=1, center.pch=as.character(clusnums[i]),col=1)
            } 
        }
      }
        
    }
    
  
        #title of the group
        #print(tmpclustrData[1,c("mxdgroupnum","archive_name","ct_groupttl","br_groupttl","sp_groupttl")])
        #testing clusters from PCs
        #tmpEM <- clusteringwithEM(tmpclustrData[,grep("[PC]", names(tmpclustrData), value=TRUE)])
        #summary(tmpEM)
        #tmpclustrData$classification = tmpEM$classification
        #tmpbinom <- createGroupClusterMatrixUsingBinomialTest(metagrpsList[j], tmpclustrData)
        #singleMetadatabinomList <- append(singleMetadatabinomList,list(metagrpinfo,tmpbinom))
  }
  #colnames(datamatrix)
  #distMatrix <- distMatrix[order(distMatrix$PC1, distMatrix$PC2, decreasing=T), ]
  #print(distMatrix)
  }
  return (grpswith2clusters)
}


#find the feature that shows greatest distance between 2 clusters through computing the distance between centroids
findAxisWLargestSeparation <- function(featureDf,featureType){
  #make sure only principal components are considered
  Df <- featureDf[,grep("[PC]", names(featureDf), value=TRUE)]
  featureValues <- unique(featureDf[,featureType])
  print(featureValues)
  retlist <- data.frame(g1='',g2='',PC='',diff=0)
  retlist <- retlist[-1,]
  #the feature must have atlest two values
  if(length(featureValues)>=2){
    #check each component
    for(p in 1:ncol(Df)){
      name <- colnames(Df)[p]
      #compute the distance between two groups for each component
      g1 <- subset(featureDf[,name], featureDf[,featureType] == featureValues[1])
      g2 <- subset(featureDf[,name], featureDf[,featureType] == featureValues[2])
      diff <- abs(mean(g1) - mean(g2))
      #print(diff)
      retlist <- rbind(retlist,data.frame(g1=featureValues[1],g2=featureValues[2],PC=name,diff=diff))
    }
  }
 #sort the data frame on diff
  retlist <- retlist[order(-retlist$diff),]
  #print(retlist)
 return (retlist)
}

#check for confounds, returns a list with group and its confounded metatypes
checkforConfounds <- function(groupdf,metatypes){
  metatypes2 <- metatypes
  print(paste("evaluating metatypes for group:",unique(groupdf$mxdgroupnum)))
  skipmetatypes <- c("Not reported", "M/F") 
  #keep metatypes that has atleast 2 meta values
  for(i in 1:length(metatypes)){
    #print(length(unique(groupdf[,metatypes[i]])))
    #keep metatypes that have alteast two values and none of them are 'Not Reported'
    skip <- skipmetatypes %in% unique(groupdf[,metatypes[i]])
    if(length(unique(groupdf[,metatypes[i]]))>1 & length(skip[skip==TRUE])==0){
      print(paste("include metatype:", metatypes[i]))
    }else{#remove metatypes that has only one meta value
      ind <- which(metatypes2==metatypes[i])
      metatypes2 <- metatypes2[-ind]
    }
  }
  
  #the return list
  grps2ret <- list()
  confmetatypes <- list()
  #check for confounds 
  print(length(metatypes2))
  if(length(metatypes2)>0){
    print(metatypes2)
    for(i in 1:length(metatypes2)){
      maintype <- metatypes2[i]
      print(maintype)
      cnfndList <- list()
      #if(all(metatypes2 %in% chkdvals)==F){
      #metatype to the return list, only if it is not already confound with others.
      if(!(maintype %in% confmetatypes) & (i+1) <= length(metatypes2)){
        
        for(j in (i+1):length(metatypes2)){
          #print(length(metatypes2))
          pairtype <- metatypes2[j]
          print(pairtype)
          randI <- adjustedRandIndex(groupdf[,maintype],groupdf[,pairtype])
          print(paste("randI", randI))
          if(randI == 1){
            cnfndList <- append(cnfndList, list(list(conftypes = pairtype, confvals = unique(groupdf[,pairtype]))))
            #remove the confounded metatype as noticed
            ind <- which(metatypes2==pairtype)
            print(paste("confounded",ind, metatypes2[ind]))
            #metatypes2 <- metatypes2[-ind]
            confmetatypes<- append(confmetatypes,metatypes2[ind])  
            #if last element of the list and not a confound, then add to the main list
          }
        }
        #print("adding..")
        t <- list(grpN = unique(groupdf$mxdgroupnum), 
                  metaD = list(metatype = maintype, metavalues = unique(groupdf[,maintype])), 
                  confounD = cnfndList)
        grps2ret <- append(grps2ret,list(t))
      }else if(!(maintype %in% confmetatypes) & (i+1) > length(metatypes2)){#there is only one metatype with 2 metavalues
        t <- list(grpN = unique(groupdf$mxdgroupnum), 
                  metaD = list(metatype = maintype, metavalues = unique(groupdf[,maintype])), 
                  confounD = cnfndList)
        grps2ret <- append(grps2ret,list(t))
      } 
    }
  }
  return (grps2ret)
}

calculateChisqpvals <- function(metagrp, subsetDf, clustergrp){
  
  grps <- aggregate(subsetDf$neuron_name, by=list(groupnum=subsetDf[,metagrp]), FUN=length)
  print(grps)
  cols <- aggregate(subsetDf$neuron_name, by=list(clusternum=subsetDf[,clustergrp]), FUN=length)
  print(cols)
  grpvals <- grps$groupnum #unique(featPCDf[,metagrp])
  clusvals <- cols$clusternum #unique(featPCDf$classification)
  #determine the #group values in a single metadata feature
  g_num <- length(groupvals)#length(unique(featPCDf[,metagrp]))
  c_num <- length(clustervals)#nrow(clustermatrix)
  sz <- nrow(subsetDf)
  #metagrp aggregate on unique group values
  grpclus <- aggregate(subsetDf$neuron_name, by=list(groupnum=subsetDf[,metagrp],
                                                    clusternum=subsetDf[,clustergrp]), FUN=length)
  #compute observed and expected matrices
  tmp <- getObservedAndExpectedMatrices(grps,cols,grpclus,sz)
  observedM <- tmp[[1]]
  chisqtest <- chisq.test(observedM)
  print(paste("chi sql table pval:",chisqtest$p.value))
  stdres <- chisqtest$stdres        #standardized residuals 
  observedMchi <- chisqtest$observed
  expectedMchi <- floor(chisqtest$expected+0.5)
  
  probMatrixchi <- pnorm(-1 * abs(stdres)) * 2  #compute p-vals by reversing to the -ve side and multiplying by 2
  probMatrixchi[which(probMatrixchi==0,arr.ind=T)] <- pnorm(-37.5)*2 #minimum value for which p-val can be calculated
  
  chisqres <- list(pval =chisqtest$p.value, ctable =probMatrixchi)
  #probMatrixchi <- probMatrixchi * bonferroni_across #correction for number of metagrps 
  problog10Matrixchi <- converttolog10(probMatrixchi,observedMchi,expectedMchi/sz,sz)   #convert pvals to log10 values
  #print(problog10Matrixchi)
  #print(dim(problog10Matrixchi))
  #print(nrow(problog10Matrixchi))
  #print(str(problog10Matrixchi))
  #chisqres <- list(pval =chisqtest$p.value, ctable =problog10Matrixchi)
  
  return (chisqres)
}



clusterPlotFor2grps <- function(grpswith2clusters){
  ####plot individual cluster plots for groups that are divided into two parts###
  for(i in seq(from=1, to=length(grpswith2clusters), by=2)){
    metagrp <- grpswith2clusters[[i]]
    subgrpPCs <- grpswith2clusters[[i+1]]
    #tmpobserved <- binomtst$observedMat
    #tmpexpected <- binomtst$expectedMat
    #tmpprob <- binomtst$probMat
    #tmpdf <- binomtst$featPCDf
    #if(nrow(which(tmpprob>0,arr.ind=T))>1){
    #print(tmpobserved)
    #print(tmpexpected)
    #print(tmpprob)
    clusnums <- unique(subgrpPCs$classification)
    print(paste("comparing..",clusnums))
    featurecol <- subgrpPCs[,metagrp$singlemetadata]
    print(metagrp)
    dim(subgrpPCs)
    #arrange color palette
    colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
    print(unique(colarr))
    colnames(subgrpPCs)
    nrow(subgrpPCs)
    par(mfrow=c(1,2))
    aggregate(subgrpPCs[,metagrpsList],by=list(subgrpPCs$cellclass3),length)
    adjustedRandIndex(featurecol, subgrpPCs$classification)
    #plot cluster
    #plotcluster(subgrpPCs[,33:34],clvecd=subgrpPCs$classification,cex=0.8,col=mypalette[colarr],main=paste(metagrp$grpNum,metagrp$grpttl),xlab="PC1",ylab="PC2",cex.main=0.8)
    plot(subgrpPCs[,33],subgrpPCs[,34],col="white",cex=.55,pch=16,xlab=names(subgrpPCs)[33],main=paste(metagrp$grpNum,metagrp$grpttl),ylab=names(subgrpPCs)[34],cex.main=0.8)
    text(subgrpPCs[,33],subgrpPCs[,34],labels=as.character(subgrpPCs$classification),col=mypalette[colarr],cex=.7)  
    legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))
    parms <- tmpEM_95$parameters
    for(i in 1:length(clusnums)){
      ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
              radius=1, center.pch=as.character(clusnums[i]),col=1)
    }
    reclustersingleGroup(metagrp$grpNum,metagrp$singlemetadata,clusnums)
  }
}

#individual tests on single groups that are divided into two clusters
reclustersingleGroup <- function(grpnum,metafeat,clusterpair){
  t <- subset(wholeReduced, wholeReduced$mxdgroupnum==grpnum & wholeReduced$classification %in% clusterpair)
  print(dim(t))
  t$mxdgroupnum
  print(unique(t[,metagrpsList]))  
  grp6em <- Mclust(t[,33:64],G=1:5)
  print(grp6em)
  t$classification <- grp6em$classification
  colnames(t)
  paste("comparing",paste(unique(t[,metafeat],collapse=","), paste(unique(t$classification),collapse=",")))
  print(adjustedRandIndex(t[,metafeat], t$classification))
  colarr <- colorcode(groupnum=t[,metafeat], length(unique(t[,metafeat])))
  plot(t[,33],t[,34],col="white",cex=.6,pch=16,xlab=names(t)[33],main=paste("grp#",grpnum,"in clusters:",paste(unique(grp6em$classification),collapse=",")),ylab=names(t)[34],cex.main=0.9)
  text(t[,33],t[,34],labels=as.character(t$classification),col=mypalette[colarr],cex=.7)
  legend("bottomleft", text.width=strwidth("1,000"), legend = unique(t[,metafeat]), pch=19,cex=0.7,col=unique(mypalette[colarr]))
  for(i in 1:length(unique(grp6em$classification))){
    ellipse(grp6em$parameters$mean[1:2,i],shape=grp6em$parameters$variance$sigma[1:2,1:2,i],radius=1, center.pch=as.character(i),col=1)
  }
}







