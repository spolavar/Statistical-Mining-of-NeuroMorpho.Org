library(tourr)
library(fpc)
library(cluster)

source("C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/plotutilfunctions.R")

#arbor type = 3 & 4
dim(metricDataApiBas)
colnames(metricDataApiBas)
metricDataApiBas[which(metricDataApiBas$neuron_name == "LY30-RGC28.CNG.swc"),]
apicalbasal <- metricDataApiBas

metricDendrites <- cntNRmvNAs(apicalbasal,113)
colnames(metricDendrites)
dim(metricDendrites)

#### Dendrite only data preprocessing ####
#check for zero depth data or flat 2D neurons and remove them from analysis
flatDendneuritespos <-  with(metricDendrites, which(is.na(metricDendrites$Depth_total_sum)))
flatDendneurites <- metricDendrites[flatDendneuritespos, c("neuron_name","N_bifs_total_sum")]
print(paste("removing ",length(flatDendneuritespos)," from", nrow(metricDendrites), " neurons"))
metricDendrites <- metricDendrites[-flatDendneuritespos,]

#remove CNG.swc from neuron_name
metricDendrites <- clnData(metricDendrites)

dim(metricDendrites)
#a subset of metrics to be included in pca
#dendpcaMetrics <- pcaMetrics[-2]
#merge metric data with meta data
length(dendpcaMetrics)
metricDendrites <- getPlotSubset(metricDendrites,metrics = dendpcaMetrics)
colnames(metricDendrites)
unique(metricDendrites$cellclass2)
dim(metricDendrites[which(metricDendrites$cellclass2 == "Basket cell" &  metricDendrites$archive_name == "Markram"),])
unique(metricDendrites$archive_name)
#usedMetaData[which(usedMetaData$archive_name == "Rodger"),]
#remove extra column that may have been inserted from merge process
metricDendrites <- metricDendrites[, !(colnames(metricDendrites) %in% c("X"))]
#remove neurons with < 4 bifurcations
Nbranchoutlierpos <- with(metricDendrites, which(metricDendrites$N_bifs_total_sum < 4))
Nbranchoutlier <- metricDendrites[Nbranchoutlierpos,c("neuron_name","N_bifs_total_sum","archive_name")]

#Nbranchoutlierpos <- with(metricDendrites, which(metricDendrites$N_bifs_total_sum + metricDendrites$N_tips_total_sum < 5))
#Nbranchoutlier <- metricDendrites[Nbranchoutlierpos,c("neuron_name","N_bifs_total_sum","archive_name")]
str(Nbranchoutlier)
metricDendrites <- metricDendrites[-Nbranchoutlierpos,]
#remove Larkman archive because his neurons are dendrograms and hence aren't valid biological metrics.
larkmanpos <- with(metricDendrites, which(metricDendrites$archive_name=="Larkman"))
larkmanoutlier <- metricDendrites[larkmanpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(larkmanoutlier)
metricDendrites <- metricDendrites[-larkmanpos,]
#remove neurons with fractal==0
zerofractalpos <- with(metricDendrites, which(metricDendrites$Fractal_Dim_avg == 0))
zerofractaloutlier <- metricDendrites[zerofractalpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(zerofractaloutlier)
metricDendrites <- metricDendrites[-zerofractalpos,]
#check for zero Bif ampl remote angle or Bif torque local angle
zeroBifamplremote <- with(metricDendrites, which(metricDendrites$Bif_torque_local_avg == 0 ))
zeroBifamplremote <-  metricDendrites[zeroBifamplremote,c("neuron_name","archive_name")]
#nrow(zerosomaoutlier)
#metricDendrites <- metricDendrites[-zerosomasurface,]
#hist(metricDendrites$Helix_max, breaks=30, xlab="Helix")
#hist(metricDendrites$Fractal_Dim_avg, breaks=30, xlab="fractal")
min(metricDendrites$Helix_max)
helix_zeropos <- with(metricDendrites, which(metricDendrites$Helix_max == 0.00))
helixzero <- metricDendrites[helix_zeropos,c("neuron_name","Fractal_Dim_avg","Contraction_avg","Branch_pathlength_max","archive_name")]
helixzero
##########Data transformation & standardization#############
plotDf <- plotTesting(metricDendrites)
 
dim(plotDf)
colnames(plotDf)
plotDf$Bif_ampl_remote_max
metricDendrites$Bif_ampl_remote_max
#standardize metric dataframe before computing PCA
for(i in 33:63){
  #get non-NA values for normalization
  colv <- plotDf[,i][!is.na(plotDf[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  plotDf[,i][!is.na(plotDf[,i])] <- colv/msd
  print(paste(names(plotDf)[i],msd))
}

#### create cell types groups for dendrites ####
aggregate(plotDf$expercond, by = list(plotDf$expercond), length)
#celltype hierarchy
primary <- list(expercond="Control")
x <- getSubset2(plotDf, primary, typeMetaData=usedMetaData)
level1 <- "cellclass1"
level2 <- "cellclass2"
level3 <- "cellclass3"
#make groups following hierarchical ordering of metadata
ctMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData, level1, level2)#level3) #,minSize=100)
printGrps(ctMetricHierarchy)
colnames(plotDf)
unique(plotDf$archive_name)
#group neuron_name (primariy id) according the hierarchy lists
Grplst <- addgrps2MetricDF(plotDf,ctMetricHierarchy)
for(i in 1:length(Grplst)){
  print(Grplst[[i]]$GroupSpec)
}
#assigning a groupnum and groupttl to each neuron according to the hierarchy list
plotDf_ct <- metricDf_plotting(plotDf,Grplst)

#### create brain regions groups for dendrites ####
primary <- list(expercond = "Control")

level1 <- "region1"
level2 <- "lobes"
level3 <- "region2"
level4 <- "region3"
#unique(usedMetaData$lobes)
unique(usedMetaData$lobes)
brMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData, level1, level2,level3) #level4)#minSize=300)

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

#### create species groups for dendrites ####
primary <- list(expercond = "Control")
level1 <- "order"
level2 <- "species"
level3 <- "strain"
spstrainMetricall <- makeHierarchyGroups(primary, plotDf,usedMetaData,level1,level2)#,minSize=55)
printGrps(spstrainMetricall)
Grplst <- addgrps2MetricDF(plotDf,spstrainMetricall)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}
aggregate(plotDf$expercond, by = list(plotDf$expercond), length)
plotDf_sp <- metricDf_plotting(plotDf,Grplst)
dim(plotDf_sp)
colnames(plotDf_sp)

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
colnames(plotDf_ct)
colnames(plotDf_ct[,33:63])
#handles missing values by removing rows that have NAs - remove soma_surface from PCA analysis as not all data have soma_surface
pcametrics <- na.omit(plotDf_ct[,33:63])
dim(pcametrics)
pcametrics <- subset(plotDf_ct, plotDf_ct$neuron_name%in%rownames(pcametrics))
colnames(pcametrics)
names(pcametrics)[64] <- "ct_groupnum"
names(pcametrics)[65] <- "ct_groupttl"

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
write.table(fourdims, file="DendriteOnlyGroupsBeforeEliminating.txt",quote=FALSE, sep=",")

#groups from 4D space combining archives/species/celltypes/brain regions
mixedgrpsDf <- subset(fourdims,fourdims$x>=40)#55)
dim(mixedgrpsDf)
mixedgrpsDf
#prepare data matrix for computing pca and for cluster analysis
Dendritedatamatrix <- addmixedgrpcol(mixedgrpsDf, pcametrics)
dim(Dendritedatamatrix)

mxdgrpctlg <- aggregate(Dendritedatamatrix$neuron_name, as.lists(Dendritedatamatrix[,c("archive_name","ct_groupttl","br_groupttl","mxdgroupnum")]), FUN = length)
write.table(mxdgrpctlg, file="mxdgrpTable4Dendrites.txt",quote=TRUE, sep=",")
colnames(Dendritedatamatrix)
#check for correlation before running PCA
prePCAcorr_dend <- mostcorrelated(Dendritedatamatrix[,34:59],threshold=0.5)
str(prePCAcorr_dend)
colnames(prePCAcorr_dend)
#eliminate one of the highly correlated pairs
removeFeatures <- eliminateCorrFeat(prePCAcorr_dend)
#remove the selected correlated features from data matrix
Dendritedatamatrix <- Dendritedatamatrix[ , -which(names(Dendritedatamatrix) %in% removeFeatures$feature)]
colnames(Dendritedatamatrix)
dim(Dendritedatamatrix)

#get pca matrix from the log transformed and normalized dendrite metrics
pcamat <- prcomp(x=Dendritedatamatrix[,33:59],center=T,scale=T)
pcasumm <- summary(pcamat)
#make scree plot to choose most contributing PCs
plot(pcamat$sdev^2/sum(pcamat$sdev^2)*100, type = "b", ylim = c(0,100),xlab="PCA Eigenvalues",ylab="% of variance") 
     #main = "Scree plot for dendrite (apical & basal) PCs",xlab="PCA Eigenvalues",ylab="% of variance")
#plot cumulative variance
lines(pcasumm$importance[3,]*100,main = "whole arbor metrics",type = "b", col = 2, xlab="cumulated PCA eigen values",ylab="% of variance")
legend(10,50,c("cumulative variance","variance"),lty=c(1,1),col=c("red","black"))

#95% of cumulative variance is at
pcasumm$importance[3,]*100

#six components were covering ~75% of variance
#choose eigenvalues >= 1
pcamat$sdev ^ 2

#writing the highest loadings on each component into a file
colnames(pcametrics)
loaddend <- abs(pcamat$rotation)
dim(loaddend)
PCloadingsarr <- matrix(nrow=27,ncol=3)
dim(PCloadingsarr)
for(c in 1:27){
  rowid <- which(loaddend == max(loaddend[,c]),arr.ind=TRUE)[1]
  colid <- which(loaddend == max(loaddend[,c]),arr.ind=TRUE)[2]
  metricname <- names(loaddend[,colid])[rowid]
  PCname <- names(loaddend[rowid,])[colid]
  PCloadingsarr[c,1] <- PCname
  PCloadingsarr[c,2] <- metricname
  PCloadingsarr[c,3] <- pcamat$rotation[rowid, colid]
}

PCloadingsarr
write.table(PCloadingsarr, file="PCloadingsarr_dend.txt",quote=FALSE, sep=",")

#loadings of PC1 and PC2
pc1dendloadings <-sort(pcamat$rotation[,1],decreasing=F)
pc1contr <- pc1dendloadings[pc1dendloadings > 0.2 | pc1dendloadings < -0.2]
dotchart(pc1contr,main="PC1 loadings",cex=0.7, xlab="variable loadings",col=2)
pc2dendloadings <-sort(pcamat$rotation[,2],decreasing=F)
dotchart(pc2dendloadings,main="PC2 loadings",cex=0.7, xlab="variable loadings",col=2)
pc2contr <- pc2dendloadings[pc2dendloadings < -0.32]
dotchart(pc2contr,main="PC2 loadings",cex=0.7, xlab="variable loadings",col=2)
pc3dendloadings <-sort(pcamat$rotation[,3],decreasing=F)
dotchart(pc3dendloadings,main="PC3 loadings",cex=0.7, xlab="variable loadings",col=2)
pc3contr <- pc3dendloadings[pc3dendloadings > 0.30 | pc3dendloadings < -0.35]
dotchart(pc3contr,main="PC3 loadings",cex=0.7, xlab="variable loadings",col=2)
pc4dendloadings <-sort(pcamat$rotation[,4],decreasing=F)
dotchart(pc4dendloadings,main="PC4 loadings",cex=0.7, xlab="variable loadings",col=2)
pc4contr <- pc4dendloadings[pc4dendloadings > 0.5]
dotchart(pc4contr,main="PC4 loadings",cex=0.7, xlab="variable loadings",col=2)
pc5dendloadings <-sort(pcamat$rotation[,5],decreasing=F)
pc5contr <- pc5dendloadings[pc5dendloadings > 0.3 | pc5dendloadings < -0.3]
dotchart(pc5dendloadings,main="PC5 loadings",cex=0.7, xlab="variable loadings",col=2)
pc6dendloadings <-sort(pcamat$rotation[,6],decreasing=F)
dotchart(pc6dendloadings,main="PC6 loadings",cex=0.7, xlab="variable loadings",col=2)
pc6contr <- pc6dendloadings[pc6dendloadings > 0.28 | pc6dendloadings < -0.38]
dotchart(pc6contr,main="PC6 loadings",cex=0.7, xlab="variable loadings",col=2)

pc7dendloadings <-sort(pcamat$rotation[,7],decreasing=F)
dotchart(pc7dendloadings,main="PC7 loadings",cex=0.7, xlab="variable loadings",col=2)
pc7contr <- pc7dendloadings[pc7dendloadings > 0.3 | pc7dendloadings < -0.3]
dotchart(pc7contr,main="PC7 loadings",cex=0.7, xlab="variable loadings",col=2)

pc8dendloadings <-sort(pcamat$rotation[,8],decreasing=F)
#pc8loadings <-sort(abs(pcamat$rotation[,8]),decreasing=F)
dotchart(pc8dendloadings,main="PC8 loadings",cex=0.7, xlab="variable loadings",col=2)

#create a data frame that has Principal components and metadata groups
dim(Dendritedatamatrix)
colnames(Dendritedatamatrix)
dim(pcamat$x)
DendReduced <- data.frame(Dendritedatamatrix[,1:32], pcamat$x, Dendritedatamatrix[,60:66])
dim(DendReduced)
unique(DendReduced$ct_groupnum)
length(unique(DendReduced$ct_groupnum))
unique(DendReduced$br_groupnum)
length(unique(DendReduced$br_groupnum))
unique(DendReduced$sp_groupnum)
length(unique(DendReduced$sp_groupnum))
unique(DendReduced$mxdgroupnum)

#### clustering the PCA space using EM ####
colnames(DendReduced)
dim(DendReduced)
unique(DendReduced$mxdgroupnum)
#cluster PCs of feature matrix on 18 features (upto 95% of cumulative variancce)
#Mclust(DendReduced[,grep("[PC]", names(DendReduced), value=TRUE)],G=1:15)
tmpdendEM <- Mclust(DendReduced[,grep("[PC]", names(DendReduced), value=TRUE)],G=1:15)
  #Mclust(DendReduced[,33:50],G=1:15, model ='VVV')
tmpdendEM_95 <- Mclust(DendReduced[,33:49],G=1:15)
summary(tmpdendEM_95)
#tmpdendEM_95 <- Mclust(DendReduced[,33:49],G=1:15,model=c('VEV','VVV'))
#summary(tmpdendEM_95)
plot(tmpdendEM_95, data=wholeReduced[,33:49], what='BIC')
tmpdendEM_95$bic
tmpdendEM_95$BIC

#tmpdendEM_95 <- Mclust(DendReduced[,33:49],G=1:15,model='VVV')
#summary(tmpdendEM_95)
#tmpdendEM_95$bic
#tmpdendEM_95$BIC

#emResult <- Mclust(clstrMatrix,G=1:15)
summary(tmpdendEM)
x11()
#choose only first 18 components for clustering
plot(tmpdendEM_95, data=DendReduced[,33:50], what='BIC')
plot(tmpdendEM_95, what="classification",dimens=1:2)
plot(tmpdendEM_95, what="uncertainty",dimens=1:2)
#plot(tmpdendEM_95, what="density",dimens=1:2)
unique(tmpdendEM_95$classification)
DendReduced$classification = tmpdendEM_95$classification
Dendritedatamatrix$classification = DendReduced$classification
#############Find significant group-cluster associations from binomial test##############

alpha <- 0.05
#takes input parameters: single metadata feature and the original data frame with PCs.
grpfeature <- 'mxdgroupnum'
#add cluster classification to featPCDf
groupmatrix <- aggregate(DendReduced$neuron_name, by=list(groupnum=DendReduced[,grpfeature]), FUN=length)
print(groupmatrix)
clustermatrix <- aggregate(DendReduced$neuron_name, by=list(clusternum=DendReduced$classification), FUN=length)
print(clustermatrix)

grpclustermatrix <- aggregate(DendReduced$neuron_name, by=list(groupnum=DendReduced[,grpfeature],
                                                                clusternum=DendReduced$classification), FUN=length)
groupvals <- groupmatrix$groupnum #unique(featPCDf[,grpfeature])
clustervals <- clustermatrix$clusternum #unique(featPCDf$classification)

row_names <- groupvals
col_names <- clustervals
#determine the #group values in a single metadata feature
g_num <- length(groupvals)#length(unique(featPCDf[,grpfeature]))
#run EM based clustering on featPCDf only on PC columns
#determine the #clusters from that
c_num <- length(clustervals)#nrow(clustermatrix)
dataN <- nrow(DendReduced)
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
#check if atleast 80% of expected table have counts > 5
g_num*c_num*0.8
length(expectedMat_chi[which(expectedMat_chi > 5)])

probMatrix_chi <- pnorm(-1 * abs(stdresMatrix_chi)) * 2  #compute p-vals by reversing to the -ve side and multiplying by 2
probMatrix_chi[which(probMatrix_chi==0,arr.ind=T)] <- pnorm(-37.5)*2 #minimum value for which p-val can be calculated

probMatrix_chi <- withinTableCorrectionBf(probMatrix_chi)
problog10Matrix_chi <- converttolog10(probMatrix_chi,observedMat_chi,expectedMat_chi/dataN,dataN)   #convert pvals to log10 values

###run binom test on #groups and #columns
#probMatrix <- runBinomtest(observedMat, expectedMat, dataN, g_num, c_num)
###convert pvals to log10(pvals)
#probMatrixlog10 <- converttolog10(probMatrix,observedMat,expectedMat,dataN)

###check for overshooting groups and undershooting on chisq table###
sigcondDend_over <- which(problog10Matrix_chi>0, arr.ind=TRUE)
#get groups and clusters indicies for all significant differences both over and under representation
sigcond_dend <- which(!is.na(problog10Matrix_chi), arr.ind=TRUE)
nrow(sigcond_dend)
str(sigcondDend_over)
sigcondDend_over[,'row']
sigcondDend_over[,'col']
nrow(sigcondDend_over)
sigMatrixDend <- getSignificantGroupsNClusters(sigcondDend_over,DendReduced,"classification")
sigMatrix_dend <- getSignificantGroupsNClusters(sigcond_dend,DendReduced,"classification")


###analyzing by groups: the ones that split in two groups###
library("car")
library("tourr")
sigMatrixDend <- calculateGCProportions(sigMatrixDend,observedMat_chi)


groupoccurences <- aggregate(sigMatrixDend$row, by=list(sigMatrixDend$row), FUN=length)
multiplegrps <- groupoccurences[with(groupoccurences, which(groupoccurences$x==2)),"Group.1"]
#select groups that are divided into 2 clusters

sigMatrix_s <- groupoccurences[groupoccurences$x==2,]

write.table(sigMatrixDend, file="OverRepresentedDendClusters.txt",quote=TRUE, sep=",",row.names=FALSE)
###print matrices to file####
#write the chisq test pval matrix to file
write.table(round(problog10Matrix_chi,2), file="Dendchisq_log10pval.txt",quote=FALSE, sep=",",row.names=TRUE)
#write data matrices into output files for analysis
write.table(sigMatrix_s, file="sigMatrix_s.txt",quote=FALSE, sep=",",row.names=FALSE)
#The observed matrix
observedMatrix
write.table(observedMat_chi, file="ObservedDendMatrix.txt",quote=FALSE, sep=",",row.names=FALSE)
#the expected matrix converted to integer for ease of comparison
expectedMatrix
write.table(expectedMat_chi, file="expectedDendMatrix.txt",quote=FALSE, sep=",",row.names=FALSE)
#the delta matrix
deltaMatrix <- observedMat_chi - expectedMat_chi
deltaMatrix
write.table(deltaMatrix, file="deltaDendMatrix.txt",quote=FALSE, sep=",",row.names=FALSE)
#the ratio matrix
dataN
ratioMatrix <- round(observedMat_chi/expectedMat_chi,2)
ratioMatrix
write.table(ratioMatrix, file="ratioDendMatrix.txt",quote=FALSE, sep=",",row.names=FALSE)

#plot all 7 clusters for whole arbor data
dim(DendReduced)
unique(DendReduced$mxdgroupnum)
colnames(DendReduced)
unique(DendReduced$archive_name)
colarr <- colorcode(groupnum=DendReduced$classification)
zscralphconvers <- c("f","e","a","c","d","b")
DendReduced[,"alphclusters"] <- zscralphconvers[DendReduced$classification]
DendReduced$alphclusters[1:10]
DendReduced$classification[1:10]


DendReduced_dup <- DendReduced
DendReduced_dup$PC2 <- -1 * DendReduced_dup$PC2

plot(DendReduced_dup[,33],DendReduced_dup[,34],col="white",cex=.55,pch=16,xlab=names(DendReduced_dup)[33],
     ylab=names(DendReduced_dup)[34],cex.main=0.7)#main="4904 dendrites only data classified into 6 models",
text(DendReduced_dup[,33],DendReduced_dup[,34],labels=DendReduced_dup$alphclusters,col=mypalette[colarr],cex=.8)

#change PC2 direction to position ellipses appropriately 
parms <- tmpdendEM_95$parameters
parms$mean[2,] <- -1 * parms$mean[2,]

clusnums <- unique(DendReduced_dup$classification)
clusalphas <- unique(DendReduced_dup$alphclusters)
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.cex = 2.5, center.pch=as.character(clusalphas[i]),col=1,lwd=2)
}

#########Not used#################

plot(tmpdendEM_95,what="classification",dimens=1:2,cex =0.8,pch=1)
summary(tmpdendEM)
tmpdendEM$classification
clusnums <- unique(DendReduced$classification)
parms <- tmpdendEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}


#plotting clusters with 
colarr <- colorcode(sigMatrix_s$Group.1,length(unique(sigMatrix_s$Group.1)))
grpscluster <- subset(DendReduced, DendReduced$mxdgroupnum%in%unique(sigMatrix_s$Group.1))
dim(grpscluster)
unique(grpscluster$mxdgroupnum)
#plot PCs in cluster plotting
colarr <- colorcode(groupnum=grpscluster$mxdgroupnum, length(unique(grpscluster$mxdgroupnum)))
plot(grpscluster[,33],grpscluster[,34],col="white",cex=.55,pch=16,xlab=names(grpscluster)[33],
     main="11 out of 37 groups that are divided between 2 clusters",ylab=names(newsubcolPCDf)[34],cex.main=0.8)
text(grpscluster[,33],grpscluster[,34],labels=as.character(grpscluster$classification),col=mypalette[colarr],cex=.7)
legend("topleft",text.width = strwidth("10"), legend = unique(grpscluster$mxdgroupnum), pch=19, cex=0.7, col=unique(mypalette[colarr]))
clusnums <- unique(grpscluster$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}




#test k=5 & 8 clusters based on the normalized scatter (atleast 2% drop)
kmeansObj <- kmeans(pcamat$x,centers=5,nstart=20,iter.max=30)
summary(kmeansObj)
str(kmeansObj)
str(pcamat$x)
unique(pcametrics[,c("groupttl","groupnum")])
rownames(pcamat$x)
#pos <- psp$groupttl%in%c("Control/Carnivora/Cat","Control/Amphibians/Salamander")
#plot all celltypes groups 
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=kmeansObj$cluster,main="k-means scatterplot Matrix")
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=nonpcaDf$groupnum,main="k-means scatterplot Matrix")
#sub sample of pyramidal cell others and interneuron others
pairs(~PC1+PC2+PC3+PC4,data=subset(pcamat$x,rownames(pcamat$x) %in% rownames(subset(nonpcaDf,nonpcaDf$groupnum%in%c(8,9,10)))),col=subset(nonpcaDf,nonpcaDf$groupnum%in%c(8,9,10))$groupnum,main="k-means scatterplot Matrix")
dendriteReduced <- data.frame(pcamat$x, cluster=kmeansObj$cluster, groupnum=nonpcaDf$groupnum, groupttl=nonpcaDf$groupttl)
colnames(dendriteReduced)

length(nonpcaDf$groupnum)
dim(pcamat$x)
dim(pcametrics)
colnames(pcametrics)
dim(nonpcaDf)
pcasub <- subset(nonpcaDf, rownames(nonpcaDf) %in% rownames(pcamat$x))

#plot model based plots
emResult <- Mclust(dendriteReduced[,1:31])
plot(emResult)

emResult <- Mclust(dendriteReduced[,1:31],G=1:31,modelNames="EII")
plot(emResult)

emResult <- Mclust(dendriteReduced[,1:3],G=1:5,modelNames="EII")
summary(emResult)
plot(emResult)


#grplst <- c(2,4,6,8,9,10)
#select a subgroup of cell types to plot on cluster plot
grplst <- c(2,4,6,8,9,10)
grplst2 <- c(5,1,7)#interneuron and other principal cells
pyrlst <- c(3)#pyramidal cells only
f <- subset(dendriteReduced,rownames(pcamat$x) %in% rownames(subset(nonpcaDf,nonpcaDf$groupnum%in%grplst)))
f2 <- subset(dendriteReduced,rownames(pcamat$x) %in% rownames(subset(nonpcaDf,nonpcaDf$ct_groupnum%in%grplst2)))
colnames(nonpcaDf)
unique(as.character(f2$groupttl))
unique(f2$groupnum)
dim(f)
dim(f2)
dim(nonpcaDf)
rownames(f)
metagrps <- subset(nonpcaDf,nonpcaDf$groupnum%in%grplst)$groupnum
unique(metagrps)
metalbls <- subset(nonpcaDf,nonpcaDf$groupnum%in%grplst)$groupttl
unique(metalbls)
plot.new()
animate(f[,1:3], grand_tour(), display_xy(col=factor(f$groupnum)))
legend("bottomleft", legend = sapply(strsplit(as.character(unique(metalbls)),"/"),tail,1), pch = 19, col=unique(factor(f$groupnum)))
f2$groupttl <- as.character(f2$groupttl)
f2<-editCelltypeLabels(f2,"groupttl")
animate(f2[,1:3], grand_tour(), display_xy(col=factor(f2$groupnum)))
legend("bottomleft", legend = as.character(unique(f2$groupttl)), pch = 19, col=unique(factor(f2$groupnum)))


#labels=length(unique(f$groupnum))
colnames(f)
colnames(hybrid3pcaDf)
unique(f$groupnum)
#distr of PC1
hist(f[,1], breaks=30, xlab="PC1")
#distr of PC2
hist(f[,2], breaks=30, xlab="PC2")
#distr of PC3
hist(f[,3], breaks=30, xlab="PC3")
#distr of PC4
hist(f[,4], breaks=30, xlab="PC4")
#loadings of PC1 and PC2
pc1loadings <-sort(pcamat$rotation[,1],decreasing=F)
dotchart(pc1loadings,main="PC1 loadings",cex=0.7, xlab="variable loadings",col=2)
pc2loadings <-sort(pcamat$rotation[,2],decreasing=F)
dotchart(pc2loadings,main="PC2 loadings",cex=0.7, xlab="variable loadings",col=2)
#plot celltype cluster and groups
plotcluster(f[,1:4],f$cluster,cex=0.8,col=factor(f$groupnum),main="clustering and groups projection plots")
legend("bottomright", legend = sapply(strsplit(as.character(unique(metalbls)),"/"),tail,1), pch=19, cex=0.7, col=unique(factor(f$groupnum)))
#plot subgroups for celltypes on brain regions
colnames(f)
unique(ff$lobes)

unique(ff1$mxdgroupnum)
aggregate(ff$mxdgroupnum,by=list(ff$lobes),length)
dim(ff)
dim(ff1)
dim(ff2)
dim(ff3)

ff1 <- subset(ff,ff$lobes=="Frontal lobe")
ff2 <- subset(ff,ff$lobes=="Parietal lobe")
ff3 <- subset(ff,ff$lobes=="Occipital lobe")
ff4 <- subset(ff,ff$lobes=="Hippocampus")
ff5 <- subset(ff,ff$lobes=="Insula")

plotcluster(ff1[,1:4],ff1$cluster,cex=0.8,col=factor(ff1$mxdgroupnum),main="clustering and groups projection plots")
plotcluster(ff2[,1:4],ff2$cluster,cex=0.8,col=factor(ff2$mxdgroupnum),main="clustering and groups projection plots")
plotcluster(ff3[,1:4],ff3$cluster,cex=0.8,col=factor(ff3$mxdgroupnum),main="clustering and groups projection plots")
plotcluster(ff5[,1:4],ff5$cluster,cex=0.8,col=factor(ff5$mxdgroupnum),main="clustering and groups projection plots")


legend("bottomright", legend = sapply(strsplit(as.character(unique(metalbls)),"/"),tail,1), pch=19, cex=0.7, col=unique(factor(f$groupnum)))


#plot subgroups for celltypes on brain regions
colnames(f)



keepcols <- c(1:4,32:34)

#choose same dataset as celltypes to plot brain region groups 
ctsubset <- subset(f[,keepcols],f$groupnum%in%grplst)
#choose first 4 PCs and celltype groupnum for remaing cell types
ctsubset2 <- f2[,keepcols]
colnames(f2)
#plot celltype groups for interneurons and other principal cells
plotcluster(f2[,1:4],f2$cluster,cex=0.8,col=factor(f2$groupnum),main="clustering and groups projection plots")
legend("bottomright", legend = unique(factor(f2$groupttl)), pch=19, cex=0.7, col=unique(factor(f2$groupnum)))

dim(ctsubset2)
dim(ctsubset)
colnames(ctsubset)
brcode <- subset(hybridpcaDf,hybridpcaDf$neuron_name%in%rownames(ctsubset))
brcode2 <- subset(hybrid3pcaDf,hybrid3pcaDf$neuron_name%in%rownames(ctsubset2))
colnames(hybrid3pcaDf)
unique(hybrid3pcaDf$br_groupttl)
dim(brcode)
dim(brcode2)
colnames(brcode)
ctsubset <- cbind(ctsubset,brgroupnum=brcode$groupnum,brgroupttl=brcode$groupttl)
colnames(ctsubset2)
unique(ctsubset2$groupttl)
unique(ctsubset2$brgroupttl)
tpos <- is.na(ctsubset2$brgroupttl)
length(tpos[tpos==T])
ctsubset2 <- cbind(ctsubset2,brgroupnum=brcode2$br_groupnum,brgroupttl=brcode2$br_groupttl)
dim(ctsubset2)
rownames(ctsubset)
colnames(ctsubset2)
unique(factor(ctsubset2$brgroupttl))
aggregate(ctsubset$brgroupttl,by=list(ctsubset$brgroupttl, ctsubset$brgroupnum),length)
aggregate(ctsubset2$brgroupttl,by=list(ctsubset2$brgroupttl, ctsubset2$brgroupnum),length)
aggregate(ctsubset2$groupttl,by=list(ctsubset2$groupttl, ctsubset2$groupnum),length)

#change the brain region labels
ctsubset$groupttl = as.character(ctsubset$groupttl)
ctsubset$brgroupttl = as.character(ctsubset$brgroupttl)
ctsubset <- editBrainRegionLabels(ctsubset,"brgroupttl")
#edit cell type labels
colnames(ctsubset2)
unique(ctsubset2$groupttl)
ctsubset2 <- editCelltypeLabels(ctsubset2,"groupttl")


#ctsubset <- addmocklabels(ctsubset,"brgroupttl")
#library(RColorBrewer)
#display.brewer.pal(7,"BrBG")
#display.brewer.pal(11,"Spectral")
#colpalette <- colorRampPalette(brewer.pal(11, "RdBu"))(diff(range(ctsubset$brgroupnum)))
#col=sample(colors()[-1], factor(ctsubset$brgroupnum), replace = FALSE) 
#factor(ctsubset$brgroupnum)
#require(grid)
#col=grid.raster(heat.colors(max(ctsubset$brgroupnum)), 
#                width=unit(1,"npc"),
#                height=unit(1,"npc"), int=FALSE),

#rgbPal <- colorRampPalette(c('red','green','blue'))
#set a manual color code

colcode <- unique(ctsubset$brgroupnum)
colcode
mypalette <- c("black", "red","green3","blue","cyan","magenta","yellow","gray","seagreen4","thistle4","violetred4","darkorange","darkmagenta","darkgreen")
length(mypalette)
plot(x=1:14,y=rep(1,14),pch=19,col=mypalette)
dim(ctsubset)
colarr <- c()
#loop and assign colors to each data point 
for(i in 1:14){
  colarr[ctsubset$brgroupnum==colcode[i]] <- i
}
unique(colarr)
plotcluster(ctsubset[,1:4],ctsubset$cluster,cex=0.8,
            col = mypalette[colarr],
            main="clusters and groups")
rjust <- legend("bottomright",text.width = strwidth("1,000,000"),
               legend = sapply(strsplit(as.character(unique(ctsubset$brgroupttl)),"Control/"),tail,1),
               pch=19,xjust = 0.5, yjust = 1,
               cex=0.7, col=unique(mypalette[colarr]))

colnames(ctsubset2)
unique(ctsubset2$brgroupnum)
unique(ctsubset2$brgroupttl)
unique(ctsubset2$groupnum)
unique(ctsubset2$groupttl)
dim(ctsubset2)
ctsubset2$brgroupttl <- as.character(ctsubset2$brgroupttl)
ctsubset2 <- editBrainRegionLabels(ctsubset2,"brgroupttl")
colarr2 <- clustercolorcode(ctsubset2,"brgroupnum")
plotcluster(ctsubset2[,1:4],ctsubset2$cluster,cex=0.8,
            col = mypalette[colarr2],
            main="clusters and groups")
rjust <- legend("bottomleft",text.width = strwidth("1,000,000"),
                legend = as.character(unique(ctsubset2$brgroupttl)),
                pch=19,xjust = 0.5, yjust = 1,
                cex=0.7, col=unique(mypalette[colarr2]))


#choose the same dataset as celltypes but plot species groups
ctsubset <- subset(f[,keepcols],f$groupnum%in%grplst)
dim(ctsubset)
colnames(ctsubset)
spcode <- subset(hybrid3pcaDf,hybridpcaDf$neuron_name%in%rownames(ctsubset))
dim(spcode)
NApos <- is.na(spcode$groupnum)
length(NApos[NApos==T])
spcode[NApos,"species"]
#spcode[NApos,"groupnum"] <- 99
#spcode[NApos,"groupttl"] <- "Control/amphibians/frog"
colnames(spcode)
unique(spcode$groupttl) 
ctsubset <- cbind(ctsubset,spgroupnum=spcode$groupnum,spgroupttl=spcode$groupttl)
dim(ctsubset)
rownames(ctsubset)
colnames(ctsubset)
unique(factor(ctsubset$spgroupttl))
aggregate(ctsubset$spgroupttl,by=list(ctsubset$spgroupttl, ctsubset$spgroupnum),length)

colcode <- unique(ctsubset$spgroupnum)
colcode
mypalette <- c("black", "red","green3","blue","cyan","magenta","yellow","gray","seagreen4","thistle4","violetred4")
length(mypalette)
plot(x=1:11,y=rep(1,11),pch=19,col=mypalette)
dim(ctsubset)
ctsubset <- subset(ctsubset, !is.na(ctsubset$spgroupnum))
colarr <- c()
#loop and assign colors to each data point 
for(i in 1:11){
  colarr[ctsubset$spgroupnum==colcode[i]] <- i
}
unique(colarr)

plotcluster(ctsubset[,1:4],ctsubset$cluster,cex=0.8,
            col = mypalette[colarr],
            main="clusters and groups")
rjust <- legend("bottomright",text.width = strwidth("1,000"),
                legend = sapply(strsplit(as.character(unique(ctsubset$spgroupttl)),"/"),tail,1),
                pch=19,xjust = 0.5, yjust = 1,
                cex=0.7, col=unique(mypalette[colarr]))
colnames(ctsubset2)
#bivariate density estimate plot
k <- kde2d(testDataNocor[,1],testDataNoCor[,2])
filled.contour(k)

#unique(factor(f$groupnum))
#col.txt=sapply(strsplit(as.character(unique(metalbls)),"/"),tail,1),
clusplot(f[,1:3],f$cluster,cex=0.8,color=TRUE,shade=TRUE, 
         labels=5,
         col.clus=unique(f$groupnum),
         lines=0,main="testing")
clusplot(pcamat$x, pcasub$groupnum, color=TRUE, shade=TRUE, labels=length(pcasub$groupnum),
         lines=0,main="testing")
clusplot(f, f$metagrps, color=length(unique(metagrps)),labels=4,
         lines=0,main="testing")
#sub sample of rat and mouse
pairs(~PC1+PC2+PC3+PC4,data=subset(pcawhole_matrix$x,rownames(pcawhole_matrix$x) %in% rownames(subset(psp,psp$groupnum%in%c(11,12)))),col=subset(psp,psp$groupnum%in%c(11,12))$groupnum,main="k-means scatterplot Matrix")
#sub sample of human and monkey
pairs(~PC1+PC2+PC3+PC4,data=subset(pcawhole_matrix$x,rownames(pcawhole_matrix$x) %in% rownames(subset(psp,psp$groupnum%in%c(9,8)))),col=subset(psp,psp$groupnum%in%c(9,8))$groupnum,main="k-means scatterplot Matrix")
#sub sample of human and mouse
pairs(~PC1+PC2+PC3+PC4,data=subset(pcawhole_matrix$x,rownames(pcawhole_matrix$x) %in% rownames(subset(psp,psp$groupnum%in%c(9,12)))),col=subset(psp,psp$groupnum%in%c(9,12))$groupnum,main="k-means scatterplot Matrix")
summary(pcamat)
str(pcamat)


#using ward's method on PCA matrix
#pcawhole_matrix <- prcomp(x=pcawhole[,2:33],center=T,scale=T)
#str(pcawhole_matrix$x)
plot(pcamat$x)
distxy <- dist(pcamat$x)
str(distxy)
hclustering <- hclust(distxy,method="ward")
plot(hclustering)
str(hclustering)
rect.hclust(hclustering,k=4, border="red")#display dendogram with red boarders around 4 clusters
nclus <- cutree(hclustering,k=4)
plot(pcamat$x,col=nclus)#scatter plot to see 6 clusters
for (i in 1:25) height[i] <- hclustering$height[i]
#height <- hclustering$height[sort(hclustering$order[1:100],decreasing=T)]
height <- sort(hclustering$height,decreasing=T)[1:15]

length(height)
plot(1:15,height,type="b",main="hierarchical agglomerative method",ylim=range(hclustering$height))#,xlab="Number of Cluster",ylab="height for clustering")
barplot()#a scree plot!?

#Determine number of clusters using ward's method
d <- dist(pcawhole[,2:33], method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=3) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=3, border="red")