#arbor type == 2
dim(metricDataAxon)
metricAxon <- metricDataAxon
colnames(metricAxon)
dim(metricAxon)
#remove NA rows
metricAxon <- cntNRmvNAs(metricAxon,113)
metricAxon$N_bifs_total_sum
###############Data cleaning######################
#check for zero depth data or flat 2D neurons
depthNA <- subset(metricAxon,is.na(metricAxon$Depth_total_sum))
dim(depthNA)
#depthNA$neuron_name
#remove CNG.swc from neuron_name
metricAxon <- clnData(metricAxon)

dim(metricAxon)
#a subset of metrics to be included in pca
#dendpcaMetrics <- pcaMetrics[-2]
dendpcaMetrics
#merge metric data with meta data
metricAxon <- getPlotSubset(metricAxon,dendpcaMetrics)
colnames(metricAxon)
#remove extra column that may have been inserted from merge process
metricAxon <- metricAxon[, !(colnames(metricAxon) %in% c("X"))]
#remove neurons with < 4 bifurcations
Nbranchoutlierpos <- with(metricAxon, which(metricAxon$N_bifs_total_sum < 4))
Nbranchoutlier <- metricAxon[Nbranchoutlierpos,c("neuron_name","N_bifs_total_sum","archive_name")]
str(Nbranchoutlier)
metricAxon <- metricAxon[-Nbranchoutlierpos,]
#remove Larkman archive because his neurons are dendrograms and hence aren't valid biological metrics.
larkmanpos <- with(metricAxon, which(metricAxon$archive_name=="Larkman"))
larkmanoutlier <- metricAxon[larkmanpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(larkmanoutlier)
#metricAxon <- metricAxon[-larkmanpos,]
#remove neurons with fractal==0
zerofractalpos <- with(metricAxon, which(metricAxon$Fractal_Dim_avg == 0))
zerofractaloutlier <- metricAxon[zerofractalpos,c("neuron_name","Fractal_Dim_avg","archive_name")]
str(zerofractaloutlier)
#metricAxon <- metricAxon[-zerofractalpos,]
hist(metricAxon$Helix_max, breaks=30, xlab="Helix")
hist(metricAxon$Fractal_Dim_avg, breaks=30, xlab="fractal")

min(metricAxon$Helix_max)
helix_zeropos <- with(metricAxon, which(metricAxon$Helix_max == 0.00))
helixzero <- metricAxon[helix_zeropos,c("neuron_name","Fractal_Dim_avg","Contraction_avg","Branch_pathlength_max","archive_name")]
dim(helixzero)
#log transform necessary columns
plotDf <- plotTesting(metricAxon)
colnames(plotDf)
dim(plotDf)
#standardize metric dataframe before computing PCA
for(i in 33:63){
  #get non-NA values for normalization
  colv <- plotDf[,i][!is.na(plotDf[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  plotDf[,i][!is.na(plotDf[,i])] <- colv/msd
  print(paste(names(plotDf)[i],msd))
}

###########add cell types groups##############
aggregate(plotDf$expercond, by = list(plotDf$expercond), length)
#celltype hierarchy
primary <- list(expercond="Control")
x <- getSubset2(plotDf, primary, typeMetaData=usedMetaData)
level1 <- "cellclass1"
level2 <- "cellclass2"
level3 <- "cellclass3"
#make groups following hierarchical ordering of metadata
ctMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData, level1, level2, level3,minSize=100)
printGrps(ctMetricHierarchy)
#group neuron_name (primariy id) according the hierarchy lists
Grplst <- addgrps2MetricDF(plotDf,ctMetricHierarchy)
for(i in 1:length(Grplst)){
  print(Grplst[[i]]$GroupSpec)
}
#assigning a groupnum and groupttl to each neuron according to the hierarchy list
plotDf_ct <- metricDf_plotting(plotDf,Grplst)

###########add brain region groups############
primary <- list(expercond = "Control")

level1 <- "region1"
level2 <- "lobes"
level3 <- "region2"
level4 <- "region3"
#unique(usedMetaData$lobes)
unique(usedMetaData$lobes)
brMetricHierarchy <- makeHierarchyGroups(primary, plotDf, usedMetaData, level1, level2,level3,level4,minSize=300)

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

###########add species groups#################
primary <- list(expercond = "Control")
level1 <- "order"
level2 <- "species"
level3 <- "strain"
spstrainMetricall <- makeHierarchyGroups(primary, plotDf,usedMetaData,level1,level2,minSize=55)
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

###########pca and determine #clusters for apical######

source("C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/plotutilfunctions.R")
#add rownames for three dimensions
rownames(plotDf_ct) <- plotDf_ct$neuron_name
rownames(plotDf_br) <- plotDf_br$neuron_name
rownames(plotDf_sp) <- plotDf_sp$neuron_name
dim(plotDf_ct)
dim(plotDf_br)
dim(plotDf_sp)
colnames(plotDf_ct)

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
#groups from 4D space combining archives/species/celltypes/brain regions
mixedgrpsDf <- subset(fourdims,fourdims$x>=55)
dim(mixedgrpsDf)
mixedgrpsDf
#prepare data matrix for computing pca and for cluster analysis
Axondatamatrix <- addmixedgrpcol(mixedgrpsDf, pcametrics)
dim(Axondatamatrix)
mxdgrpctlg <- aggregate(Axondatamatrix$neuron_name, as.list(Axondatamatrix[,c("archive_name","ct_groupttl","br_groupttl","mxdgroupnum")]), FUN = length)
write.table(mxdgrpctlg, file="mxdgrpTable4Axon.txt",quote=FALSE, sep=",")
colnames(Axondatamatrix)

###########perform kmeans clustering on PCs##########
pcamat <- prcomp(x=pcametrics[,33:63],center=T,scale=T)
pcasumm <- summary(pcamat)
#make scree plot to choose most contributing PCs
plot(pcamat$sdev^2/sum(pcamat$sdev^2)*100, type = "b", ylim = c(0,100), main = "Scree plot for dendrites (basal only) PCs",xlab="PCA Eigenvalues",ylab="% of variance")
#plot cumulative variance
lines(pcasumm$importance[3,]*100,main = "whole arbor metrics",type = "b", col = 2, xlab="cumulated PCA eigen values",ylab="% of variance")
legend(10,50,c("cumulative variance","variance"),lty=c(1,1),col=c("red","black"))

#six components were covering ~75% of variance
#choose eigenvalues >= 1
pcamat$sdev ^ 2
rownames(pcametrics)
colnames(pcametrics)

#test k=5 & 8 clusters based on the normalized scatter (atleast 2% drop)
kmeansObj <- kmeans(pcamat$x,centers=5,nstart=20,iter.max=30)
summary(kmeansObj)
str(kmeansObj)
dim(pcamat$x)

#pos <- psp$groupttl%in%c("Control/Carnivora/Cat","Control/Amphibians/Salamander")
#plot all celltypes groups 
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=kmeansObj$cluster,main="k-means scatterplot Matrix")
unique(pcametrics$ct_groupnum)
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=pcametrics$ct_groupnum,main="celltype groups scatterplot Matrix")
unique(pcametrics$br_groupnum)
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=pcametrics$br_groupnum,main="brain region groups scatterplot Matrix")
unique(pcametrics$sp_groupnum)
pairs(~PC1+PC2+PC3+PC4,data=pcamat$x,col=pcametrics$sp_groupnum,main="species groups scatterplot Matrix")
#sub sample of pyramidal cell others and interneuron others
pairs(~PC1+PC2+PC3+PC4,data=subset(pcamat$x,rownames(pcamat$x) %in% rownames(subset(nonpcaDf,nonpcaDf$groupnum%in%c(8,9,10)))),col=subset(nonpcaDf,nonpcaDf$groupnum%in%c(8,9,10))$groupnum,main="k-means scatterplot Matrix")
axonReduced <- data.frame(pcamat$x, cluster=kmeansObj$cluster, 
                           ct_groupnum=pcametrics$ct_groupnum, 
                           ct_groupttl=pcametrics$ct_groupttl,
                           br_groupnum=pcametrics$br_groupnum,
                           br_groupttl=pcametrics$br_groupttl,
                           sp_groupnum=pcametrics$sp_groupnum,
                           sp_groupttl=pcametrics$sp_groupttl)
colnames(axonReduced)
rownames(axonReduced)
unique(pcametrics$sp_groupttl)
unique(pcametrics$br_groupttl)
ps <- is.na(pcametrics$sp_groupttl)
length(ps[ps==TRUE])
unique(pcametrics$ct_groupttl)
#list combination of cell types and brain regions

#edit groupttl for all 3 categories
dim(pcametrics)
pcametrics <- subset(pcametrics, !is.na(pcametrics$sp_groupnum))# && !is.na(pcametrics$ct_groupnum) && !is.na(pcametrics$sp_groupnum)
pcametrics <- subset(pcametrics, !is.na(pcametrics$br_groupnum))# && !is.na(pcametrics$ct_groupnum) && !is.na(pcametrics$sp_groupnum)
pcametrics <- editBrainRegionLabels(pcametrics,"br_groupttl")
pcametrics <- editCelltypeLabels(pcametrics,"ct_groupttl")
pcametrics <- editspeciesLabels(pcametrics,"sp_groupttl")
fourdims <- aggregate(pcametrics$neuron_name, as.list(pcametrics[,c("archive_name","ct_groupttl","br_groupttl","sp_groupttl")]), FUN = length)
fourdims
nrow(fourdims)
dim(fourdims)
min(fourdims$x)
max(fourdims$x)
#groups from 4D space combining archives/species/celltypes/brain regions
mixedgrpsDf <- subset(fourdims,fourdims$x>=55)
dim(mixedgrpsDf)
mixedgrpsDf

ff <- addmixedgrpcol(mixedgrpsDf, pcametrics, 
                     subset(axonReduced,rownames(axonReduced)%in%pcametrics$neuron_name))
dim(ff)
colnames(fourdims)
colnames(ff)

unique(ff$mxdgroupnum)
unique(ff$ct_groupttl)
unique(ff$br_groupttl)
unique(ff[,c("sp_groupnum","sp_groupttl")])
unique(ff$sp_groupnum)
unique(ff$archive_name)

library(cluster)
library(fpc)
dim(ff)
#primate basal dendrites

grpnames <- aggregate(ff$mxdgroupnum,as.list(ff[,c("mxdgroupnum","archive_name","ct_groupttl","br_groupttl","sp_groupttl")]),length)
grpperc <- aggregate(ff$cluster,as.list(ff[,c("mxdgroupnum","cluster")]),length)
grpnames
grpperc
#grpperc <- data.frame(grpperc,"perc")
colnames(grpperc)
perc <- list()
for(i in 1:nrow(grpperc)){
  for(j in 1:nrow(grpnames)){
    if(grpperc[i,'mxdgroupnum']==grpnames[j,'mxdgroupnum']){
      print(grpnames[j,])
      print(round(grpperc[i,'x']*100/grpnames[j,'x'],2))
      perc<- append(perc, list(round(grpperc[i,'x']*100/grpnames[j,'x'],2)))
    }
  }
}
perc <- unlist(perc)
length(perc)
grpperc <- data.frame(grpperc, perc)
grpperc

mypalette <- c("black", "darkgrey", "red", "darkred","green","darkgreen","blue","darkorchid","cyan","darkcyan","magenta","darkmagenta","yellow","seagreen4","violetred4","darkorange")
length(mypalette)
plot(x=1:16,y=rep(1,16),pch=19,col=mypalette)
colcode <- unique(ff$mxdgroupnum)
colarr <- c()
#loop and assign colors to each data point 
for(i in 1:12){
  colarr[ff$mxdgroupnum==colcode[i]] <- i
}
unique(colarr)
#first 8 groups
pairs(~PC1+PC2+PC3+PC4,data=ff[,1:31],col=mypalette[colarr],main="mixed groups scatterplot Matrix")

plotcluster(ff[,1:4],ff$cluster,cex=0.8,col=mypalette[colarr],main="axons clustering and groups")
legend("bottomleft", legend = unique(ff$mxdgroupnum), pch=19, cex=0.7, col=unique(mypalette[colarr]))
aggregate(ff$cluster,as.list(ff[,c("mxdgroupnum","cluster")]),length)

#plot interneuron subgroup
interaxn <- subset(pcametrics, pcametrics$cellclass1=="Interneuron")
colnames(interaxn)
fourdims <- aggregate(interaxn$neuron_name, as.list(interaxn[,c("archive_name","ct_groupttl","br_groupttl","sp_groupttl")]), FUN = length)
fourdims
submixedgrpsDf <- subset(fourdims,fourdims$x>=55)
submixedgrpsDf

ff <- addmixedgrpcol(submixedgrpsDf, pcametrics, 
                     subset(axonReduced,rownames(axonReduced)%in%pcametrics$neuron_name))
dim(ff)
colnames(fourdims)
colnames(ff)

grpnames <- aggregate(ff$mxdgroupnum,as.list(ff[,c("mxdgroupnum","archive_name","ct_groupttl","br_groupttl","sp_groupttl")]),length)
grpperc <- aggregate(ff$cluster,as.list(ff[,c("mxdgroupnum","cluster")]),length)
grpnames
grpperc
#grpperc <- data.frame(grpperc,"perc")
colnames(grpperc)
perc <- list()
for(i in 1:nrow(grpperc)){
  for(j in 1:nrow(grpnames)){
    if(grpperc[i,'mxdgroupnum']==grpnames[j,'mxdgroupnum']){
      print(grpnames[j,])
      print(round(grpperc[i,'x']*100/grpnames[j,'x'],2))
      perc<- append(perc, list(round(grpperc[i,'x']*100/grpnames[j,'x'],2)))
    }
  }
}
perc <- unlist(perc)
length(perc)
grpperc <- data.frame(grpperc, perc)
grpperc

plotcluster(ff[,1:4],ff$cluster,cex=0.8,col=ff$mxdgroupnum,main="Interneuron axons clustering and groups")
legend("bottomleft", legend = unique(ff$mxdgroupnum), pch=19, cex=0.7, col=unique(ff$mxdgroupnum))
