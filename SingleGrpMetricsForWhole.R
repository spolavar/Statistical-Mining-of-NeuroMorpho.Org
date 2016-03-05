###step1: clusters on PC with over represented groups without the non-significant data###

sigdf_OG <- data.frame()
for(i in 1:nrow(sigMatrix)){
  sigdf_OG <- rbind(sigdf_OG, subset(wholeReduced, 
                               wholeReduced$mxdgroupnum == sigMatrix$row[i] & 
                                 wholeReduced[,clusfeature] == sigMatrix$col[i]))
}
dim(sigdf_OG)

#get the metric data for only significant groups
dim(datamatrix)
colnames(datamatrix)
sigmetricdf_OG <- subset(datamatrix, datamatrix$neuron_name %in% sigdf_OG$neuron_name)
dim(sigmetricdf_OG)

colnames(sigmetricdf_OG)
rownames(sigmetricdf_OG)
rownames(sigdf_OG)

###step2: analyze the clusters as determined by the EM###
sigdf <- wholeReduced
rownames(sigdf) <- rownames(wholeReduced)

#add pmid column to sigdf 
pmidfile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/pmid.txt"
pmidcol <- read.csv(pmidfile,header = TRUE, sep="\t",na.strings=c("NA",""),quote="\"",stringsAsFactors=FALSE)
dim(pmidcol)
colnames(pmidcol)
sigdf <- merge(x=sigdf,y=pmidcol,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(sigdf)
colnames(sigdf)

sigdenddf <- merge(x=sigdenddf,y=pmidcol,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(sigdenddf)
colnames(sigdenddf)

###step3: clusters on metric data in comparison with PC data.###
sigmetricdf <- datamatrix
dim(sigmetricdf)
colnames(sigmetricdf)
#sigmetricdf <- sigmetricdf[order(rownames(sigmetricdf)),]
#sigdf <- sigdf[order(rownames(sigdf)),]
#sigmetricdf$classification <- sigdf$classification

#add pmid col to sigmetricdf
sigmetricdf <- merge(x=sigmetricdf,y=pmidcol,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(sigmetricdf)
colnames(sigmetricdf)
#flip the -ve loadings on PC1 (size)
sigdf$PC1 <- -1 * sigdf$PC1
sigdf_OG$PC1 <- -1 * sigdf_OG$PC1

###step4: separate PC and metric parameters by clusters###
cluster4 <- subset(sigdf, sigdf$classification == 4)
dim(cluster4)
cluster5 <- subset(sigdf, sigdf$classification == 5)
dim(cluster5)
cluster6 <- subset(sigdf, sigdf$classification == 6)
dim(cluster6)
cluster7 <- subset(sigdf, sigdf$classification == 7)
dim(cluster7)
cluster2<- subset(sigdf, sigdf$classification == 2)
dim(cluster2)
cluster1 <- subset(sigdf, sigdf$classification == 1)
dim(cluster1)
cluster3 <- subset(sigdf, sigdf$classification == 3)
dim(cluster3)

clusterM7 <- subset(sigmetricdf,sigmetricdf$classification == 7)
#colnames(clusterM7)
clusterM4 <- subset(sigmetricdf,sigmetricdf$classification == 4)
#colnames(clusterM4)
clusterM6 <- subset(sigmetricdf,sigmetricdf$classification == 6)
#colnames(clusterM6)
clusterM3 <- subset(sigmetricdf,sigmetricdf$classification == 3)
#colnames(clusterM3)
unique(clusterM3$mxdgroupnum)
clusterM5 <- subset(sigmetricdf,sigmetricdf$classification == 5)
#colnames(clusterM5)
clusterM2 <- subset(sigmetricdf,sigmetricdf$classification == 2)
dim(clusterM2)
clusterM1 <- subset(sigmetricdf,sigmetricdf$classification == 1)
dim(clusterM1)

#density plot of all clusters along PC1
plotDensitiesOnPCsForAllClusters(sigdf, 'PC1')


#check for features that are unique in clusterM1 through PCA
cluster1EM <- Mclust(cluster1[,33:50],G=1:15)
summary(cluster1EM)
plot(cluster1EM, what="classification",dimens=1:2)
plot(cluster1EM, what="classification",dimens=5:7)

N_cluster1 <- reEMonSubsets(cluster1)

dim(N_cluster1)
colnames(N_cluster1)
unique(N_cluster1$classification)
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdf)
sigdf[which(sigdf$classification %in% c(3) & sigdf$PC2 > 4), c('neuron_name','archive_name','mxdgroupnum','cellclass2','cellclass3','pmid')]
PlotDoubleMetricDistrOnClusters('PC1','PC2',N_cluster1)
#S1 pyramidal cells with large axons
N_cluster1[which(N_cluster1$classification==9),c("neuron_name","archive_name","species","cellclass3","region3")]
#Gulyas CA1 interneurons
N_cluster1[which(N_cluster1$classification==8),c("neuron_name","archive_name","species","cellclass3","region3")]
#interneurons from perirhinal cortex with a few others from other non-cortical regions
N_cluster1[which(N_cluster1$classification==4),c("neuron_name","archive_name","species","cellclass2","region2")]
#neurogliaform cells
N_cluster1[which(N_cluster1$classification==7),c("neuron_name","archive_name","species","cellclass2","region2")]
N_cluster1[which(N_cluster1$classification==5),c("neuron_name","archive_name","species","cellclass2","region2")]
cluster1[which(cluster1$mxdgroupnum==35),c("neuron_name","archive_name","species","cellclass2","region2")]


#Brown archive data
cluster1[which(cluster1$mxdgroupnum==31),c("neuron_name","age_scale","age_class","species","cellclass2","region2")]
#gulyas archive data
cluster1[which(cluster1$mxdgroupnum==23),c("neuron_name","age_scale","age_class","species","cellclass2","region2")]

PlotDoubleMetricDistrOnClusters('Bif_ampl_local_avg','Bif_tilt_local_avg',clusterM1)
PlotDoubleMetricDistrOnClusters('Bif_ampl_remote_avg','Bif_tilt_local_avg',clusterM1)
PlotDoubleMetricDistrOnClusters('Bif_ampl_local_avg','Bif_tilt_local_avg',clusterM4)
PlotDoubleMetricDistrOnClusters('Bif_ampl_remote_avg','Bif_tilt_local_avg',clusterM4)

N_cluster7 <- reEMonSubsets(cluster7)
PlotDoubleMetricDistrOnClusters('PC1','PC2',N_cluster7)
#human magnopyramidal cells
N_cluster7[which(N_cluster7$classification==1),c("neuron_name","archive_name","species","cellclass3","region3","pmid")]
#monkey pyramidal
N_cluster7[which(N_cluster7$classification==2),c("neuron_name","archive_name","species","cellclass3","region3")]
#human basal pyramidal
N_cluster7[which(N_cluster7$classification==3),c("neuron_name","archive_name","species","cellclass3","region3")]

dim(cluster5)
N_cluster5 <- reEMonSubsets(cluster5)
PlotDoubleMetricDistrOnClusters('PC1','PC2',N_cluster5)
colnames(N_cluster5)
unique(N_cluster5[,c('cellclass2','mxdgroupnum','pmid')])
aggregate(N_cluster5$neuron_name, by=list(cluster=N_cluster5$classification,group=N_cluster5$mxdgroupnum,
                              cellclass=N_cluster5$cellclass2,archive=N_cluster5$archive_name),FUN=length)
#basket cells and tangential cells
N_cluster5[which(N_cluster5$classification %in% c(8,9,12,11)),c('neuron_name','classification','age_class','age_scale')]
#layer 2/3 interneurons
N_cluster5[which(N_cluster5$classification %in% c(8,10)),c('neuron_name','classification','age_class','age_scale')]
#Notreported S1 interneurons
N_cluster5[which(N_cluster5$mxdgroupnum ==19),c('neuron_name','classification','age_class','age_scale','cellclass2')]
#S1 pyramidal cells
N_cluster5[which(N_cluster5$mxdgroupnum %in% c(34,16,14)),c('neuron_name','classification','age_class','age_scale','cellclass2')]

N_cluster5[which(N_cluster5$classification==7),c("neuron_name","archive_name","species","cellclass3","region3")]
N_cluster5[which(N_cluster5$classification==9),c("neuron_name","archive_name","species","cellclass3","region3")]
N_cluster5[which(N_cluster5$classification==12),c("neuron_name","archive_name","species","cellclass3","region3")]
unique(N_cluster5[which(N_cluster5$mxdgroupnum == 33),'classification'])

PlotDoubleMetricDistrOnGroups('Bif_ampl_local_avg','Bif_tilt_local_avg',subset(clusterM5,clusterM5$mxdgroupnum%in%c(33,32,3,24)))
PlotDoubleMetricDistrOnGroups('Bif_ampl_remote_avg','Bif_tilt_local_avg',subset(clusterM5,clusterM5$mxdgroupnum%in%c(33,32,3,24)))

N_cluster3 <- reEMonSubsets(cluster3)
PlotDoubleMetricDistrOnClusters('PC1','PC2',N_cluster3)
aggregate(N_cluster3$neuron_name, by=list(cluster=N_cluster3$classification,group=N_cluster3$mxdgroupnum,
                                          cellclass=N_cluster3$cellclass2,archive=N_cluster3$archive_name),FUN=length)
#yuste interneurons
N_cluster3[which(N_cluster3$mxdgroupnum ==19),c('neuron_name','classification','age_class','age_scale','cellclass2')]


N_cluster4 <- reEMonSubsets(cluster4)
PlotDoubleMetricDistrOnClusters('PC1','PC2',N_cluster4)

aggregate(N_cluster3$neuron_name, by=list(cluster=N_cluster3$classification,group=N_cluster3$mxdgroupnum,
                                          cellclass=N_cluster3$cellclass2,archive=N_cluster3$archive_name),FUN=length)


###step3: plot z-scores (mean/SD) of each component for each cluster
colnames(sigdf)
zscr <- computezscore(sigdf[,c(33:60,68)], clusfeature)
colnames(zscr)
EucDistFrmGlobalmean <- computezscoreDist(zscr)
#compute zscore using mean/clusterSD for each cluster
zscr_1 <- computezscore2(sigdf[,c(33:60,68)], clusfeature)
zscr_1[1,2]
colnames(sigdf[,c(33:60,68)])
EucDistFrmGlobalmean_1 <- computezscoreDist(zscr_1)

varscr <- computegrpvar(sigdf[,c(33:60,68)], clusfeature)
dim(varscr)
varscr[which(zscr$c==7),'sd']
varscr[which(zscr$c==7),'mu']
varscr[which(zscr$c==7),'meanbysd']
varianceMagnitude <- computevarscrMag(varscr)
varzscr <- merge(EucDistFrmGlobalmean_1, varianceMagnitude, by=0, all=TRUE)
varzscr <- varzscr[order(varzscr$EucdistanceFromMu,decreasing=T),]
#write.table(EucDistFrmGlobalmean, file="clusterzscoreDist.txt",quote=FALSE, sep=",")
#str(EucDistFrmGlobalmean)
#zscr[which(zscr$c==7),'m']
#barplot with x-axis labels
#mp <- barplot(EucDistFrmGlobalmean$EucdistanceFromMu,main="Euclidean distance from clusters to global mean",cex.main=1,
#        xlab = "clusters")
#axis(1,at=mp,labels=rownames(EucDistFrmGlobalmean))

#barplot with x-axis labels
plotTop <- max(varzscr$EucdistanceFromMu)
mp <- barplot(varzscr$EucdistanceFromMu,main="z-score distance from clusters to global mean",cex.main=1,
              xlab = "clusters",ylim=c(0,plotTop))
axis(1,at=mp,labels=rownames(varzscr))
#arrows(mp, varzscr$EucdistanceFromMu-varzscr$magnitude, 
#       mp, varzscr$EucdistanceFromMu+varzscr$magnitude,
#       lwd=2, angle=90, code=3)
plotTop <- max(varzscr$magnitude)
mp <- barplot(varzscr$magnitude,main="SDs of each cluster",cex.main=1,
              xlab = "clusters",ylim=c(0,plotTop))
axis(1,at=mp,labels=rownames(varzscr))

#compute
DavidBouldinDist <- interclusterdist_DB(sigdf[,c(33:60,68)], clusfeature)
plotdistances(DavidBouldinDist)
write.table(DavidBouldinDist, file="DavidBouldinDist.txt",quote=FALSE, sep=",")

#check cluster 4 'soma_surface' values 
dim(clusterM4)
colnames(clusterM4)
aggregate(clusterM4[,c(3,7, 11,12,33)], by=list(clusterM4$mxdgroupnum), FUN=min)
aggregate(clusterM4[,c(3,33)], by=list(clusterM4$mxdgroupnum), FUN=max)
DavidBouldinDist <- interclusterdist_DB(sigdf[,c(33:64,72,73)], clusfeature)
write.table(DavidBouldinDist, file="DavidBouldinDist.txt",quote=FALSE, sep=",")

#input Dmatrix data frame
plotdistances <- function(Dmatrix){
  grpVals <- row.names(Dmatrix)
  topY <- max(Dmatrix, na.rm=TRUE)
  s <- Dmatrix[which(row.names(Dmatrix)==grpVals[1]),] # || names(Dmatrix)==grpVals[1]),]
  #s[which(is.na(s))] <- 0
  str(s)
  #names(s[1,])
  # specify the data
  x <- c(1:length(grpVals)) #; y <- x; z <- 10/x
  plot(as.numeric(s[1,]),type="b",ylim = c(0,topY),xlab="clusters",ylab = "pairwise distances", xaxt = 'n')
  #ylab="pairwise distances using std. mean/sd",main="comparing pairwise distances",add=TRUE)
  axis(1,at=x,labels=paste(colnames(s[1,])))
#loop for plotting other clusters
  for(i in 2:length(grpVals)){
    s <- Dmatrix[which(row.names(Dmatrix)==grpVals[i]),] # || names(Dmatrix)==grpVals[i]),]
    #s[which(is.na(s))] <- 0
    print(grpVals[i])
    print(as.numeric(s[1,]))
    lines(as.numeric(s[1,]),type="b",lty=i,pch=i)
  }
  #legend("topright",text.width = strwidth("1"),legend = grpVals, pch=c(1:length(grpVals)),cex=0.7)
}
 

#Davies–Bouldin index algo.
interclusterdist_DB <- function(coordsdf,grpTyp){
  grpVals <- unique(coordsdf[,grpTyp])
  dfsz <- length(coordsdf)-1
  print(dfsz)
  m <- matrix(nrow=length(grpVals),ncol=length(grpVals))
  DBmatrix <- as.data.frame(m)
  row.names(DBmatrix) <- grpVals
  names(DBmatrix) <- grpVals
  for(i in 1:length(grpVals)){
    #means of each cluster
    print(paste("reading for", grpVals[i]))
    x <- subset(coordsdf[,1:dfsz],coordsdf[,grpTyp] == grpVals[i])
    print(dim(x))
    clusmeanX <- sapply(x, mean)
    clusscatterX <- avgscatter(x)
    for(j in 1:length(grpVals)){
      #make diagonal values zeros
      if(i == j)
        DBmatrix[i,j] <- 0
      #populate the upper triangle
      else if(i < j){
        print(paste("to ", grpVals[j]))
        y <- subset(coordsdf[,1:dfsz],coordsdf[,grpTyp] == grpVals[j])
        print(dim(y))
        clusmeanY <- sapply(y, mean)
        clusscatterY <- avgscatter(y)
        meanD <- pairwiseMDdist(clusmeanX, clusmeanY) 
        print(meanD)
        DBmatrix[i,j] <- meanD/(clusscatterX + clusscatterY)
        DBmatrix[j,i] <- DBmatrix[i,j]
        print( DBmatrix[j,i])
      }
    }
  }
  return (DBmatrix)
}

#compute the average scatter of each point from its centroid
avgscatter <- function(clusdf){
  c <- sapply(clusdf, mean)
  S <- 0
  sctr <- -1
  if(length(clusdf[i,] == length(c))){
    for(i in 1:nrow(clusdf)){
      S <- S + pairwiseMDdist(clusdf[i,],c)
    }   
    #calculate avg scatter of each point in the cluster
    sctr <- S/nrow(clusdf) 
  }
  return (sctr)
}

#coordsdf: PC or metric data along with group information
#grpTyp: either has 'classification' or 'mxdgroupnum' as grouping criterion
#grpVals: array of group values according to grpTyp e.g. c(23, 31, 6) or c(1,2,3,4,5,6,7) 
computezscore <- function(coordsdf, grpTyp){
  grpVals <- unique(coordsdf[,grpTyp])
  #SD and mean of the entire data
  dfsz <- length(coordsdf)-1
  print(dfsz)
  globalSD <- sapply(coordsdf[,1:dfsz], sd)
  globalMean <- sapply(coordsdf[,1:dfsz], mean)
  print(dim(globalMean))
  #a single array containing z-scores of all groups
  zscr <- data.frame(c=character(0), m=I(vector('list', 0)),stringsAsFactors=FALSE)
    
  #print(zscr)
  for(i in 1:length(grpVals)){
    print("**************")
    #print(length(grpVals))
    x <- subset(coordsdf,coordsdf[,grpTyp] == grpVals[i])
    dim(x)
    clusterMean <- aggregate(x[,1:(length(coordsdf)-1)],by=list(x[,grpTyp]),FUN=mean)
    clusterMean <- as.matrix(clusterMean[,2:length(clusterMean)])

    S <- (clusterMean-globalMean)/globalSD
    print("Z-score")
    zscrpercluster <- matrix(0, nrow=1,ncol=dfsz)
    for(j in 1:length(S)){
      zscrpercluster <- zscrpercluster +S^2
    }
    zscrpercluster <- sqrt(zscrpercluster)
    print(zscrpercluster)
    zscr[i,'c'] <- as.character(grpVals[i])
    zscr[[i,'m']] <- list(zscrpercluster)
  }
  #colnames(zscr) <- "mean/globalSD"
  return (zscr)
}

#coordsdf: PC or metric data along with group information
#grpTyp: either has 'classification' or 'mxdgroupnum' as grouping criterion
#compute z-score by normalizing the distance by SD of EACH cluster
computezscore2 <- function(coordsdf, grpTyp){
  grpVals <- unique(coordsdf[,grpTyp])
  print(grpVals)
  dfsz <- length(coordsdf)-1
  print(dfsz)
  #SD and mean of the entire data
  #globalSD <- sapply(coordsdf[,1:32], sd)
  globalMean <- sapply(coordsdf[,1:dfsz], mean)
  print(length(globalMean))
  #a single array containing z-scores of all groups
  zscr <- data.frame(c=character(0), m=I(vector('list', 0)),stringsAsFactors=FALSE)
  
  #print(zscr)
  for(i in 1:length(grpVals)){
    print("**************")
    #print(length(grpVals))
    x <- subset(coordsdf,coordsdf[,grpTyp] == grpVals[i])
    print(dim(x))
    clusterMean <- aggregate(x[,1:dfsz],by=list(x[,grpTyp]),FUN=mean)
    print(length(clusterMean))
    clusterMean <- as.matrix(clusterMean[,2:length(clusterMean)])
    clusterSD <- sapply(x[,1:dfsz], sd)
    S <- (clusterMean-globalMean)/clusterSD
    print(dim(S))
    print("Z-score")
    zscrpercluster <- matrix(0, nrow=1,ncol=dfsz)
    
    for(j in 1:length(S)){
      zscrpercluster <- zscrpercluster +S^2
    }
    zscrpercluster <- sqrt(zscrpercluster)
    print(zscrpercluster)
    zscr[i,'c'] <- as.character(grpVals[i])
    zscr[[i,'m']] <- list(zscrpercluster)
  }
  #colnames(zscr) <- "mean/globalSD"
  return (zscr)
}

#z-score function that calculates the z-score of a group from the centroid of the entire data#
#plot the normalized Z-scores
ymax <- max(unlist(zscr_1[[which(zscr_1$c==7),2]])[1:18])
plot(unlist(zscr_1[[which(zscr_1$c==1),2]])[1:18],type="b",ylim = c(0,ymax),xlim=c(1,18),xlab="principal components",ylab="distances from the data centroid",main="z-scores of each cluster on whole data")
lines(unlist(zscr_1[[which(zscr_1$c==2),2]])[1:18],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(zscr_1[[which(zscr_1$c==3),2]])[1:18],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(zscr_1[[which(zscr_1$c==4),2]])[1:18],type="b",lty=4,pch=12,col="red")
lines(unlist(zscr_1[[which(zscr_1$c==5),2]])[1:18],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(zscr_1[[which(zscr_1$c==6),2]])[1:18],type="b",lty=6,pch=13)#,col="lightgreen")
lines(unlist(zscr_1[[which(zscr_1$c==7),2]])[1:18],type="b",lty=7,pch=2)#,col="red")
legend("topright",text.width = strwidth("1"),legend = sort(zscr_1$c), pch=c(1,0,6,12,7,13,2), cex=0.7, col=c(rep(1,3),2,rep(1,3)))
#legend("topright",text.width = strwidth("1"),legend = c(1:params$variance$G), pch=c(1,0,6,12,7,13,2), cex=0.7, col=c(rep(1,3),2,rep(1,3)))

#compute single magnitude of the variance of each group
computevarscrMag <- function(v){
  length(v)
  varmag <- data.frame(magnitude = numeric(nrow(v)),row.names=v$c)
  for(i in 1:nrow(v)){
    dist <- 0
    for(j in 1:length(v$sd[[i]])){
      dist <- dist + unlist(v$sd[[i]])[j]^2
    }
    dist <- sqrt(dist)
    #print(dist)
    varmag[i,"magnitude"] <- dist
  }
  #rownames(EucDistFrmGlobalmean) <- zscrV$c
  
  #EucDistFrmGlobalmean <- sort(EucDistFrmGlobalmean,decreasing=TRUE)
  varmag <- varmag[order(varmag[,1],decreasing=TRUE),,drop=FALSE]
  print(varmag)
  return (varmag)
}

#compute distances from cluster z-scores to globalmean
computezscoreDist <- function(zscrV){
  length(zscrV)
  EucDistFrmGlobalmean <- data.frame(EucdistanceFromMu = numeric(nrow(zscrV)),row.names=zscrV$c)
  for(i in 1:nrow(zscrV)){
    dist <- 0
    for(j in 1:length(zscrV$m[[i]])){
      dist <- dist + unlist(zscrV$m[[i]])[j]^2
    }
    dist <- sqrt(dist)
    #print(dist)
    EucDistFrmGlobalmean[i,"EucdistanceFromMu"] <- dist
  }
  #rownames(EucDistFrmGlobalmean) <- zscrV$c
 
  #EucDistFrmGlobalmean <- sort(EucDistFrmGlobalmean,decreasing=TRUE)
  EucDistFrmGlobalmean <- EucDistFrmGlobalmean[order(EucDistFrmGlobalmean[,1],decreasing=TRUE),,drop=FALSE]
  print(EucDistFrmGlobalmean)
  return (EucDistFrmGlobalmean)
}


###Step4: plot the variance of each cluster along each component###
computegrpvar <- function(coordsdf, grpTyp){
  grpVals <- unique(coordsdf[,grpTyp])
  dfsz <- length(coordsdf)-1
  print(dfsz)
  varscr <- data.frame(c=character(0), sd=I(vector('list', 0)),mu= I(vector('list',0)),meanbysd=I(vector('list',0)),stringsAsFactors=FALSE)
  #stdMeanMatrix <- matrix(data=NA, nrow= length(grpVals),ncol=length(grpVals))
  #rownames(varscr) <- c(1:7)
  #colnames(varscr) <- "mean/SD" 
  #colnames(subset(sigdf[,c(33:64,72)],sigdf$classification == 1))
  for(i in 1:length(grpVals)){
    x <- subset(coordsdf,coordsdf[,grpTyp] == grpVals[i])
    dim(x)
    clusterMean <- aggregate(x[,1:dfsz],by=list(x[,grpTyp]),FUN=mean)
    clusterMean <- as.matrix(clusterMean[,2:length(clusterMean)])
    clusSD <- sapply(x[,1:dfsz], sd)
    #compute the individual clusters mean/SD 
    #stdMeanMatrix[i] <- list(clusterMean/clusSD)
    #compute the variance from ratio of SDs
    varscr[i,'c'] <- as.character(grpVals[i])
    varscr[[i,'sd']] <- list(clusSD)
    varscr[[i,'mu']] <- list(clusterMean)
    varscr[[i,'meanbysd']] <- list(clusterMean/clusSD)
    #varscr[[i,'sd']] <- list(clusSD/globalSD)
  }
  return (varscr)
}

#plot the normalized and standardized variances
base <- (unlist(zscr_1$m[[1]])+unlist(zscr_1$m[[2]])+unlist(zscr_1$m[[3]])+
         unlist(zscr_1$m[[4]])+unlist(zscr_1$m[[5]])+unlist(zscr_1$m[[6]])+unlist(zscr_1$m[[7]]))/7

for(i in 1:6){
  for(j in (i+1):7){
   
    ymax <- max(unlist(zscr_1$m[[i]]), unlist(zscr_1$m[[j]]))
    plot(base[1:18],type="b",cex=0.7, ylim=c(0,ymax),xlim=c(1,18),col="red",
         xlab="principal components",ylab="zscore",main=paste("zscores between cluster pairs:",varscr$c[i],varscr$c[j]))
    lines(unlist(zscr_1$m[[i]])[1:18],type="b",lty=2,pch=zscr_1$c[i])#,col="darkred")
    lines(unlist(zscr_1$m[[j]])[1:18],type="b",lty=3,pch=zscr_1$c[j])#,col="darkorange")
  }
}

#plot SDs of all clusters 
ymax <- max(unlist(varscr[[which(varscr$c==1),2]])[1:18])
plot(unlist(varscr$sd[[1]])[1:18],type="b",cex=0.7, ylim=c(0,ymax),xlim=c(1,18),xlab="principal components",ylab="SDs",main="variances of each cluster")
lines(unlist(varscr$sd[[2]])[1:18],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(varscr$sd[[3]])[1:18],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(varscr$sd[[4]])[1:18],type="b",lty=4,pch=12,col="red")
lines(unlist(varscr$sd[[5]])[1:18],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(varscr$sd[[6]])[1:18],type="b",lty=6,pch=13)#,col="lightgreen")
lines(unlist(varscr$sd[[7]])[1:18],type="b",lty=7,pch=2)#,col="red")
legend("topright",text.width = strwidth("1"),legend = c(1:params$variance$G), pch=c(1,0,6,12,7,13,2), cex=0.7, col=c(rep(1,3),2,rep(1,3)))

#plot means of all clusters
ymax <- max(unlist(varscr[[which(varscr$c==7),3]])[1:18])
ymin <- min(unlist(varscr[[which(varscr$c==5),3]])[1:18])
plot(unlist(varscr$mu[[1]])[1:18],type="b",cex=0.7, ylim=c(ymin,ymax),xlim=c(1,18),xlab="principal components",ylab="means",main="means of each cluster")
lines(unlist(varscr$mu[[2]])[1:18],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(varscr$mu[[3]])[1:18],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(varscr$mu[[4]])[1:18],type="b",lty=4,pch=12,col="red")
lines(unlist(varscr$mu[[5]])[1:18],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(varscr$mu[[6]])[1:18],type="b",lty=6,pch=13)#,col="lightgreen")
lines(unlist(varscr$mu[[7]])[1:18],type="b",lty=7,pch=2)#,col="red")
legend("topright",text.width = strwidth("1"),legend = c(1:params$variance$G), pch=c(1,0,6,12,7,13,2), cex=0.7, col=c(rep(1,3),2,rep(1,3)))


###step5: plot all clusters from the multidimensional space of the entire data###
#plot all clusters with just the over represented groups
library("car")
library("tourr")
featurecol <- sigdf$classification
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
plot(sigdf[,33],sigdf[,34],col="white",cex=.55,pch=16,xlab=names(sigdf)[33],
     main="4196 whole arbor data classified into 7 models",ylab=names(sigdf)[34],cex.main=0.8)
text(sigdf[,33],sigdf[,34],labels=as.character(sigdf$classification),col=mypalette[colarr],font=2,cex=.7)
legend("bottomleft",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))


clusnums <- unique(tmpEM$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[1:2,clusnums[i]],shape=parms$variance$sigma[1:2,1:2,clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}
sigdf[,c(33:35,72)]
animate(sigdf[,33:37], grand_tour(), display_xy(col=mypalette[colarr]))

#plot density
require(MASS)
oldPar <- par()
x11()
dens2D <- kde2d(sigdf[,33], sigdf[,34])
kde2dplot(dens2D)
par(oldPar)


###step7: plots highlighting the cell groups within each cluster###
colnames(sigdf)
unique(sigdf$ct_groupttl)
#create a data frame with groups
cellgrpsdf <- data.frame()
#create a subset of basal pyramidal cells and its comparative groups
cellgrpMeans <- data.frame(c=character(0), mu=I(vector('list', 0))
                           ,sd=I(vector('list', 0)),stringsAsFactors=FALSE)
#find centroid of basal pyramidal dendrites
baspyr <- subset(sigdf[,c(33:64,72)],sigdf$classification == 7)
dim(baspyr)
baspyr$cellgrp <- apply(baspyr,1,function(r) "basal pyramidal")
colnames(baspyr)
cellgrpsdf <- rbind(cellgrpsdf, baspyr)
cellgrpMeans[1,] <- add2cellgrpMeans(baspyr)
cellgrpMeans$mu
cellgrpMeans[[1,'mu']][,'PC4']

CA1_intr <-  subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum == 23)
CA1_intr$cellgrp <- apply(CA1_intr,1,function(r) "CA1 INs" )
dim(CA1_intr)
cellgrpsdf <- rbind(cellgrpsdf, CA1_intr)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(CA1_intr)

blwfly_intr <-  subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum == 3)
blwfly_intr$cellgrp <- apply(blwfly_intr,1,function(r) "sensory INs")
cellgrpsdf <- rbind(cellgrpsdf, blwfly_intr)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(blwfly_intr)

CA3_intr <- subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum == 31)
CA3_intr$cellgrp <- apply(CA3_intr,1,function(r) "CA3 INs (4)" )
cellgrpsdf <- rbind(cellgrpsdf, CA3_intr)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(CA3_intr)

SS_pyr_4 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum %in% c(26,36,20))
SS_pyr_4$cellgrp <- apply(SS_pyr_4,1,function(r) "SS pyramidal (4)")
cellgrpsdf <- rbind(cellgrpsdf, SS_pyr_4)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(SS_pyr_4)

othr_pyr_4 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum == 18)
othr_pyr_4$cellgrp <- apply(othr_pyr_4,1,function(r) "cortical pyramidal (4)")
cellgrpsdf <- rbind(cellgrpsdf, othr_pyr_4)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(othr_pyr_4)

mnky_pfc_4 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 4 & sigdf$mxdgroupnum == 12)
mnky_pfc_4$cellgrp <- apply(mnky_pfc_4,1,function(r) "monkey frontal pyramidal")
cellgrpsdf <- rbind(cellgrpsdf,mnky_pfc_4)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(mnky_pfc_4)

ctx_bskt_5 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 5 & sigdf$mxdgroupnum %in% c(24,33))
ctx_bskt_5$cellgrp <- apply(ctx_bskt_5,1,function(r) "basket cells(cortical)")
cellgrpsdf <- rbind(cellgrpsdf,ctx_bskt_5)
cellgrpMeans[nrow(cellgrpMeans)+1,] <- add2cellgrpMeans(ctx_bskt_5)

vis_pyr_5 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 5 & sigdf$mxdgroupnum == 14)
vis_pyr_5$cellgrp <- apply(vis_pyr_5,1,function(r) "V1 pyramidal")
cellgrpsdf <- rbind(cellgrpsdf,vis_pyr_5)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(vis_pyr_5)

ss_pyr_5 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 5 & sigdf$mxdgroupnum %in% c(16,20,26,36))
ss_pyr_5$cellgrp <- apply(ss_pyr_5,1,function(r) "SS pyramidal (5)")
cellgrpsdf <- rbind(cellgrpsdf,ss_pyr_5)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(ss_pyr_5)

ss_intr_5 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 5 & sigdf$mxdgroupnum %in% c(19,34))
ss_intr_5$cellgrp <- apply(ss_intr_5,1,function(r) "SS INs (5)")
cellgrpsdf <- rbind(cellgrpsdf,ss_intr_5)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(ss_intr_5)

othr_pyr_5 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 5 & sigdf$mxdgroupnum == 18)
othr_pyr_5$cellgrp <- apply(othr_pyr_5,1,function(r) "cortical pyramidal (5)")
cellgrpsdf <- rbind(cellgrpsdf,othr_pyr_5)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(othr_pyr_5)

frontal_pyr_2 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 2 & sigdf$mxdgroupnum %in% c(13,25))
frontal_pyr_2$cellgrp <- apply(frontal_pyr_2,1,function(r) "Frontal pyramidal")
cellgrpsdf <- rbind(cellgrpsdf,frontal_pyr_2)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(frontal_pyr_2)

MSN_2 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 2 & sigdf$mxdgroupnum == 17)
MSN_2$cellgrp <- apply(MSN_2,1,function(r) "non-cortical non-pyramidal (2)")
cellgrpsdf <- rbind(cellgrpsdf,MSN_2)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(MSN_2)

olf_pyr_2 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 2 & sigdf$mxdgroupnum == 28)
olf_pyr_2$cellgrp <- apply(olf_pyr_2,1,function(r) "OB pyramidal (2)")
cellgrpsdf <- rbind(cellgrpsdf,olf_pyr_2)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(olf_pyr_2)

olf_npyr_2 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 2 & sigdf$mxdgroupnum == 27)
olf_npyr_2$cellgrp <- apply(olf_npyr_2,1,function(r) "OB non-pyramidal (2)")
cellgrpsdf <- rbind(cellgrpsdf,olf_npyr_2)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(olf_npyr_2)

human_PFC_6 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 6 & sigdf$mxdgroupnum == 8)
human_PFC_6$cellgrp <- apply(human_PFC_6,1,function(r) "human PFC")
cellgrpsdf <- rbind(cellgrpsdf,human_PFC_6)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(human_PFC_6)

mnky_PFC_6 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 6 & sigdf$mxdgroupnum == 12)
mnky_PFC_6$cellgrp <- apply(mnky_PFC_6,1,function(r) "monkey PFC")
cellgrpsdf <- rbind(cellgrpsdf,mnky_PFC_6)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(mnky_PFC_6)

PR_intr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 32)
PR_intr_1$cellgrp <- apply(PR_intr_1,1,function(r) "PR interneurons")
cellgrpsdf <- rbind(cellgrpsdf,PR_intr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(PR_intr_1)

Nctx_npyr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum %in% c(29,30))
Nctx_npyr_1$cellgrp <- apply(Nctx_npyr_1,1,function(r) "non-cortical non-pyramidal (1)")
cellgrpsdf <- rbind(cellgrpsdf,Nctx_npyr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(Nctx_npyr_1)

ss_pyr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 15)
ss_pyr_1$cellgrp <- apply(ss_pyr_1,1,function(r) "SS pyramidal (1)")
cellgrpsdf <- rbind(cellgrpsdf,ss_pyr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(ss_pyr_1)

sal_gang <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 1)
sal_gang$cellgrp <- apply(sal_gang,1,function(r) "Retinal ganglion (1)")
cellgrpsdf <- rbind(cellgrpsdf,sal_gang)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(sal_gang)

olf_pyr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 28)
olf_pyr_1$cellgrp <- apply(olf_pyr_1,1,function(r) "OB pyramidal (1)")
cellgrpsdf <- rbind(cellgrpsdf,olf_pyr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(olf_pyr_1)

olf_npyr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 27)
olf_npyr_1$cellgrp <- apply(olf_npyr_1,1,function(r) "OB non-pyramidal (1)")
cellgrpsdf <- rbind(cellgrpsdf,olf_npyr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(olf_npyr_1)

CA3_intr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 31)
CA3_intr_1$cellgrp <- apply(CA3_intr_1,1,function(r) "CA3 INs (1)")
cellgrpsdf <- rbind(cellgrpsdf,CA3_intr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(CA3_intr_1)

ss_intr_1 <- subset(sigdf[,c(33:64,72)],sigdf$classification == 1 & sigdf$mxdgroupnum == 19)
ss_intr_1$cellgrp <- apply(ss_intr_1,1,function(r) "SS INs (1)")
cellgrpsdf <- rbind(cellgrpsdf,ss_intr_1)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(ss_intr_1)

ret_gang <- subset(sigdf[,c(33:64,72)],sigdf$classification == 3)
ret_gang$cellgrp <- apply(ret_gang,1,function(r) "Retinal ganglion (3)")
cellgrpsdf <- rbind(cellgrpsdf,ret_gang)
cellgrpMeans[nrow(cellgrpMeans)+1,] <-add2cellgrpMeans(ret_gang)

nrow(cellgrpMeans)
unique(cellgrpsdf$cellgrp)
colnames(cellgrpsdf)

#compute z-scores
colnames(cellgrpsdf)
z1 <- computezscore(cellgrpsdf[,c(1:32,34)],"cellgrp")
length(z1)
dim(z1)
write.table(z1$c, file="cellgrpNames.txt",quote=FALSE, sep=",")

EucDistFrmGlobalmean <- computezscoreDist(z1)
dim(EucDistFrmGlobalmean)
#compute variance
v1 <- computegrpvar(cellgrpsdf[,c(1:32,34)],"cellgrp")
length(v1)
dim(v1)
grpvarmag <- computevarscrMag(v1)



## Increase bottom margin to make room for rotated labels
par(mar = c(7, 4, 4, 2) - 0.5)
#DD <- table(rpois(100, lambda=5))
#names(DD) <- paste("long", names(DD), sep="_")
## Plot, but suppress the labels
#midpts <- barplot(DD,names.arg="")
## Create plot with no x axis and no x axis label
mp <- barplot(EucDistFrmGlobalmean$EucdistanceFromMu,
              main="Euclidean distance from clusters to global mean",cex.main=1,
              xaxt="n")
## Set up x axis with tick marks alone
axis(1, labels = FALSE)
## Create some text labels
labels <- paste(rownames(EucDistFrmGlobalmean))
## Plot x axis labels at default tick marks
#
text(1:32, y = par("usr")[3] - 0.25,cex=0.6,srt =45,adj=1,
     labels = labels, xpd = TRUE)

#barplot of variances
par(mar = c(7, 4, 4, 1) + 0.2)
## Create plot with no x axis and no x axis label
mp <- barplot(grpvarmag$magnitude,
              main="variance of each group",cex.main=1,
              xaxt="n")
## Set up x axis with tick marks alone
axis(1, labels = FALSE)
## Create some text labels
labels <- paste(rownames(grpvarmag))
## Plot x axis labels at default tick marks
#
text(1:32, y = par("usr")[3] - 0.05,cex=0.6,srt =45,adj=1,
     labels = labels, xpd = TRUE)

#plot z-scores
plot(unlist(z1[[1,2]]),type="b",ylim = c(1,24), xlim=c(1,32),
     xlab="principal components",ylab="distances from the data centroid",
     main="z-scores from the centroid",col="red")
lines(unlist(z1[[2,2]]),type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(z1[[3,2]]),type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(z1[[4,2]]),type="b",lty=4,pch=12)
legend("topright",text.width = strwidth("1,000000000"),legend = z1$c, pch=c(1,0,6,12), cex=0.7, col=c(2,rep(1,3)))

#plot barplot of z-scores


#populate the cellgrpMeans with group means and SD
add2cellgrpMeans <- function(grp){
  t <- data.frame(c=character(0), mu=I(vector('list', 0))
                             ,sd=I(vector('list', 0)),stringsAsFactors=FALSE)
  t[1,'c'] <- grp$cellgrp[1]
  print(t[1,'c'])
  t[[1,'mu']] <- as.matrix(aggregate(grp[,1:32],by=list(grp$cellgrp),FUN=mean)[2:33])
  t[[1,'sd']] <-  as.matrix(aggregate(grp[,1:32],by=list(grp$cellgrp),FUN=sd)[2:33])
  return (t)
}


#pick the top two dimensions with largest zscr values on the main group
#extracting the 'PC2' and 'PC4' columns from the array of means
a <- sapply(cellgrpMeans$mu,function(i) i[c(2,4)])
minpc2 <- floor(min(a[1,]))
maxpc2 <- ceiling(max(a[1,]))
minpc4 <- floor(min(a[2,]))
maxpc4 <- ceiling(max(a[2,]))

plot(globalMean['PC2'],globalMean['PC4'],xlim=c(minpc2,maxpc2),ylim=c(minpc4,maxpc4),
     pch=1,xlab='branching density',ylab='local bifurcation torque',main="",cex.main=0.8)
segments(globalMean['PC2'], globalMean['PC4']-globalSD['PC4']/2,globalMean['PC2'],globalMean['PC4']+globalSD['PC4']/2)
segments(globalMean['PC2']-globalSD['PC2'],globalMean['PC4'],globalMean['PC2']+globalSD['PC2']/2,globalMean['PC4'])
plot(cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC2'],
     cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC4'],
     labels="basal pyramidal(cortical)",col='red',cex=0.7)
segments(cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC2'],
         cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC4']+
           cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),3]][,'PC4']/2,
         cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC2'],
         cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),2]][,'PC4']+
           cellgrpMeans[[which(cellgrpMeans$c == "basal pyramidal"),3]][,'PC4']/2)
text(cellgrpMeans[[which(cellgrpMeans$c == "CA1 interneurons"),2]][,'PC2'],
     cellgrpMeans[[which(cellgrpMeans$c == "CA1 interneurons"),2]][,'PC4'],
     labels="CA1 interneurons",cex=0.7)
text(cellgrpMeans[[which(cellgrpMeans$c == "sensory interneurons"),2]][,'PC2'],
     cellgrpMeans[[which(cellgrpMeans$c == "sensory interneurons"),2]][,'PC4'],
     labels="sensory interneurons",cex=0.7)
text(cellgrpMeans[[which(cellgrpMeans$c == "CA3 interneurons"),2]][,'PC2'],
     cellgrpMeans[[which(cellgrpMeans$c == "CA3 interneurons"),2]][,'PC4'],
     labels="CA3 interneurons",cex=0.7)

clusnums <- unique(tmpEM$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[c(1,4),clusnums[i]],shape=parms$variance$sigma[c(1,4),c(1,4),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}




#plot highlighting each cluster
plot(sigdf[,33],sigdf[,34],col="white",cex=.55,pch=16,xlab=names(sigdf)[33],
     main="Cluster 7 within component space of whole 4196 neurons",ylab=names(sigdf)[34],cex.main=0.8)
text(sigdf[,33],sigdf[,34],labels=as.character(sigdf$classification),col=rgb(0,100,0,50,maxColorValue=255),cex=1)
#text(cluster4[,33],cluster4[,36],labels=as.character(cluster4$classification),col="black",font=2,cex=0.8)
#text(cluster5[,33],cluster5[,34],labels=as.character(cluster5$classification),col="black",font=2,cex=.7)
#text(cluster6[,33],cluster6[,34],labels=as.character(cluster6$classification),col="black",font=2,cex=.7)
text(cluster7[,33],cluster7[,34],labels=as.character(cluster7$classification),col="black",font=2,cex=.7)
#text(cluster2[,33],cluster2[,34],labels=as.character(cluster2$classification),col="black",font=2,cex=.7)
#text(cluster1[,33],cluster1[,34],labels=as.character(cluster1$classification),col="black",font=2,cex=.7)

clusnums <- unique(tmpEM$classification)
parms <- tmpEM$parameters
for(i in 1:length(clusnums)){
  ellipse(parms$mean[c(1,4),clusnums[i]],shape=parms$variance$sigma[c(1,4),c(1,4),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

animate(cluster4[,33:37], grand_tour(), display_xy(col=1))

###Step#7: Choose groups with data significantly being in only one cluster###
clus7Grps <- c(4,5,6,7,9,10,11)
cluss6Grps <- c(8,12)
singleGrps <- c(3,13,14,15,16,17,21,22,23,24,25,37)
#metagrpList_1 <- c("neuron_name","cellclass3","region3","age_class","gender", "protocol","slice_thickness","slicing_direction","stain","pmid")



#impact of only pyramidal cells in all clusters
colnames(sigdf)
unique(sigdf$cellclass2)
sigdf_pyr <- subset(sigdf, sigdf$cellclass2=="Pyramidal cell")
dim(sigdf_pyr)

colnames(sigmetricdf)
unique(sigmetricdf$cellclass2)
sigmetricdf_pyr <- subset(sigmetricdf, sigmetricdf$cellclass2=="Pyramidal cell")
dim(sigmetricdf_pyr)


#plot cluster7 on size metrics
sigdf_1 <- sigdf#[,c(33:64,71,72)]
sigdf_1$PC1  <- -1 * sigdf_1$PC1
#sigdf_1$PC2  <- -1 * sigdf_1$PC2
colnames(sigdf_1)
dim(sigdf_1)
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdf_1)
sigdf_1[which(sigdf_1$PC1 <= -10),c("neuron_name")]
sigdf_1[which(sigdf_1$PC1 >= 7),c("neuron_name")]
PlotDoubleMetricDistrOnClusters('PC2','PC4',sigdf_1)
sigdf_1[which(sigdf_1$PC2>10),c("neuron_name","cellclass2","archive_name","classification")]
sigdf_1[which(sigdf_1$PC2>5),c("neuron_name","cellclass2","archive_name","classification")]
sigdf_1 <- subset(sigdf_1,sigdf_1$classification == 7)

sigMdf_1 <- sigmetricdf[,c(33:64,72)]
sigMdf_1$PC1  <- -1 * sigmetricdf$PC1
#sigdf_1$PC2  <- -1 * sigdf_1$PC2
colnames(sigMdf_1)

PlotDoubleMetricDistrOnGroups('PC1','PC2',sigdf_1)
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdf)
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdf_pyr)
unique(sigdf_pyr$classification)
PlotDoubleMetricDistrOnClusters('ln(Length_total_sum)','ln(Branch_pathlength_avg)',sigmetricdf)
PlotDoubleMetricDistrOnGroups('ln(Length_total_sum)','ln(Branch_pathlength_avg)',clusterM4)

PlotDoubleMetricDistrOnGroups('ln(Length_total_sum)','ln(Branch_pathlength_avg)',clusterM7)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clusterM7)
PlotDoubleMetricDistrOnGroups('Bif_ampl_remote_avg','ln(Length_total_sum)',clusterM7)

unique(clusterM2[,c("archive_name","mxdgroupnum")])
clusterM27 <- subset(sigmetricdf, sigmetricdf$classification == 7 | sigmetricdf$mxdgroupnum %in% c(6,13,17))
dim(clusterM27)
colnames(clusterM27)
PlotDoubleMetricDistrOnClusters('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clusterM27)

clusterM27_all <- subset(sigmetricdf, sigmetricdf$classification %in% c(7,2))
dim(clusterM27_all)
PlotDoubleMetricDistrOnClusters('Bif_ampl_remote_avg','ln(Branch_pathlength_avg)',clusterM27_all)



plotDensitiesForAllClusters('ln(EucDistance_avg)')

plotDensitiesOnPCsForAllClusters(sigdf, 'PC1')
plotDensitiesOnPCsForAllClusters('PC2')
plotDensitiesOnPCsForAllClusters('PC7')

plotDensitiesForAllClusters('ln(Length_total_sum)')
plotDensitiesForAllClusters('Bif_torque_local_avg')
plotDensitiesForAllClusters('ln(Contraction_avg)')
plotDensitiesForAllClusters('Bif_ampl_remote_avg')
plotDensitiesForAllClusters('ln(Soma_Surface_total_sum)')
plotDensitiesForAllClusters('Partition_asymmetry_avg')
metricCor <- cor(sigmetricdf[,33:64])
dim(metricCor)
write.table(metricCor, file="metricCor.txt",quote=FALSE, sep=",")

plotDensitiesOnPCsForAllClusters('PC2')

cluster6$PC24 <- -1 * cluster6$PC24
PlotDoubleMetricDistrOnGroups('PC2','PC24',cluster6)
PlotDoubleMetricDistrOnGroups('PC2','PC24',cluster6)
#2 monkey pyramidal cells >2.0
cluster6[which(cluster6$PC2 > 2),]
#2 magnopyramidal cells >-5
cluster6[which(cluster6$PC2 < -5),]

sigdf_pyr$PC24 <- -1 * sigdf_pyr$PC24
#sigmetricdf_pyr$ln(Branch_pathlength_avg) <- -1 * sigmetricdf_pyr$ln(Branch_pathlength_avg)

#only pyramidal cells from cluster 6
clustermetric6_pyr <- subset(sigmetricdf_pyr,sigmetricdf_pyr$classification %in% c(6))
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clustermetric6_pyr)

#only pyramidal cells from cluster 2
cluster2_pyr <- subset(sigdf_pyr,sigdf_pyr$classification %in% c(2))
#separate basal only groups
cluster2_basalpyr <- subset(cluster2_pyr, cluster2_pyr$mxdgroupnum %in% c(6,13))
#separate basal & apical groups
cluster2_pyr <- subset(cluster2_pyr, cluster2_pyr$mxdgroupnum %in% c(28,25))
dim(cluster2_pyr)
PlotDoubleMetricDistrOnGroups('PC2','PC24',cluster2_pyr)
PlotDoubleMetricDistrOnGroups('PC2','PC24',cluster2_basalpyr)
#cluster pyramidal cell metrics
clustermetric2_pyr <- subset(sigmetricdf_pyr,sigmetricdf_pyr$classification %in% c(2))
#test
clustermetric2_pyr[which(clustermetric2_pyr$mxdgroupnum %in% c(28,25)),c("neuron_name","mxdgroupnum","archive_name","region3","pmid")]
clustermetric2_pyr[which(clustermetric2_pyr$mxdgroupnum %in% c(6,13)),c("neuron_name","mxdgroupnum","archive_name","region3","pmid")]
unique(clustermetric2_pyr$mxdgroupnum)
unique(cluster2_pyr$mxdgroupnum)
unique(clustermetric2_basalpyr$mxdgroupnum)
#separate basal only groups (6 and 13)
clustermetric2_basalpyr <- subset(clustermetric2_pyr, clustermetric2_pyr$mxdgroupnum %in% c(6,13))
#separate apical & basal groups (28 and 25)
clustermetric2_pyr <- subset(clustermetric2_pyr, clustermetric2_pyr$mxdgroupnum %in% c(28,25))
dim(clustermetric2_pyr)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clustermetric2_pyr)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clustermetric2_basalpyr)

unique(cluster2_pyr[,c("mxdgroupnum","archive_name","region3")])
cluster2_pyr[which(cluster2_pyr$PC24>0.9),c("neuron_name","mxdgroupnum","archive_name","region3","pmid")]
cluster2_pyr[which(cluster2_pyr$PC2>5),c("neuron_name","mxdgroupnum","archive_name","region3","pmid")]
PlotDoubleMetricDistrOnClusters('PC2','PC24',sigdf_pyr)

cluster67_size <- subset(sigdf_1, sigdf_1$classification %in% c(6,7))
PlotDoubleMetricDistrOnClusters('PC1','PC32',cluster67_size)

cluster6732_size <- subset(sigdf_1, sigdf_1$classification %in% c(6,7,3,2))
PlotDoubleMetricDistrOnClusters('PC1','PC2',cluster6732_size)

clusterM67 <- subset(sigmetricdf,sigmetricdf$classification %in% c(6,7))
dim(clusterM67)
PlotDoubleMetricDistrOnClusters('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM67)
PlotDoubleMetricDistrOnClusters('ln(Width_total_sum)','ln(N_tips_total_sum)',clusterM67)
PlotDoubleMetricDistrOnClusters('ln(EucDistance_avg)','ln(N_tips_total_sum)',clusterM67)
PlotDoubleMetricDistrOnClusters('ln(EucDistance_avg)','ln(N_bifs_total_sum)',clusterM67)
PlotDoubleMetricDistrOnClusters('ln(Branch_Order_max)','ln(N_bifs_total_sum)',clusterM67)
colnames(clusterM7)
cluster67 <- subset(sigdf,sigdf$classification %in% c(6,7))
N_cluster67 <- reEMonSubsets(cluster67)
N_cluster7 <- reEMonSubsets(cluster7)

#plot cluster6 on PC2 and PC3
colnames(clusterM6)
PlotDoubleMetricDistrOnClusters('ln(Branch_pathlength_avg)','ln(Contraction_avg)',sigmetricdf)
PlotDoubleMetricDistrOnClusters('Bif_ampl_remote_avg','ln(Contraction_avg)',sigmetricdf)
PlotDoubleMetricDistrOnClusters('ln(Branch_pathlength_avg)','ln(Contraction_avg)',sigmetricdf)
PlotDoubleMetricDistrOnClusters('Bif_ampl_remote_avg','ln(Branch_pathlength_avg)',sigmetricdf)
PlotDoubleMetricDistrOnClusters('PC24','PC2',sigdf)
cluster4$PC1<- -1 * cluster4$PC1
#cluster4$PC2<- -1 * cluster4$PC2
PlotDoubleMetricDistrOnGroups('PC1','PC2',cluster4)
PlotDoubleMetricDistrOnGroups('PC24','PC2',cluster4)
unique(cluster4[,c('archive_name','cellclass2','cellclass3')])
cluster4[which(cluster4$PC7 == max(cluster4$PC7)),]

       #highly dense and greater length
cluster4[which(cluster4$PC1 > 7.75),]
#sparse and greater length
cluster4[which(cluster4$PC1 > 6 & cluster4$mxdgroupnum == 36),]
#sparse and small length
cluster4[which(cluster4$PC1 < -4 & cluster4$mxdgroupnum == 23),]
#in between 23 and 36
cluster4[which(cluster4$PC1 > 0 & cluster4$PC1 < 1 & cluster4$mxdgroupnum == 12),]


N_cluster6 <- reEMonSubsets(cluster6)

#plot cluster3 on PC1 and PC2
cluster3$PC6 <- -1 * cluster3$PC6
cluster3$PC1 <- -1 * cluster3$PC1
unique(cluster3[,c('mxdgroupnum','species','archive_name')])
PlotDoubleMetricDistrOnGroups('PC1','PC6',cluster3)
PlotDoubleMetricDistrOnClusters('PC1','PC6',sigdf)
PlotDoubleMetricDistrOnClusters('PC1','PC3',sigdf_1)

PlotDoubleMetricDistrOnClusters('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM3)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clusterM3)
PlotDoubleMetricDistrOnGroups('Partition_asymmetry_avg','ln(EucDistance_max)',clusterM3)
#PlotDoubleMetricDistrOnClusters('Partition_asymmetry_avg','ln(EucDistance_max)',sigmetricdf)

#plot cluster2 on PC1 and PC2
PlotDoubleMetricDistrOnGroups('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM2)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clustermetric2_pyr)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','ln(Contraction_avg)',clusterM2)
PlotDoubleMetricDistrOnGroups('Bif_torque_local_avg','Bif_ampl_remote_avg', clusterM2)
PlotDoubleMetricDistrOnGroups('ln(Soma_Surface_total_sum)','Partition_asymmetry_avg', clusterM2)
N_cluster2 <- reEMonSubsets(cluster2)

#plot cluster#7 and just Insular cells in cluster#2
clusterM72 <- subset(sigmetricdf,sigmetricdf$classification == 7 | 
  sigmetricdf$classification == 2 & sigmetricdf$lobes == 'Insula')
dim(clusterM72)
PlotDoubleMetricDistrOnClusters('ln(Branch_Order_max)','ln(Branch_pathlength_avg)', clusterM72)
PlotDoubleMetricDistrOnGroups('Bif_torque_local_avg','Bif_ampl_remote_avg', clusterM72)

#plot cluster5 on PC1 and PC2
PlotDoubleMetricDistrOnGroups('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM5)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','ln(Contraction_avg)',clusterM5)
PlotDoubleMetricDistrOnGroups('Bif_torque_local_avg','Bif_ampl_remote_avg', clusterM5)
PlotDoubleMetricDistrOnGroups('ln(Soma_Surface_total_sum)','Partition_asymmetry_avg', clusterM5)
PlotDoubleMetricDistrOnClusters('ln(Soma_Surface_total_sum)','Partition_asymmetry_avg', sigmetricdf)
cluster5 <- subset(cluster5, cluster5$cellclass2 == "Pyramidal cell")
cluster5 <- subset(cluster5, cluster5$mxdgroupnum %in%  c(16,26,20,35,36))
clusterM5_pyr <- subset(clusterM5, clusterM5$cellclass2 == "Pyramidal cell")
clusterM5_pyr <- subset(clusterM5, clusterM5$mxdgroupnum %in% c(16,26,20,35,36))
dim(clusterM5_pyr)
PlotDoubleMetricDistrOnGroups('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM5_pyr)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','ln(Contraction_avg)',clusterM5_pyr)
PlotDoubleMetricDistrOnGroups('Bif_torque_local_avg','Bif_ampl_remote_avg', clusterM5_pyr)
PlotDoubleMetricDistrOnGroups('ln(N_bifs_total_sum)','Partition_asymmetry_avg', clusterM5_pyr)
PlotDoubleMetricDistrOnGroups('ln(N_bifs_total_sum)','ln(EucDistance_max)', clusterM5_pyr)
dim(cluster5)
unique(cluster5$cellclass2)
N_cluster5 <- reEMonSubsets(cluster5)

#plot density distributions for cluster#1
PlotDoubleMetricDistrOnGroups('ln(Branch_Order_max)','ln(Branch_pathlength_avg)',clusterM1)

PlotDoubleMetricDistrOnClusters('ln(Width_total_sum)','ln(N_bifs_total_sum)', sigmetricdf)
PlotDoubleMetricDistrOnClusters('ln(Width_total_sum)','Bif_torque_local_avg', sigmetricdf)
PlotDoubleMetricDistrOnClusters('Fractal_Dim_max','Bif_tilt_remote_avg', subset(origdenddf,origdenddf$classification %in% c(4,5)))

sigmetricdf[which(sigmetricdf$'ln(Width_total_sum)' < 4 & sigmetricdf$classification == 1),]
sigmetricdf[which(sigmetricdf$'ln(Width_total_sum)' > 6 & sigmetricdf$classification == 7),]
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_torque_local_avg',clusterM1)

#plot for cluster4
PlotDoubleMetricDistrOnGroups('ln(Width_total_sum)','ln(Soma_Surface_total_sum)',clusterM4)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','Bif_ampl_remote_avg',clusterM4)
PlotDoubleMetricDistrOnClusters('ln(Width_total_sum)','ln(Soma_Surface_total_sum)',sigmetricdf)
PlotDoubleMetricDistrOnGroups('ln(Width_total_sum)','Partition_asymmetry_avg', clusterM4)
PlotDoubleMetricDistrOnGroups('ln(Branch_pathlength_avg)','ln(Contraction_avg)',clusterM4)
unique(clusterM4[,c("archive_name","cellclass3","mxdgroupnum","region3")])


#plot density distributions of single metric
plotDensitiesForAllClusters <-  function(sigdf, metric, group){
  plot(density(sigmetricdf[,metric]),xlim=c(min(sigmetricdf[,metric]), max(sigmetricdf[,metric])), 
       xlab=metric,cex.main=0.8, main=paste(metric, "on all 4196 neurons"))
  #lines(density(clusterM7[,metric]), lwd=2, col=rgb(0, 1, 0, 0.5))#,freq=FALSE)
  lines(density(clusterM4[,metric]), lwd=2, col=rgb(0, 0, 1, 0.5))#,freq=FALSE)
  lines(density(clusterM6[,metric]), lwd=2, col=rgb(1, 0, 1, 0.5))#,freq=FALSE)
  lines(density(clusterM3[,metric]), lwd=2, col=rgb(1, 0, 0, 0.5))#,freq=FALSE)
  lines(density(clusterM5[,metric]), lwd=2, col=rgb(1, 1, 0, 0.5))#,freq=FALSE)
  lines(density(clusterM2[,metric]), lwd=2, col=rgb(0, 1, 1, 0.5))#,freq=FALSE)
  lines(density(clusterM1[,metric]), lwd=2, col=rgb(0.7, 0.3, 0.9, 0.5))#,freq=FALSE)
  
  legend('topleft',c('cluster 7','cluster 4','cluster 6','cluster 3','cluster 5','cluster 2','cluster 1'),
         fill = c(rgb(0, 1, 0, 0.5),rgb(0, 0, 1, 0.5),rgb(1, 0, 1, 0.5),rgb(1, 0, 0, 0.5),rgb(1, 1, 0, 0.5),rgb(0, 1, 1, 0.5),rgb(0.7, 0.3, 0.9, 0.5)), bty = 'n',
         border = NA)
}

plotDensitiesOnCategories <- function(sigdf, metric, group, grpname,unts=NULL,highlight=NULL){
  if(!missing(unts)){
    unts = paste(metric,'(',unts,')')
  }else{
    unts = metric
  }
  #first plot the density plot for entire data
  plot(density(sigdf[,metric]),xlim=c(min(sigdf[,metric]), max(sigdf[,metric])), 
       xlab=unts,cex.main=0.8, lwd=3, main=paste(metric, "on ",nrow(sigdf),"neurons"))
  clusgrps <- unique(sigdf[,group])
  labelcol <- list()
  lgnd <- list()
  #overlay the group plots on the top
  for(i in 1:length(clusgrps)){
    tmpclus <- subset(sigdf, sigdf[,group] == clusgrps[i])
    lines(density(tmpclus[,metric]), lwd=2, col=mypalette_50[i])
    labelcol <- append(labelcol,mypalette_50[i])
    lgnd <- append(lgnd,tmpclus[,grpname][1])
  }
  print(unlist(labelcol))
  #add legend
  legend("topleft",text.width = strwidth("1,0"), 
         legend = lgnd, pch=19, cex = 0.7, col=unlist(labelcol))
}

PlotDoubleMetricDistrOnClusters <- function(m1, m2, df,highlight = NULL){
  Xname <- m1
  Yname <- m2
  print(highlight)
  if(!missing(highlight)){
    hdf <- subset(df, df$classification %in% highlight)
    print(unique(hdf$classification))
    #clusnums <- highlight
  }else{
    hdf <- df
    #clusnums <- unique(hdf$classification)
  }
  featurecol <- hdf$classification
  clusnums <- unique(featurecol)
  #generate alphabet code for original data set
  alphaind <- alphabetcode(groupnum=df$classification)
  #grpttl <- paste(length(clusnums),'clusters') #toString(unique(df$classification, sep=","))
  plot(df[,Xname],df[,Yname],
       col="white",cex=.55,pch=16,xlab=Xname,ylab=Yname,cex.main=0.8)#main=grpttl,
  text(df[,Xname],df[,Yname],labels=alphaind,col="grey",cex=.7)
  #color code for highlighted data sets
  colarr <- colorcode(groupnum=featurecol)
  #re-assign alphabet codes for highlighted data sets
  alphaind <- alphabetcode(groupnum=featurecol)
  text(hdf[,Xname],hdf[,Yname],labels=alphaind,col=mypalette[colarr],cex=.7)
  #legend("bottomright",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))
}

PlotDoubleMetricDistrOnClustersNGroups <- function(m1, m2,m1desc,m2desc,df,highlight=NULL){
  Xname <- m1
  Yname <- m2
  if(!missing(highlight)){
    hdf <- subset(df, df$cmpoundGrps %in% highlight)
  }else{
    hdf <- df
  }
  featurecol <- hdf$cmpoundGrps
  grplabels <- hdf$cmpoundGrpName
  labelnum <- df$classification
  alphaind <- alphabetcode(groupnum=labelnum)
  clusnums <- unique(labelnum)
  #grpttl <- paste(length(unique(labelnum)),' clusters and ',length(unique(featurecol)),' groups') #toString(unique(df$classification, sep=","))
  colarr <- colorcode(groupnum=featurecol, highlight)
  ymax <- max(df[,m2])
  print(ymax)
  print(unique(colarr))
  plot(df[,Xname],df[,Yname],
       col="white",cex=.55,pch=16,xlab=paste(Xname,'(',m1desc,')'),ylab=paste(Yname,'(',m2desc,')'),cex.main=0.7)#main=ttl)
  #plot all points in grey
  text(df[,Xname],df[,Yname],labels=alphaind,col="grey",cex=.7)#col=mypalette[colarr]    alphaind <- alphabetcode(groupnum=hdf$classification, length(hdf$classification))
  #plot particular groups to highlight
  labelnum <- hdf$classification
  alphaind <- alphabetcode(groupnum=labelnum)
  print(length(alphaind))
  text(hdf[,Xname],hdf[,Yname],labels=alphaind,col=mypalette[colarr],cex=.9)
  #legend("bottomright",text.width = strwidth("1,0000"), legend = unique(grplabels), pch=19, cex=0.7, col=unique(mypalette[colarr]))
}

plotmetricHist <- function(origdf,metric,xlabel,groups,unts=NULL,binsize=30,logT=FALSE,round2=0,scaleupby=NULL,legendpos='topright'){
  
  if(logT==TRUE){
    origdf[,metric] <- logTransformData(origdf[,metric])
  }
  #get the xlim range from the log transformation of the original data for the given metric
  numTics <- 11
  #axpos <- seq(from=min(origdf[,metric]),to=max(origdf[,metric]),by=diff(range(origdf[,metric]))/numTics)
  #axpos <- origdf[,metric]
  #print(min(origdf[,metric]))
  #print(max(origdf[,metric]))
  #axpos <- c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7)
  axpos <- c(0,0.15,0.30,0.45,0.60,0.75)
  #axpos <- seq(from=40,to=180,by=20)
  print(axpos)
  #print(exp(axpos))
  if(logT){
    axval <- signif(exp(axpos),0)
    #axval <- c(20,40,60,80,100,250,400,550,700)
  }else{
    axval <- axpos
  }
  #axval <- axpos
  print("axval")
  print(axval)
  
  if(!missing(unts)){
    xttl <- paste(metric, ' (',unts,')')
    print(xttl)
  }else{
    xttl <- metric
  }
 
  hist(origdf[,metric], breaks=binsize,xaxt='n', xlab=xlabel, ylab = "Count", col="grey",cex.main=0.6,main="")#group distributions")
  #set specific values when needed override the following array
  #axval <- seq(from=40,to=180,by=20)
  #axval <- c(200,300,400,500,600,700,800,900,1000,1100)
  #ticks
  #axpos <- log(axval)
  axis(side=1,at=axpos,labels=round(axval,round2),cex=0.6)
  #minor.tick(nx=5,tick.ratio=0.20)
  
  #highlight the groups of interest by overlaying them on the background
  hdf <- subset(origdf, origdf$cmpoundGrps %in% groups)

  cols <- list()
  colnames <- list()
  for(i in 1:length(groups)){
   tmp <- subset(origdf, origdf$cmpoundGrps==groups[i])
   
   name <- origdf[which(origdf$cmpoundGrps == groups[i]),"cmpoundGrpName"][1]
   
   
   #carr <- col2rgb(mypalette[i])
   col1 <- mypalette[i]
   cols <- append(cols,col1)

   newp<- hist(tmp[,metric],breaks = binsize,plot=FALSE)#add = T,col = rgb(t(col2rgb(col1)/255), alpha = 0.5))
   print("counts & breaks")
   print(newp$counts)
   print(newp$breaks)
   if(!missing(scaleupby)){
     print(paste("scaling by factor ",scaleupby[i]," for ",name))
     print(scaleupby[i])
     print(unique(tmp[,metric]))
     
     print(aggregate(tmp[,metric],
                     by=list(cut(tmp[,metric],newp$breaks,include.lowest=TRUE)),
                     FUN=length))
     print(newp$counts[newp$counts > 0])
     #print(list(cut(tmp[,metric],newp$breaks)))
     #scale the counts for each group by its own scaling factor
     newp$counts[newp$counts > 0] <- aggregate(tmp[,metric],
                                               by=list(cut(tmp[,metric],newp$breaks,include.lowest=TRUE)),
                                               FUN=length)$x*scaleupby[i]
     
     colnames <- append(colnames,paste(name,'x',scaleupby[i]))
   }else{
     #re-calculates the counts without changing their values
     newp$counts[newp$counts > 0] <- aggregate(tmp[,metric],
                                               by=list(cut(tmp[,metric],newp$breaks)),
                                               FUN=length)$x
     colnames <- append(colnames,name)
   }
  
   print(newp$counts)
   #the new ticks values
   #axpos <- seq(from=floor(min(newp$counts)),to=ceiling(max(newp$counts)),
   #              by=diff(range(newp$counts))/(length(newp$breaks)-1))
   #axpos <- seq(from=40,to=180,by=20)
  
   plot(newp,add=T,col=rgb(t(col2rgb(col1)/255), alpha = 0.5),axes=FALSE)
   #d <- density(newp$counts)
   #lines(x=d$x,y=d$y,lwd=2,col=mypalette_50[i])
   #plot(x=newp$mids,y=newp$density,add=T,type="l", xaxt="n")
   
   #axis(4)
  }
  #axis(side=4,at=axpos,labels=round(axpos,0),cex=0.7)
  print(unlist(colnames))
  print(unlist(cols))
  #legend(legendpos,unlist(colnames),fill = unlist(cols), bty = 'n',border = NA) 
}
library(Hmisc)
PlotDoubleMetricDistrOnGroups <- function(m1, m2, df){
  Xname <- m1
  Yname <- m2
  grpttl <- toString(unique(df$classification, sep=","))
  featurecol <- df$mxadgroupnum`
  colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
  plot(df[,Xname],df[,Yname],
       col="white",cex=.55,pch=16,xlab=Xname,main=paste("clusters:", grpttl, "on metrics"),ylab=Yname,cex.main=0.8)
  text(df[,Xname],df[,Yname],labels=as.character(df$mxdgroupnum),col=mypalette[colarr],cex=.7)
  legend("bottomright",text.width = strwidth("1,0"),
         legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))
}

#re-do EM on given subset and find how many new models are discovered
reEMonSubsets <- function(df){
  dim(df)
  reEM <- Mclust(df[,33:50],G=1:15)
  summary(reEM)
  plot(reEM, data=df[,33:50], what='BIC')
  df$classification <-sdfdsfsdf reEM$classification
  unique(df$classification)
  length(df$mxdgroupnum)
  dim(df)
  aggregate(df,by=list(df[,c("mxdgroupnum")]),FUN=length)
  tabtmp <- aggregate(df$neuron_name, by=sdfslist(row=df$mxdgroupnum,
                                               col=df$classification), FUN=length)
  tabtmp <- tabtmp[order(tabtmp$x,decreasing=T),]
  print(tabtmp)
  return (df)
}


#plot histogram of entire data and overlay cluster7 on branch_order_max
plot(density(sigmetricdf[,'ln(Branch_Order_max)']),xlim=c(min(sigmetricdf$'ln(Branch_Order_max)'), max(sigmetricdf$'ln(Branch_Order_max)')), 
     xlab="ln(Branch_Order_max)",cex.main=0.8, main="Branch order on all 4196 neurons")#, freq = FALSE)
lines(density(clusterM7$'ln(Branch_Order_max)'), lwd=2, col=rgb(0, 1, 0, 0.5))#,freq=FALSE)
lines(density(clusterM4$'ln(Branch_Order_max)'), lwd=2, col=rgb(0, 0, 1, 0.5))#,freq=FALSE)
lines(density(clusterM6$'ln(Branch_Order_max)'), lwd=2, col=rgb(1, 0, 1, 0.5))#,freq=FALSE)
lines(density(clusterM3$'ln(Branch_Order_max)'), lwd=2, col=rgb(1, 0, 0, 0.5))#,freq=FALSE)
lines(density(clusterM5$'ln(Branch_Order_max)'), lwd=2, col=rgb(1, 1, 0, 0.5))#,freq=FALSE)
lines(density(clusterM2$'ln(Branch_Order_max)'), lwd=2, col=rgb(0, 1, 1, 0.5))#,freq=FALSE)
lines(density(clusterM1$'ln(Branch_Order_max)'), lwd=2, col=rgb(0.7, 0.3, 0.9, 0.5))#,freq=FALSE)

legend('topright',c('cluster 7','cluster 4','cluster 6','cluster 3','cluster 5','cluster 2','cluster 1'),
       fill = c(rgb(0, 1, 0, 0.5),rgb(0, 0, 1, 0.5),rgb(1, 0, 1, 0.5),rgb(1, 0, 0, 0.5),rgb(1, 1, 0, 0.5),rgb(0, 1, 1, 0.5),rgb(0.7, 0.3, 0.9, 0.5)), bty = 'n',
       border = NA)

#plot histogram of entire data and overlay cluster7 on branch_pathlength_avg
plot(density(sigmetricdf$'ln(Branch_pathlength_avg)'), xlim=c(min(sigmetricdf$'ln(Branch_pathlength_avg)'), max(sigmetricdf$'ln(Branch_pathlength_avg)')), 
     xlab="ln(Branch_pathlength_avg)",cex.main=0.8, main="Average branch pathlength on all 4196 neurons")#,freq = FALSE)
lines(density(clusterM7$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(0, 1, 0, 0.5))#,freq = FALSE)
lines(density(clusterM4$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(0, 0, 1, 0.5))#,freq = FALSE)
lines(density(clusterM6$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(1, 0, 1, 0.5))#,freq = FALSE)
lines(density(clusterM3$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(1, 0, 0, 0.5))#,freq = FALSE)
lines(density(clusterM5$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(1, 1, 0, 0.5))#,freq = FALSE)
lines(density(clusterM2$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(0, 1, 1, 0.5))#,freq = FALSE)
lines(density(clusterM1$'ln(Branch_pathlength_avg)'), lwd=2,col=rgb(0.7, 0.3, 0.9, 0.5))#,freq = FALSE)

legend('topleft',c('cluster 7','cluster 4','cluster 6','cluster 3','cluster 5','cluster 2','cluster 1'),
       fill = c(rgb(0, 1, 0, 0.5),rgb(0, 0, 1, 0.5),rgb(1, 0, 1, 0.5),rgb(1, 0, 0, 0.5),rgb(1, 1, 0, 0.5),rgb(0, 1, 1, 0.5),rgb(0.7, 0.3, 0.9, 0.5)), bty = 'n',
       border = NA)

###step#2: plot all groups within the same cluster with colors identifying the groups###
getSignificantGroupsNClusters(sigcond_over,cluster6)

#blowfly, CA1, CA3 interneurons 
tmp <- subset(cluster4, cluster4$mxdgroupnum %in% c(3,23,31,12,18,20,26,36))
dim(tmp)
colnames(tmp)
clus5grps <- c(14,16,18,19,20,24,26,33,34,35,36)
length(clus5grps)
tmp <- subset(cluster5, cluster5$mxdgroupnum %in% c(14,16,18,19,20,24,26,33,34,35,36))
dim(tmp)


aggregate(tmp$mxdgroupnum,by=list(tmp$mxdgroupnum),FUN=length)
#choose group number for color
featurecol <- sigdf$classification
unique(featurecol)
featurecol <- featurecol[order(featurecol)]
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))

colarr[colarr %in% c(1,2,3,5,6,7)] <- 1
colarr[colarr %in% c(4)] <- 2
#plot on PC1 and PC2
Xname <- 'PC1'
Yname <- 'PC4'
grpN <- 'all clusters'
grpttl <- 'cluster 4 vs others'
sigdf <- sigdf[order(sigdf$classification),]

plot(sigdf[,Xname],sigdf[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)

text(sigdf[,Xname],sigdf[,Yname],labels=as.character(sigdf$classification),col=mypalette[colarr],cex=.7)
#legend("bottomleft",text.width = strwidth("1,00"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#plot all on metrics
Xname <- "ln(Branch_Order_max)"#"ln(Soma_Surface_total_sum)"#"ln(Fractal_Dim_max)"#"ln(Terminal_degree_avg)"#
Yname <- "Bif_torque_local_avg"#"ln(Branch_pathlength_avg)"#"ln(PathDistance_max)"#"ln(Height_total_sum)"#"ln(Contraction_avg)"#
grpN <- 'all 7 clusters'
grpttl <- 'on metrics'
featurecol <- sigmetricdf$classification
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
colnames(sigmetricdf)
sigmetricdf[,Xname]
plot(sigmetricdf[,Xname],sigmetricdf[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(sigmetricdf[,Xname],sigmetricdf[,Yname],labels=as.character(sigmetricdf$classification),col=rgb(0,100,0,50,maxColorValue=255),cex=1)
text(clusterM4[,Xname],clusterM4[,Yname],labels=as.character(clusterM4$classification),col="black",cex=0.8)

#legend("bottomleft",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#plot cluster4 on metrics
Xname <- "ln(Branch_Order_max)"#"ln(Soma_Surface_total_sum)"#"ln(Fractal_Dim_max)"#"ln(Terminal_degree_avg)"#
Yname <-  "ln(Branch_pathlength_avg)"#"Bif_torque_local_avg"#"ln(PathDistance_max)"#"ln(Height_total_sum)"#"ln(Contraction_avg)"#
grpN <- 'cluster 4'
grpttl <- 'on metrics showing the groups'
featurecol <- clusterM4$mxdgroupnum
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))

clusterM4[,Xname]
plot(clusterM4[,Xname],clusterM4[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(clusterM4[,Xname],clusterM4[,Yname],labels=as.character(clusterM4$mxdgroupnum),col=mypalette[colarr],cex=.7)
legend("bottomright",text.width = strwidth("1,00"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#plot cluster5 on metrics
Xname <- "ln(Branch_Order_max)"#"ln(Soma_Surface_total_sum)"#"ln(Fractal_Dim_max)"#"ln(Terminal_degree_avg)"#
Yname <-  "ln(Branch_pathlength_avg)"#"Bif_torque_local_avg"#"ln(PathDistance_max)"#"ln(Height_total_sum)"#"ln(Contraction_avg)"#
grpN <- 'cluster 5'
grpttl <- 'on metrics showing the groups'
clusterM5 <- subset(sigmetricdf,sigmetricdf$classification == 5)
colnames(clusterM5)
featurecol <- clusterM5$mxdgroupnum
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
plot(clusterM5[,Xname],clusterM5[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(clusterM5[,Xname],clusterM5[,Yname],labels=as.character(clusterM5$mxdgroupnum),col=mypalette[colarr],cex=.7)
legend("bottomright",text.width = strwidth("1,00"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

###step#3: do re-EM on the subset###
library(fpc)
library(car)

getSignificantGroupsNClusters(tabtmp[,c("row","col")],tmp)

dim(Axondatamatrix)
colnames(Axondatamatrix)
dim(tmp)
tmpAx <- subset(Axondatamatrix, Axondatamatrix$neuron_name %in% tmp$neuron_name)
dim(tmpAx)
tmpAx$neuron_name
tmpWhole <- subset(tmp,tmp$neuron_name %in% tmpAx$neuron_name)
dim(tmpWhole)
unique(tmpWhole$mxdgroupnum)
colnames(tmpAx)
colnames(tmpWhole)
rownames(tmpWhole) <- tmpWhole$neuron_name
rownames(tmpAx)
replaceCols <- c("ct_groupnum","ct_groupttl","br_groupnum","br_groupttl","sp_groupnum","sp_groupttl","mxdgroupnum","classification","pmid")
length(replaceCols)
tmpAx <- tmpAx[order(tmpAx$neuron_name),]
tmpWhole <- tmpWhole[order(tmpWhole$neuron_name),]
#tmpAx[,"mxdgroupnum"] <- tmpWhole$mxdgroupnum
tmpAx[,replaceCols] <- tmpWhole[,replaceCols]
tmpAx[with(tmpAx,which(tmpAx$mxdgroupnum==36)),c("neuron_name","archive_name")]
tmpWhole[with(tmpWhole,which(tmpWhole$mxdgroupnum==36)),c("neuron_name","archive_name")]

aggregate(tmpAx$neuron_name,by=list(tmpAx[,c("mxdgroupnum")]),FUN=length)
aggregate(tmpWhole$neuron_name,by=list(tmpWhole[,c("mxdgroupnum")]),FUN=length)

#plot with whole arbors
tmp$neuron_name

plot(tmp[,Xname],tmp[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(tmp[,Xname],tmp[,Yname],labels=as.character(tmp$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

clusnums <- unique(reEM$classification)
parms <- reEM$parameters

for(i in 1:length(clusnums)){
  XYpcs <- tmp[with(tmp,which(tmp$classification == clusnums[i])),c(Xname,Yname)]
  XYmeans <- colMeans(XYpcs)
  
  #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
  #       radius=1, center.pch=as.character(clusnums[i]),col=1)
  ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

#plot with axons only
featurecol <- tmpAx$mxdgroupnum
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
plot(tmpAx[,Xname],tmpAx[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(tmpAx[,Xname],tmpAx[,Yname],labels=as.character(tmpAx$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))


#plot with dendrite arbors only
dim(dendriteReduced)
colnames(dendriteReduced)
dim(tmp)
colnames(tmp)
tmpWhole <- subset(tmp,tmp$neuron_name %in% rownames(dendriteReduced))
dim(tmpWhole)
dendritePCA <- subset(dendriteReduced, rownames(dendriteReduced)%in% tmpWhole$neuron_name)
dim(dendritePCA)
tmp$neuron_name
rownames(dendritePCA)
colnames(dendritePCA)
dendritePCA <- data.frame(tmpWhole[,1:30],dendritePCA[,1:31],tmpWhole[,65:73])
dim(dendritePCA)
tmpDend <- dendritePCA
dim(tmpDend)
colnames(tmpDend)

tmpDend <- tmpDend[order(tmpDend$neuron_name),]
tmpWhole <- tmpWhole[order(tmpWhole$neuron_name),]
#tmpAx[,"mxdgroupnum"] <- tmpWhole$mxdgroupnum
tmpDend[,replaceCols] <- tmpWhole[,replaceCols]
tmpDend[with(tmpDend,which(tmpDend$mxdgroupnum==36)),c("neuron_name","archive_name")]
tmpWhole[with(tmpWhole,which(tmpWhole$mxdgroupnum==36)),c("neuron_name","archive_name")]
Xname <- "PC1"#ln(Branch_Order_max)"
Yname <- "PC2"#ln(Branch_pathlength_avg)"#"Bif_ampl_remote_avg"#
featurecol <- tmpDend$mxdgroupnum
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
plot(tmpDend[,Xname],tmpDend[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(tmpDend[,Xname],tmpDend[,Yname],labels=as.character(tmpDend$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

clusnums <- unique(reEM$classification)
parms <- reEM$parameters

for(i in 1:length(clusnums)){
  XYpcs <- tmp[with(tmp,which(tmp$classification == clusnums[i])),c(Xname,Yname)]
  XYmeans <- colMeans(XYpcs)
  
  #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
  #       radius=1, center.pch=as.character(clusnums[i]),col=1)
  ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}

unique(subset(wholeReduced[,c("region2","region3","cellclass2","cellclass3")], wholeReduced$neuron_name %in% wholeReduced[with(wholeReduced, 
                                                                     which(wholeReduced$mxdgroupnum == 17)), 
                                                                c("neuron_name")]))

###step#3: check loadings on PC1 vs PC2 for this group###
dotchart(pc1loadings,main="PC1 loadings",cex=0.7, xlab="variable loadings",col=2)
dotchart(pc2loadings,main="PC2 loadings",cex=0.7, xlab="variable loadings",col=2)

###step#4: plot cluster plot on metrics with highest loadings###

###step#5: load data from entire cluster###

###step#6: plot cluster plot on metrics with color code for groups###

###step#7: plot cluster plot on metrics with color code for meta-metrics###

###check single groups -- grp#23###

featurecol <- grp23$classification
colarr <- colorcode(groupnum=featurecol, length(unique(featurecol)))
Xname <- 'PC1'
Yname <- 'PC2'
grpN <- 23
grpttl <- paste(grp23$ct_groupttl[1], grp23$br_groupttl[1], grp23$sp_groupttl[1])
#plotting the original morphometrics
plot(grp23[,Xname],grp23[,Yname],
     col="white",cex=.55,pch=16,xlab=Xname,main=paste(grpN,grpttl),ylab=Yname,cex.main=0.8)
text(grp23[,Xname],grp23[,Yname],labels=as.character(grp23$classification),col=mypalette[colarr],cex=.7)  
legend("bottomleft",text.width = strwidth("1,000"), legend = unique(featurecol), pch=19, cex=0.7, col=unique(mypalette[colarr]))

#overlay ellipsoids on the clusters
clusnums <- unique(grp23$classification)
parms <- reEM$parameters

for(i in 1:length(clusnums)){
  XYpcs <- grp12MetricData[with(grp12MetricData,which(grp12MetricData$classification == clusnums[i])),c(Xname,Yname)]
  XYmeans <- colMeans(XYpcs)
  
  #ellipse(parms$mean[c(Xname,Yname),clusnums[i]],shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
  #       radius=1, center.pch=as.character(clusnums[i]),col=1)
  ellipse(XYmeans,shape=parms$variance$sigma[c(Xname,Yname),c(Xname,Yname),clusnums[i]],
          radius=1, center.pch=as.character(clusnums[i]),col=1)
}