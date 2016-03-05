#### get PC data for dendrites ####
sigdenddf <- DendReduced
rownames(sigdenddf) <- rownames(DendReduced)
#to correlate with the pc space add metric space also for ease in interpretation
dim(Dendritedatamatrix)
sigdendmetric <- Dendritedatamatrix
rownames(sigdendmetric)

#### get PC data for dendrites from overrepresented groups ####
dim(DendReduced)
sigdenddf <- data.frame()
for(i in 1:nrow(sigMatrixDend)){
  sigdenddf <- rbind(sigdenddf, subset(DendReduced, 
                                       DendReduced$mxdgroupnum == sigMatrixDend$row[i] & 
                                         DendReduced$classification == sigMatrixDend$col[i]))
}
#### add pmids to dendrite data ####
sigdenddf <- merge(x=sigdenddf,y=pmidcol,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(sigdenddf)
colnames(sigdenddf)
dim(Dendritedatamatrix)
sigdendmetric <- Dendritedatamatrix
sigdendmetric <- merge(x=sigdendmetric,y=sigdenddf[,c("neuron_name","pmid")],by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(sigdendmetric)
colnames(sigdendmetric)
unique(sigdendmetric$classification)
unique(sigdenddf$classification)

#### original metric dataframe ####
colnames(sigdenddf)
dim(sigdenddf)
colnames(metricDendrites)

origdenddf <- merge(x=metricDendrites[,c(1,33:63)],y=sigdenddf[,c(1:32,60:68)], by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(origdenddf)
colnames(origdenddf)

#flip the -ve loadings on PC2(total length), PC3(contraction) and avg remote amplitude(PC5) for dendrites which is equivalent of size
parmsforplot <- tmpdendEM_95$parameters
str(parmsforplot)
names(parmsforplot)
sigdenddf$PC2 <- -1 * sigdenddf$PC2

#parmsforplot$mean[,2] <- -1 *parmsforplot$mean[,2]
#flip sign along rows
#parmsforplot$variance$sigma[2,,] <- -1 * parmsforplot$variance$sigma[2,,]
#flip sign along column
#parmsforplot$variance$sigma[,2,] <- -1 * parmsforplot$variance$sigma[,2,]
#PC3 contraction is flipped to have tortuosity
sigdenddf$PC3 <- -1 * sigdenddf$PC3
#parmsforplot$PC3 <- -1 *parmsforplot$PC3
#PC4
#sigdenddf$PC4 <- -1 * sigdenddf$PC4
#PC5
#sigdenddf$PC5 <- -1 * sigdenddf$PC5
#parmsforplot$PC5 <- -1 *parmsforplot$PC5
#PC6
sigdenddf$PC6 <- -1 * sigdenddf$PC6
#PC7
#sigdenddf$PC7 <- -1 * sigdenddf$PC7
#### dendrite clusters ####
unique(sigdenddf$classification)
cluster4 <- subset(sigdenddf, sigdenddf$classification == 4)
dim(cluster4)
cluster5 <- subset(sigdenddf, sigdenddf$classification == 5)
dim(cluster5)
cluster6 <- subset(sigdenddf, sigdenddf$classification == 6)
dim(cluster6)
#cluster7 <- subset(sigdenddf, sigdenddf$classification == 7)
#dim(cluster7)
#cluster8 <- subset(sigdenddf, sigdenddf$classification == 8)
#dim(cluster8)
cluster2<- subset(sigdenddf, sigdenddf$classification == 2)
dim(cluster2)
cluster1 <- subset(sigdenddf, sigdenddf$classification == 1)
dim(cluster1)
cluster3 <- subset(sigdenddf, sigdenddf$classification == 3)
dim(cluster3)

#### dendrite clusters for metrics ####
clusterM4 <- subset(sigdendmetric,sigdendmetric$classification == 4)
clusterM6 <- subset(sigdendmetric,sigdendmetric$classification == 6)
clusterM3 <- subset(sigdendmetric,sigdendmetric$classification == 3)
clusterM5 <- subset(sigdendmetric,sigdendmetric$classification == 5)
clusterM2 <- subset(sigdendmetric,sigdendmetric$classification == 2)
clusterM1 <- subset(sigdendmetric,sigdendmetric$classification == 1)
#clusterM7 <- subset(sigdendmetric,sigdendmetric$classification == 7)
#clusterM8 <- subset(sigdendmetric,sigdendmetric$classification == 8)


#plotDensitiesOnCategories(sigdenddf, 'PC1','classification')


#### plot z-scores (mean/SD) for each cluster ####
#colnames(sigdenddf)
#zscrdend <- computezscore(sigdenddf[,c(33:49,67)], clusfeature)
#colnames(zscrdend)
#EucDistFrmGlobalmean <- computezscoreDist(zscrdend)
#compute zscore using mean/clusterSD for each cluster
zscrdend_1 <- computezscore2(sigdenddf[,c(33:49,67)], clusfeature)
zscrdend_1[1,2]
colnames(sigdenddf[,c(33:49,67)])
EucDistFrmGlobalmean_dend <- computezscoreDist(zscrdend_1)

varscr_dend <- computegrpvar(sigdenddf[,c(33:49,67)], clusfeature)
dim(varscr_dend)
varscr_dend[which(zscrdend_1$c==6),'sd']
varscr_dend[which(zscrdend_1$c==6),'mu']
varscr_dend[which(zscrdend_1$c==6),'meanbysd']
varianceMagnitude_dend <- computevarscrMag(varscr_dend)
str(varianceMagnitude_dend)
rownames(varianceMagnitude_dend)[1]

#### plot zscore distances from global mean ####
plotTop <- max(EucDistFrmGlobalmean_dend$EucdistanceFromMu)
mp <- barplot(EucDistFrmGlobalmean_dend$EucdistanceFromMu,cex.main=1,
              ylim=c(0,plotTop))#xlab = "clusters",main="z-score distance from clusters to global mean",
clusteraxis <- zscralphconvers[as.numeric(rownames(EucDistFrmGlobalmean_dend))]
axis(1,at=mp,labels=clusteraxis)

#### plot variances of each cluster for dendrite ####
plotTop <- max(varianceMagnitude_dend$magnitude)
mp <- barplot(varianceMagnitude_dend$magnitude,main="SDs of each cluster",cex.main=1,
              xlab = "clusters",ylim=c(0,plotTop))
axis(1,at=mp,labels=rownames(varianceMagnitude_dend))

#### compute intercluster distances ####
DavidBouldinDist <- interclusterdist_DB(sigdenddf[,c(33:49,67)], clusfeature)
plotdistances(DavidBouldinDist)
write.table(DavidBouldinDist, file="DavidBouldinDist.txt",quote=FALSE, sep=",")

#### distances (z-scores) of each cluster from the mean ####

#plot the normalized Z-scores
#get maximum value on Y-axis from the cluster with highest zscore distance
maxpos <- with(EucDistFrmGlobalmean_dend, which(EucDistFrmGlobalmean_dend$EucdistanceFromMu == max(EucDistFrmGlobalmean_dend$EucdistanceFromMu)))
maxclus <- rownames(EucDistFrmGlobalmean_dend)[maxpos]
ymax <- max(unlist(zscrdend_1[[which(zscrdend_1$c==maxclus),2]])[1:17])
plot(unlist(zscrdend_1[[which(zscrdend_1$c==1),2]])[1:17],type="b",ylim = c(0,ymax),xlim=c(1,18),xlab="principal components",ylab="distances from the data centroid",main="z-scores of each cluster on dendrites")
lines(unlist(zscrdend_1[[which(zscrdend_1$c==2),2]])[1:17],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(zscrdend_1[[which(zscrdend_1$c==3),2]])[1:17],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(zscrdend_1[[which(zscrdend_1$c==4),2]])[1:17],type="b",lty=4,pch=12,col="red")
lines(unlist(zscrdend_1[[which(zscrdend_1$c==5),2]])[1:17],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(zscrdend_1[[which(zscrdend_1$c==6),2]])[1:17],type="b",lty=6,pch=13)#,col="lightgreen")
#lines(unlist(zscrdend_1[[which(zscrdend_1$c==7),2]])[1:17],type="b",lty=7,pch=2)#,col="red")
#lines(unlist(zscrdend_1[[which(zscrdend_1$c==8),2]])[1:17],type="b",lty=7,pch=3)#,col="red")
#lines(unlist(zscrdend_1[[which(zscrdend_1$c==9),2]])[1:17],type="b",lty=7,pch=4)#,col="red")
legend("topright",text.width = strwidth("1"),legend = sort(zscrdend_1$c), pch=c(1,0,6,12,7,13), cex=0.7, col=c(rep(1,3),2,1))

#### plot SDs of each clusters along all features #### 
maxpos <- with(varianceMagnitude_dend, which(varianceMagnitude_dend$magnitude == max(varianceMagnitude_dend$magnitude)))
maxclus <- rownames(varianceMagnitude_dend)[maxpos]
ymax <- max(unlist(varscr_dend[[which(varscr_dend$c==maxclus),2]])[1:17])
plot(unlist(varscr_dend$sd[[1]])[1:17],type="b",cex=0.7, ylim=c(0,ymax),xlim=c(1,17),xlab="principal components",ylab="SDs",main="variances of each cluster")
lines(unlist(varscr_dend$sd[[2]])[1:17],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(varscr_dend$sd[[3]])[1:17],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(varscr_dend$sd[[4]])[1:17],type="b",lty=4,pch=12,col="red")
lines(unlist(varscr_dend$sd[[5]])[1:17],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(varscr_dend$sd[[6]])[1:17],type="b",lty=6,pch=13)#,col="lightgreen")
#lines(unlist(varscr_dend$sd[[7]])[1:17],type="b",lty=7,pch=2)#,col="red")
#lines(unlist(varscr_dend$sd[[8]])[1:17],type="b",lty=7,pch=3)#,col="red")
#lines(unlist(varscr_dend$sd[[9]])[1:17],type="b",lty=7,pch=4)#,col="red")
legend("topright",text.width = strwidth("1"),legend =varscr_dend$c, pch=c(1,0,6,12,7,13), cex=0.7, col=c(rep(1,3),2,1))

#### plot means of all clusters ####
ymax <- max(unlist(varscr_dend[[which(varscr_dend$c==2),3]])[1:18])
ymin <- min(unlist(varscr_dend[[which(varscr_dend$c==5),3]])[1:18])
plot(unlist(varscr_dend$mu[[1]])[1:18],type="b",cex=0.7, ylim=c(ymin,ymax),xlim=c(1,18),xlab="principal components",ylab="means",main="means of each cluster")
lines(unlist(varscr_dend$mu[[2]])[1:18],type="b",lty=2,pch=0)#,col="darkred")
lines(unlist(varscr_dend$mu[[3]])[1:18],type="b",lty=3,pch=6)#,col="darkorange")
lines(unlist(varscr_dend$mu[[4]])[1:18],type="b",lty=4,pch=12,col="red")
lines(unlist(varscr_dend$mu[[5]])[1:18],type="b",lty=5,pch=7)#,col="darkgreen")
lines(unlist(varscr_dend$mu[[6]])[1:18],type="b",lty=6,pch=13)#,col="lightgreen")
#lines(unlist(varscr_dend$mu[[7]])[1:18],type="b",lty=7,pch=2)#,col="red")
legend("topright",text.width = strwidth("1"),legend = varscr_dend$c, pch=c(1,0,6,12,7,13), cex=0.7, col=c(rep(1,3),2,1))

#### plot z-scores from pairwise clusters ####
#plot the normalized and standardized variances
base <- (unlist(zscrdend_1$m[[1]])+unlist(zscrdend_1$m[[2]])+unlist(zscrdend_1$m[[3]])+
           unlist(zscrdend_1$m[[4]])+unlist(zscrdend_1$m[[5]])+unlist(zscrdend_1$m[[6]]))/6
for(i in 1:5){
  for(j in (i+1):6){
    print(i)
    ymax <- max(unlist(zscrdend_1$m[[i]]), unlist(zscrdend_1$m[[j]]))
    plot(base[1:17],type="b",cex=0.7, ylim=c(0,ymax),xlim=c(1,17),col="red",
         xlab="principal components",ylab="zscore",main=paste("zscores between cluster pairs:",varscr_dend$c[i],varscr_dend$c[j]))
    lines(unlist(zscrdend_1$m[[i]])[1:17],type="b",lty=2,pch=zscrdend_1$c[i])#,col="darkred")
    lines(unlist(zscrdend_1$m[[j]])[1:17],type="b",lty=3,pch=zscrdend_1$c[j])#,col="darkorange")
  }
}

#### visualize dendrite clusters ####
largetorqueanlges <-sigdenddf[which(sigdenddf$PC4 < -6),c("neuron_name","archive_name","classification")]
sigdenddf <- subset(sigdenddf, !sigdenddf$neuron_name %in% largetorqueanlges$neuron_name)
dim(sigdenddf)
#highlight all 9 clusters
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdenddf)
         
sigdenddf[which(sigdenddf$PC6 < -3 & sigdenddf$classification == 3),c("neuron_name","archive_name","PC3")]
         
#highlight 4,6,3
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdenddf,c(3,1,5,4,2,6))
PlotDoubleMetricDistrOnClusters('PC3','PC1',sigdenddf,c(6,5,4))
PlotDoubleMetricDistrOnClusters('PC5','PC1',sigdenddf,c(4,5))

PlotDoubleMetricDistrOnClusters('PC5','PC1',sigdenddf,c(2,1))
         
PlotDoubleMetricDistrOnClusters('PC6','PC1',sigdenddf,c(3,6))
         
sigdenddf[which(sigdenddf$PC3 > -1 & sigdenddf$classification == 5),
          c("neuron_name","archive_name","PC3")]
        
PlotDoubleMetricDistrOnClusters('PC7','PC1',sigdenddf,c(3,6))
PlotDoubleMetricDistrOnClusters('PC1','PC5',sigdenddf,c(3,6))
         
sigdenddf[which(sigdenddf$classification == 3 & sigdenddf$PC2 < -5),c('neuron_name','archive_name','PC2')]
sigdenddf[which(sigdenddf$classification == 3 & sigdenddf$PC6 > 2),c('neuron_name','archive_name','PC2')]
sigdenddf[which(sigdenddf$PC1 > 10),c("neuron_name","archive_name")]
sigdenddf[which(sigdenddf$PC1 < -6),c("neuron_name","archive_name")]

PlotDoubleMetricDistrOnClusters('PC5','PC6',sigdenddf,c(2,6))
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdenddf,c(2,1))

PlotDoubleMetricDistrOnClusters('PC2','PC3',sigdenddf)
PlotDoubleMetricDistrOnClusters('PC3','PC4',sigdenddf)
PlotDoubleMetricDistrOnClusters('PC4','PC5',sigdenddf)
PlotDoubleMetricDistrOnClusters('PC5','PC6',sigdenddf)

#### eliminate possible outliers ####
sigdenddf[which(sigdenddf$PC4 < -15),"neuron_name"]
sigdenddf <- sigdenddf[-which(sigdenddf$neuron_name %in% c("13A1-Cell17", "P4-DEV255")),]
sigdenddf[which(sigdenddf$PC5 < -4),"neuron_name"]
unique(sigdenddf$ct_groupttl)

#### plot pairwise/limited #clusters ####
#farthest away clusters

PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdenddf[which(sigdenddf$classification %in% c(6,2)),])
PlotDoubleMetricDistrOnClusters('PC1','PC2',sigdenddf[which(sigdenddf$classification %in% c(2,3)),])
#ganglion cell clusters
highlight = c(3,21)
highlight = c(1,2,21)
alphaind <- alphabetcode(groupnum=sigdenddf$classification, length(sigdenddf$classification))
PlotDoubleMetricDistrOnClustersNGroups('PC1','PC2',sigdenddf[which(sigdenddf$classification %in% c(2,4,5)),],highlight) #,parmsforplot)

#### checking Bif_torque_local_max values ####
#values < 0.33 have daughter bifs sharing the same plane as parent
PlotDoubleMetricDistrOnClustersNGroups('ln(Bif_torque_local_max)','ln(Bif_torque_remote_max)',sigmetricdf[which(sigmetricdf$classification %in% c(2,4,5)),],highlight) #,parmsforplot)
biftorquecoplanar <- sigmetricdf[which(sigmetricdf[,'ln(Bif_torque_local_max)'] < 0.33),"neuron_name"]
           #c("archive_name","neuron_name",'ln(Bif_torque_local_max)',"ln(Depth_total_sum)")]
#check for the original angles
metricDendrites[which(metricDendrites$neuron_name == "Vn01132006-0-C"),]
metricDendrites[which(metricDendrites$neuron_name %in% biftorquecoplanar), "Bif_torque_local_max"]
#the optimal angles are > 20 and < 80 has more natural dendritic shape
biftorquevals <- metricDendrites[which(metricDendrites$Bif_torque_local_max > 20 & metricDendrites$Bif_torque_local_max < 80),c("neuron_name","Bif_torque_local_max","Depth_total_sum","Fractal_Dim_max","cellclass2","archive_name")]
dim(biftorquevals)
write.table(biftorquevals, file="biftorquevals.txt",quote=TRUE, sep=",")
#check for the correspoding natural log values on the log transformed and standardized original values
biftorquevals$neuron_name
sigdendmetric[which(sigdendmetric$neuron_name %in% biftorquevals$neuron_name),]

#### checking Bif_torque_local_avg values ####
#values < 0.33 have daughter bifs sharing the same plane as parent
PlotDoubleMetricDistrOnClustersNGroups('Bif_torque_local_avg','Bif_torque_remote_avg',
                                      sigmetricdf[which(sigmetricdf$classification %in% c(2,4,5)),],highlight) #,parmsforplot)
colnames(sigmetricdf)
biftorquelocaltest <- sigmetricdf[which(sigmetricdf[,'Bif_torque_local_avg'] >0.6),
                                  c("neuron_name","ln(Depth_total_sum)")]
dim(biftorquelocaltest)
#c("archive_name","neuron_name",'ln(Bif_torque_local_max)',"ln(Depth_total_sum)")]
#check for the original angles
metricDendrites[which(metricDendrites$neuron_name == "Vn01132006-0-C"),]
metricDendrites[which(metricDendrites$neuron_name %in% biftorquelocaltest$neuron_name),"Bif_torque_local_avg"]
#the optimal angles are > 20 and < 80 has more natural dendritic shape
biftorqueavgvals <- metricDendrites[which(metricDendrites$Bif_torque_local_avg > 20 & 
                                            metricDendrites$Bif_torque_local_avg < 80),c("neuron_name","Bif_torque_local_avg","Fractal_Dim_max","archive_name")]
dim(biftorqueavgvals)
write.table(biftorqueavgvals, file="biftorqueavgvals.txt",quote=TRUE, sep=",")


#### adding cmpoundGrps & cmpoundGrpName to sigdenddf and origdenddf####

sigdenddf <- sigdenddf[, -which(names(sigdenddf) %in% c("cmpoundGrps","cmpoundGrpName"))]
sigdendmetric <- sigdendmetric[, -which(names(sigdendmetric) %in% c("cmpoundGrps","cmpoundGrpName"))]
origdenddf <- origdenddf[, -which(names(origdenddf) %in% c("cmpoundGrps","cmpoundGrpName"))]

colnames(sigdenddf)
colnames(origdenddf)
tretinalgroups <- unique(sigdenddf[which(sigdenddf$region1 == "Retina" & !sigdenddf$species == "Salamander"),"mxdgroupnum"])
#### add cmpoundGrps and cmpoundGrpName to the original data frame ####
#cluster3 rodent retinal ganglion cells
cond <- list(classification = 3, region1 = "Retina", order = "Rodent")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
dim(tmpgrp)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 1
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Rodent Ganglion"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 1
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Rodent Ganglion"

cond <- list(classification = 3, order="Insects")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 2
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Fly Tangential"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 2
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Fly Tangential"

cond <- list(classification = 3,archive_name="Smit-Rigter")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 3
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "young rat pyramidal"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 3
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Young rat pyramidal"

cond <- list(classification = 6)
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 4
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Human basal pyramidal"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 4
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Human basal pyramidal"

#scatter plot of rodent ganglion cells with high branching density showing their space filling function
#jpeg(filename='spacefilling_1.jpg',width=8.5,height=8,units='cm',res=600,pointsize=0.7)
x11()
PlotDoubleMetricDistrOnClustersNGroups('PC6','PC1','asymmetry index','branching density',
                                       sigdenddf,c(1,4,2))
#histogram showing the difference in branching density between rodent and jacobs basal dendrites
colnames(origdenddf)
colnames(sigdendmetric)
dim(origdenddf)
####fig6####
x11()
par(mfrow=c(1,3))
plotmetricHist(origdenddf,"Bif_ampl_remote_max", "Maximum remote amplitude angle (in degrees)",c(1,4),binsize=20,legendpos='topleft')
plotmetricHist(origdenddf,"Branch_pathlength_max","Maximum branch pathlength (mu)",logT=TRUE,c(1,4),binsize=30)#,scaleupby=c(2,2))

plotmetricHist(origdenddf, "Partition_asymmetry_avg","Average asymmetry index",c(1,2),round2=2,unts='microns in log-scale',scaleupby=c(2,6),binsize=20,legendpos='topleft')


#fish and salmander ganglion cells
cond <- list(classification = c(1,3),order="Bony fishes")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 5
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Fish Ganglion"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 5
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Fish Ganglion"

cond <- list(classification = 1,order="Amphibians")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 6
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Salamander Ganglion"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 6
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Salamander Ganglion"

x11()
PlotDoubleMetricDistrOnClustersNGroups('PC2','PC3','size','tortuosity',
                                       sigdenddf,c(1,5,6))

plotmetricHist(origdenddf,"Length_total_sum",logT=TRUE,c(1,5,6),unts='microns in log-scale',scaleupby=c(4,4,4))
#reverse the contraction values to show tortuosity
origdenddf[,'Contraction_avg'] <- 1/origdenddf[,'Contraction_avg']
names(origdenddf)[16] <- 'Tortuosity_avg'
plotmetricHist(origdenddf,"Tortuosity_avg",logT=TRUE,c(1,5,6),round=2,unts='microns in log-scale',scaleupby=c(8,8,8))


#pyramidal cells from cluster#5
cond <- list(classification = 5,cellclass2="Pyramidal cell",archive_name=c("Lewis","Wearne_Hof"))#"Luebke","Smit-Rigter","Jacobs",
#cond <- list(classification = 5)
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
dim(tmpgrp)
unique(tmpgrp$archive_name)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 7
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Incomplete Neocortical pyramidal"
dim(sigdenddf[which(sigdenddf$cmpoundGrps == 7),])
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 7
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Incomplete Neocortical pyramidal"

cond <- list(classification = 5,cellclass2="Granule cell")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 8
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Dentate gyrus granule cell"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 8 
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Dentate gyrus granule cell"

cond <- list(classification = 4,ct_groupttl="Pyramidal")
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
dim(tmpgrp)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 9
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Rodent Neocortical pyramidal"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 9
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "Rodent Neocortical pyramidal"


#PlotDoubleMetricDistrOnClustersNGroups('PC5','PC1','tortuosity','branching density',
#                                       sigdenddf,c(1,8,7))
###fig7 plot####
x11()
PlotDoubleMetricDistrOnClustersNGroups('PC3','PC1','tortuosity','branching density',
                                       sigdenddf,c(7,9))
points(c(-5.134872,-4.254340), c(-5.169661,-3.166460),col = "red",pch=8)
points(c(-1.1239570,-1.1658539),c(-0.43564123,-0.96729958),col = "blue",pch=8)
points(c(0.17238324,0.97776150),c( 0.99033299,0.96358669),col="black",pch=8)
points(c(5.29684542,5.107149),c(5.275715,4.130621),col="dark green",pch=8)
plot(sigdenddf[which(sigdenddf$cmpoundGrps == 9),c('PC1','PC3')])

plot(sigdenddf[which(sigdenddf$cmpoundGrps == 7),c('PC1','PC3')])

paste(sigdenddf[which(sigdenddf$cmpoundGrps == 7 & sigdenddf$PC3<=-3 & sigdenddf$PC1 <=-3),c("neuron_name","PC1","PC3")],collapse=",")

paste(sigdenddf[which(sigdenddf$cmpoundGrps == 7 & sigdenddf$PC3<=0 & sigdenddf$PC3 >= -1),c("neuron_name","PC1","PC3")],collapse=",")

paste(sigdenddf[which(sigdenddf$cmpoundGrps == 9 & sigdenddf$PC1>=0 & sigdenddf$PC1 <= 1 & sigdenddf$PC3>=0 & sigdenddf$PC3<=1),c("neuron_name","PC1","PC3")],collapse=",")

paste(sigdenddf[which(sigdenddf$cmpoundGrps == 9 & sigdenddf$PC3>4 & sigdenddf$PC1>4),c("neuron_name","PC1","PC3")],collapse=",")


PlotDoubleMetricDistrOnClustersNGroups('PC2','PC1','size','branching density',
                                       sigdenddf,c(7,9,8,1))

PlotDoubleMetricDistrOnClustersNGroups("Bif_ampl_remote_max","Fractal_Dim_max","branching density","tortuosity",
                                       origdenddf,c(7,9,8))
#test scatter
plot(sigmetricdf$Bif_ampl_remote_max,sigmetricdf[,'ln(Fractal_Dim_max)'])

#Bif_ampl_remote_max
Branch_pathlength_max
x11()
plotmetricHist(origdf=origdenddf,metric="Branch_pathlength_avg",groups=c(1,9),
               unts='in degrees',scaleupby=c(2,2,2,2),legendpos='topright')
plotmetricHist(origdf=origdenddf,metric="Fractal_Dim_max",groups=c(7,9,8),
               logT=TRUE,round2=2,scaleupby=c(2,2,2),legendpos='topright')




dim(sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$archive_name == "Wearne_Hof" & sigdenddf$age_class == "Old"),c("neuron_name","archive_name","age_class")])
colnames(sigmetricdf)
sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$archive_name %in% c("Luebke","Smit-Rigter") & sigdenddf$PC1 > 0),c("neuron_name","archive_name","age_class","min_age","max_age","age_scale")]
X11()
plotmetricHist(origdenddf,"Branch_pathlength_max",logT=TRUE,c(1,8),unts='microns in log-scale',scaleupby=c(3,4),legendpos='topleft')
plotmetricHist(origdenddf,"Branch_pathlength_avg",logT=TRUE,c(1,7),unts='microns in log-scale',scaleupby=c(3,3))
plotmetricHist(origdenddf,"Bif_ampl_local_max",c(1,8,7),unts='in degrees',scaleupby=c(3,3,3),legendpos='topleft')

plotmetricHist(origd=origdenddf,metric="Tortuosity_avg",unts='microns in log-scale',logT=TRUE,groups=c(1,8),round=2,scaleupby=c(3,3))
plotmetricHist(origd=origdenddf,metric="Branch_pathlength_max",logT=TRUE,groups=c(1,7),round=2,scaleupby=c(3,3))
plotmetricHist(origd=origdenddf,metric="Bif_ampl_local_max",groups=c(1,7),round=2,scaleupby=c(3,3))

logTransformData(origdenddf$Bif_tilt_remote_max)
plotmetricHist(origdf=origdenddf,metric="Bif_tilt_remote_avg",groups=c(1,8,7),unts='in degrees',scaleupby=c(2,2,2),legendpos='topleft')
plotmetricHist(origdf=origdenddf,metric="Bif_tilt_local_avg",groups=c(1,8,7),unts='in degrees',scaleupby=c(2,2,2),legendpos='topleft')

plotmetricHist(origdf=origdenddf,metric="Bif_ampl_remote_max",groups=c(1,8,7),unts='in degrees',scaleupby=c(2,2,2),legendpos='topleft')
plotmetricHist(origdf=origdenddf,metric="Fractal_Dim_max",groups=c(1,8,7),logT=TRUE,round2=2,unts='in degrees',scaleupby=c(2,2,2),legendpos='topleft')

paste(origdenddf[which(origdenddf$Bif_tilt_remote_avg>145),"neuron_name"],collapse=",")
paste(sample(origdenddf[which(origdenddf$cmpoundGrps == 4),"neuron_name"],5),collapse=",")


sigdenddf[which(origdenddf$classification == 5 & origdenddf$Bif_tilt_remote_max > 45),c("neuron_name","archive_name")]
clus5pyr <- origdenddf[which(origdenddf$cmpoundGrps == 8 & origdenddf$Bif_tilt_remote_avg <145) ,c("neuron_name","Bif_tilt_remote_avg","Bif_ampl_remote_avg")]
clus3ret <- origdenddf[which(origdenddf$cmpoundGrps == 1 & origdenddf$Bif_tilt_remote_avg <145) ,c("neuron_name","Bif_tilt_remote_avg","Bif_ampl_remote_avg")]
wearnehof <- origdenddf[which(origdenddf$archive_name == "Wearne_Hof"),c("neuron_name","Bif_tilt_remote_avg","Fractal_Dim_avg","Fractal_Dim_max","Tortuosity_avg","Bif_tilt_local_avg")]
jacobs <- origdenddf[which(origdenddf$archive_name == "Jacobs" & origdenddf$cmpoundGrps == 8),c("neuron_name","Bif_tilt_remote_avg","Fractal_Dim_avg","Fractal_Dim_max","Tortuosity_avg","Bif_tilt_local_avg")]
range(wearnehof$Bif_tilt_remote_avg)
range(jacobs$Bif_tilt_remote_avg)
range(wearnehof$Fractal_Dim_avg)
range(jacobs$Fractal_Dim_avg)
range(wearnehof$Tortuosity_avg)
range(jacobs$Tortuosity_avg)
range(wearnehof$Bif_tilt_local_avg)
range(jacobs$Bif_tilt_local_avg)
sample()
dim(wearnehof)

plot(c(wearnehof$Fractal_Dim_max, jacobs$Fractal_Dim_max),c(wearnehof$Tortuosity_avg, jacobs$Tortuosity_avg))



range(clus5pyr$Bif_ampl_remote_avg)
range(clus5pyr$Bif_tilt_remote_avg)
range(clus3ret$Bif_ampl_remote_avg)
range(clus3ret$Bif_tilt_remote_avg)

sigdenddf[which(sigdenddf$classification == 6 & sigdenddf$lobe == "Insula"),"gender"]
cond <- list(classification = 2)
tmpgrp <- getSubset2(sigdenddf,specifications=cond,typeMetaData=sigdenddf)
dim(tmpgrp)
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 9
sigdenddf[which(sigdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "cluster#2"
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrps"] <- 9
origdenddf[which(origdenddf$neuron_name %in% tmpgrp$neuron_name),"cmpoundGrpName"] <- "cluster#2"

PlotDoubleMetricDistrOnClustersNGroups('PC3','PC1','tortuosity','branching density',
                                       sigdenddf,c(1,8,7,9))
plotmetricHist(origd=origdenddf,metric="Tortuosity_avg",logT=TRUE,groups=c(1,9),round=2,scaleupby=c(3,3))
plotmetricHist(origd=origdenddf,metric="Branch_pathlength_max",logT=TRUE,groups=c(1,9),round=2,scaleupby=c(3,3))
plotmetricHist(origd=origdenddf,metric="Bif_ampl_local_max",groups=c(1,9),round=2,scaleupby=c(3,3))

#add same info for metric data too
sigdendmetric[which(sigdendmetric$classification  == 3 & sigdendmetric$region1 == "Retina") ,"cmpoundGrps"] <- 1
sigdendmetric[which(sigdendmetric$classification  == 3 & sigdendmetric$region1 == "Retina") ,"cmpoundGrpName"] <- "Rodent retinal ganglion"
#add new group name to original data
origdenddf[which(origdenddf$classification  == 3 & origdenddf$region1 == "Retina") ,"cmpoundGrps"] <- 1
origdenddf[which(origdenddf$classification  == 3 & origdenddf$region1 == "Retina") ,"cmpoundGrpName"] <- "Rodent retinal ganglion"

#add new group for jacobs basal dendrites
Hbasalgroups <- unique(sigdenddf[which(sigdenddf$classification == 6),"mxdgroupnum"])
sigdenddf[which(sigdenddf$classification == 6),"cmpoundGrps"] <- 2
sigdenddf[which(sigdenddf$classification == 6),"cmpoundGrpName"] <- "jacobs basal"
sigdendmetric[which(sigdendmetric$classification == 6),"cmpoundGrps"] <- 2
sigdendmetric[which(sigdendmetric$classification == 6),"cmpoundGrpName"] <- "jacobs basal"
origdenddf[which(origdenddf$classification  == 6) ,"cmpoundGrps"] <- 2
origdenddf[which(origdenddf$classification  == 6) ,"cmpoundGrpName"] <- "jacobs basal"

#blowfly
sigdenddf[which(sigdenddf$classification  == 3 & sigdenddf$cellclass2 =="Tangential cell"),"cmpoundGrps"] <- 3
sigdenddf[which(sigdenddf$classification == 3 & sigdenddf$cellclass2 =="Tangential cell"), "cmpoundGrpName"] <- "Tangential"
origdenddf[which(origdenddf$classification  == 3 & origdenddf$cellclass2 == "Tangential cell"),"cmpoundGrps"] <- 3
origdenddf[which(origdenddf$classification  == 3  & origdenddf$cellclass2 == "Tangential cell"),"cmpoundGrpName"] <- "Tangential"

#smit-ritger
         
unique(sigdenddf$region1)
         colnames(sigdenddf)
sigdenddf[which(sigdenddf$classification  %in% c(3,5) & sigdenddf$archive_name =="Smit-Rigter"),c("neuron_name","age_class","classification")]
sigdenddf[which(sigdenddf$classification  == 3 & sigdenddf$archive_name =="Smit-Rigter"),"cmpoundGrps"] <- 4
  
sigdenddf[which(sigdenddf$classification == 3 & sigdenddf$archive_name =="Smit-Rigter"), "cmpoundGrpName"] <- "Pyramidal"
origdenddf[which(origdenddf$classification  == 3 & origdenddf$archive_name =="Smit-Rigter"),"cmpoundGrps"] <- 4
origdenddf[which(origdenddf$classification  == 3  & origdenddf$archive_name =="Smit-Rigter"),"cmpoundGrpName"] <- "Pyramidal"
       
#add fish ganglion cell class
unique(sigdenddf$order)
sigdenddf[which(sigdenddf$classification  %in% c(3,1) & sigdenddf$order == "Bony fishes" ),"cmpoundGrps"] <- 5
sigdenddf[which(sigdenddf$classification %in% c(3,1)  & sigdenddf$order == "Bony fishes" ), "cmpoundGrpName"] <- "Fish"
#sigdendmetric[which(sigdendmetric$classification  == 3 & sigdendmetric$order == "Bony fishes"),"cmpoundGrps"] <- 5
#sigdendmetric[which(sigdendmetric$classification  == 3  & sigdendmetric$order == "Bony fishes"),"cmpoundGrpName"] <- "Fish only"
origdenddf[which(origdenddf$classification %in% c(3,1) & origdenddf$order == "Bony fishes"),"cmpoundGrps"] <- 5
origdenddf[which(origdenddf$classification %in% c(3,1) & origdenddf$order == "Bony fishes"),"cmpoundGrpName"] <- "Fish"

#salamander cells only
sigdenddf[which(sigdenddf$order == "Amphibians" ),"cmpoundGrps"] <- 6
sigdenddf[which(sigdenddf$order == "Amphibians" ), "cmpoundGrpName"] <- "Salamander"
origdenddf[which(origdenddf$order == "Amphibians"),"cmpoundGrps"] <- 6
origdenddf[which(origdenddf$order == "Amphibians"),"cmpoundGrpName"] <- "Salamander"

#cluster 5 
sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$cellclass2 == "Pyramidal cell"),"cmpoundGrps"] <- 7
sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$cellclass2 == "Pyramidal cell"), "cmpoundGrpName"] <- "pyramidal"     
origdenddf[which(origdenddf$classification == 5 & origdenddf$cellclass2 == "Pyramidal cell"),"cmpoundGrps"] <- 7
origdenddf[which(origdenddf$classification == 5 & origdenddf$cellclass2 == "Pyramidal cell"),"cmpoundGrpName"] <- "Pyramidal"

         
#choose random sample of 5 examples from each group

clus1Npyr <- cluster1[sample(1:nrow(cluster1),5,replace=FALSE),"neuron_name"]
paste(clus1Npyr,collapse=',')
unique(cluster5$species)         
clusSal <- subset(cluster1,cluster1$species == "Salamander")
dim(clusSal) 
         
salsample <- clusSal[sample(1:nrow(clusSal),5,replace=FALSE),"neuron_name"]
paste(salsample,collapse=',')
         
clus3ret <- cluster3[sample(1:nrow(cluster3),5,replace=FALSE),"neuron_name"]
paste(clus3ret,collapse=',')
clusfish <- subset(sigdenddf, sigdenddf$cmpoundGrps == 5)
fishsample <- clusfish[sample(1:nrow(clusfish),5,replace=FALSE),"neuron_name"]
paste(fishsample,collapse=',')
         
blwfly <- subset(sigdenddf, sigdenddf$cmpoundGrps == 3)
blwflysample <- blwfly[sample(1:nrow(blwfly),5,replace=FALSE),"neuron_name"]
paste(blwflysample,collapse=',')    

jacobsbasal <- subset(sigdenddf, sigdenddf$cmpoundGrps == 2)
jacobsbasalsample <- jacobsbasal[sample(1:nrow(jacobsbasal),5,replace=FALSE),"neuron_name"]
paste(jacobsbasalsample,collapse=',')    


clus5pyr <- cluster5[sample(1:nrow(cluster5),5,replace=FALSE),"neuron_name"]
paste(clus5pyr,collapse=',')

cluspyronly <- subset(sigdenddf,sigdenddf$cmpoundGrps %in% c(7,9))
dim(cluspyronly)
pyronly <-  sample(cluspyronly$neuron_name,5,replace=FALSE)
paste(pyronly,collapse=',')
dim(clusSal) 
         
#cluster 5 dentate gyrus granule cell
sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$archive_name == "Claiborne"),"cmpoundGrps"] <- 8
sigdenddf[which(sigdenddf$classification == 5 & sigdenddf$archive_name == "Claiborne"), "cmpoundGrpName"] <- "DG granule"     
origdenddf[which(origdenddf$classification == 5 & origdenddf$archive_name == "Claiborne"),"cmpoundGrps"] <- 8
origdenddf[which(origdenddf$classification == 5 & origdenddf$archive_name == "Claiborne"),"cmpoundGrpName"] <- "DG granule"
         
#other cluster3 cells
sigdenddf[which(sigdenddf$classification  == 3 & !sigdenddf$cmpoundGrps %in% c(4,5,6,7)),"cmpoundGrps"] <- 8
sigdenddf[which(sigdenddf$classification == 3 & !sigdenddf$cmpoundGrps %in% c(4,5,6,7)), "cmpoundGrpName"] <- "Other flat cells"
origdenddf[which(origdenddf$classification  == 3 &  !origdenddf$cmpoundGrps %in% c(4,5,6,7)),"cmpoundGrps"] <- 8
origdenddf[which(origdenddf$classification  == 3  & !origdenddf$cmpoundGrps %in% c(4,5,6,7)),"cmpoundGrpName"] <- "Other flat cells"

unique(sigdenddf$cmpoundGrps)
unique(origdenddf[,c('cmpoundGrps','cmpoundGrpName')])
         
#unique species in retinal cells
unique(cluster3$species)
unique(cluster7$species)
unique(cluster8$species)
unique(cluster1$species)
unique(cluster2$species)
unique(cluster5$species)





sigdenddf[which(sigdenddf$mxdgroupnum %in% tretinalgroups),c("neuron_name","note","cellclass3")]

sigdenddf[which(sigdenddf$mxdgroupnum == 21 & sigdenddf$classification == 2),c("neuron_name","note","cellclass3")]

sigdenddf[which(sigdenddf$mxdgroupnum == 1 & sigdenddf$classification == 5),c("neuron_name","note", "cellclass3")]
sigdenddf[which(sigdenddf$mxdgroupnum == 1 & sigdenddf$classification == 2),c("neuron_name","note", "cellclass3")]


unique(sigdenddf[which(sigdenddf$mxdgroupnum == 21),c("archive_name","cellclass2","cellclass3","pmid")])
unique(sigdenddf[which(sigdenddf$mxdgroupnum == 2),c("archive_name","cellclass2","cellclass3","pmid")])
unique(sigdenddf[which(sigdenddf$mxdgroupnum == 36),c("archive_name","cellclass2","cellclass3","pmid")])
unique(sigdenddf[which(sigdenddf$mxdgroupnum == 22),c("archive_name","cellclass2","cellclass3","pmid")])
unique(sigdenddf[which(sigdenddf$mxdgroupnum == 1),c("archive_name","cellclass2","cellclass3","pmid")])

ganglioncells <- subset(sigdenddf, sigdenddf$cellclass2 == "Ganglion cell")
dim(ganglioncells)
tmp_gang <- Mclust(ganglioncells[,33:50],G=1:15)
summary(tmp_gang)

plot(tmp_gang, what="classification",dimens=1:2)
ganglioncells$classification = tmp_gang$classification
PlotDoubleMetricDistrOnClustersNGroups('PC1','PC2',ganglioncells)
PlotDoubleMetricDistrOnClustersNGroups('PC1','PC2',ganglioncells,c(1,2))
PlotDoubleMetricDistrOnClustersNGroups('PC1','PC2',sigdenddf[which(sigdenddf$classification %in% c(5,2)),],c(1,2,32,30,31,19,23))
#### cluster 4 cortical interneurons ####
rodentINsFrm4 <- cluster4[which(cluster4$cellclass1 == "Interneuron" & cluster4$species != "Blowfly"),]
dim(rodentINsFrm4)
unique(rodentINsFrm4$mxdgroupnum)
PlotDoubleMetricDistrOnClustersNGroups('PC1','PC2',sigdenddf[which(sigdenddf$classification %in% c(1,2,4)),],c(21,3,24))
#highlight the retinal and human basal groups in a scatter plot
PlotDoubleMetricDistrOnClustersNGroups('PC4','PC2',sigdenddf,cond,c(1,3,5,6,4))
PlotDoubleMetricDistrOnClustersNGroups('PC3','PC1',sigdenddf,c(1,7,8))
PlotDoubleMetricDistrOnClustersNGroups('PC2','PC4',sigdenddf,c(1,7,8))
#highlight the retinal and human basal dendrites on metrics
colnames(sigmetricdf)
dim(metricDendrites)
PlotDoubleMetricDistrOnClustersNGroups('Bif_tilt_local_avg','Bif_ampl_local_max',origdenddf,c(1,2))
PlotDoubleMetricDistrOnClustersNGroups('Bif_ampl_remote_max','Contraction_avg',origdenddf,c(1,2))
PlotDoubleMetricDistrOnClustersNGroups('Bif_ampl_remote_max','Bif_torque_local_avg',origdenddf,c(1,2))

colnames(sigdenddf)
#### plot histograms per each metric on groups ####

colnames(sigdendmetric)
colnames(metricDendrites)
metric <- "Partition_asymmetry_avg" #"ln(Depth_total_sum)"#"ln(Contraction_avg)"#"Bif_ampl_remote_max"
origmetric <- "Partition_asymmetry_avg" #"Depth_total_sum" #"Contraction_avg"
xmin <- min(sigdendmetric[,metric])
xmax <- max(sigdendmetric[,metric])
g1 <- subset(sigdendmetric, sigdendmetric$cmpoundGrps==1)
g2 <- subset(sigdendmetric, sigdendmetric$cmpoundGrps==2)
g1[,metric]
g2[,metric]
clusterM1[,metric]
mean(metricDendrites[which(metricDendrites$neuron_name %in% g1$neuron_name),origmetric])
mean(metricDendrites[which(metricDendrites$neuron_name %in% g2$neuron_name),origmetric])

hist(origdenddf[which(origdenddf$cmpoundGrps == 1),metric], xlim=c(xmin,xmax), xlab=metric, 
     breaks=30, col="grey",cex.main=0.8, main="group distributions")
         
axpos <- seq(from=min(origdenddf[,metric]),to=max(origdenddf[,metric]),by=diff(range(origdenddf[,metric]))/8)
print(axpos)
axval <- axpos
exp(axpos)
secondval <- 
hist(origdenddf[,metric], xaxt='n', xlab=metric, col="grey",cex.main=0.8, main="group distributions")
axis(side=1,at=axpos,labels=round(axval,0),cex=0.7)
axis(side=4,at=origdenddf[,metric]/100,labels=origdenddf[,metric]/100, cex=0.7)         
hist(origdenddf[which(origdenddf$cmpoundGrps == 4),metric],add = T, breaks=30, col=rgb(0, 1, 0, 0.5))
hist(origdenddf[which(origdenddf$cmpoundGrps == 5),metric],add = T, breaks=30, col=rgb(0, 0, 1, 0.5))
hist(origdenddf[which(origdenddf$cmpoundGrps == 6),metric],add = T, breaks=30, col=rgb(1, 0, 0, 0.5))
hist(origdenddf[which(origdenddf$cmpoundGrps == 7),metric],add = T, breaks=30, col=rgb(0.3, 0.2, 1, 0.5))
hist(origdenddf[which(origdenddf$cmpoundGrps == 8),metric],add = T, breaks=30, col=rgb(1, 0.3, 0.2, 0.5))


legend('topleft',c('retina','jacobs'),
       fill = c(rgb(0, 1, 0, 0.5),rgb(0, 0, 1, 0.5)), bty = 'n',
       border = NA)
dim(origdenddf)
colnames(origdenddf)
colnames(sigdendmetric)
#plot histogram on "Bif_ampl_remote_max" highlighting the two cmpoundgrps (1,2)
origdenddf[which(origdenddf$cmpoundGrps == 1),"cmpoundGrpName"]
origdenddf[,"Bif_ampl_remote_max"]
         
max(origdenddf[which(origdenddf$cmpoundGrps == 1),'Branch_pathlength_avg'])         
plotmetricHist(origdenddf,"Bif_ampl_remote_max",c(1,2),"degrees",scaleupby=c(2,2),legend='topleft')
plotmetricHist(origdenddf,"Bif_ampl_remote_max",c(1,5,6),"degrees",scaleupby=c(2,10,10),legend='topleft')
plotmetricHist(origdenddf, "Branch_pathlength_avg", c(1,2),"in microns",scaleupby=c(4,4))  
plotmetricHist(origdenddf, "Branch_pathlength_avg", c(1,5,6),"in microns",scaleupby=c(4,50,50))

plotmetricHist(origdenddf, "Bif_torque_remote_avg",c(1,7),"degrees")
plotmetricHist(origdenddf, "Fractal_Dim_avg",c(1,2),logT=FALSE,round=2,scaleupby=c(10,10))
         

plotmetricHist(origdenddf, "Depth_total_sum", c(1,3,4),logT=TRUE,"microns in log scale")

plotmetricHist(origdenddf, "Bif_tilt_local_avg",c(1,7),"degrees")
plotmetricHist(origdenddf, "Contraction_avg", c(1,2),logT=FALSE,round2=2,scaleupby=c(10,10))
         

plotmetricHist(origdenddf,"Partition_asymmetry_avg", c(1,7,2),"microns in log scale",round2=2)

plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,2,3,4)),
                          "Bif_ampl_remote_max","cmpoundGrps","cmpoundGrpName","in degrees")
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,2,7)),
                                   "Length_total_sum","cmpoundGrps","cmpoundGrpName","in microns")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,2,7)),
                                   "Bif_ampl_remote_max","cmpoundGrps","cmpoundGrpName","in degrees")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7)),
                                   "Branch_pathlength_avg","cmpoundGrps","cmpoundGrpName","in microns")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7)),
                                   "Bif_tilt_local_avg","cmpoundGrps","cmpoundGrpName","in degrees")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7)),
                                   "Bif_torque_local_avg","cmpoundGrps","cmpoundGrpName","in degrees")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7)),
                                   "Bif_ampl_local_avg","cmpoundGrps","cmpoundGrpName","in degrees")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7)),
                                   "Contraction_avg","cmpoundGrps","cmpoundGrpName")
         
plotDensitiesOnCategories(subset(origdenddf,origdenddf$cmpoundGrps %in% c(1,7,2)),
                                   "Partition_asymmetry_avg","cmpoundGrps","cmpoundGrpName")
origdenddf[which(origdenddf$)]
#### cluster 4 SS pyramidal cells
