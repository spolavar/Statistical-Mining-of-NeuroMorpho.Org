###plot time vs. newly added terms###
spfile <- "C:/Users/sridevi/thesismaterial/ontosearchpaper/speciesrate.txt"
brfile <- "C:/Users/sridevi/thesismaterial/ontosearchpaper/regionrate.txt"
csfile <- "C:/Users/sridevi/thesismaterial/ontosearchpaper/cellclassrate.txt"

metadatasp <- read.csv(spfile,header = FALSE, sep="\t",stringsAsFactors=FALSE)
metadatabr <- read.csv(brfile,header = FALSE, sep="\t",stringsAsFactors=FALSE)
metadatacs <- read.csv(csfile,header = FALSE, sep="\t",stringsAsFactors=FALSE)

print("finished reading...")
#add header
metadatasp <- metadatasp[,1:3]
names(metadatasp) <- c("species","date","count")
metadatabr <- metadatabr[,1:3]
names(metadatabr) <- c("regions","date","count")
metadatacs <- metadatacs[,1:3]
names(metadatacs) <- c("class","date","count")

dim(metadatasp)
dim(metadatabr)
dim(metadatacs)


#plot date vs cumulative count of reconstructions
library(ggplot2)
library(plyr)

sp <- ddply(metadatasp,.(date,species),summarize, 
            cum = cumsum(count)) 

spcloud <- aggregate(metadatasp[,3],by=list(metadatasp[,"species"]),FUN=sum)
brcloud <- aggregate(metadatabr[,3],by=list(metadatabr[,"regions"]),FUN=sum)
brcloud <- brcloud[with(brcloud, order(-x,Group.1)),]
cscloud <- aggregate(metadatacs[,3],by=list(metadatacs[,"class"]),FUN=sum)
cscloud <- cscloud[with(cscloud, order(-x,Group.1)),]

cloud <- rbind(spcloud,brcloud,cscloud)

#add a cumulative column for species
sp1 <- aggregate(metadatasp[,1],by=list(metadatasp[,"date"]),FUN=length)
names(sp1)=c("date","concepts")
sp2 <- aggregate(metadatasp[,3],by=list(metadatasp[,"date"]),FUN=sum)
names(sp2)=c("date","totalsum")
sp <- merge(sp1,sp2,"date")
sp <- sp[with(sp, order(as.Date(date), concepts,totalsum)),]
sp <- within(sp, cum <- cumsum(totalsum))
colnames(sp)

#add a cumulative column for brain regions
br1 <- aggregate(metadatabr[,1],by=list(metadatabr[,"date"]),FUN=length)
names(br1)=c("date","concepts")
br2 <- aggregate(metadatabr[,3],by=list(metadatabr[,"date"]),FUN=sum)
names(br2)=c("date","totalsum")
br <- merge(br1,br2,"date")
br <- br[with(br, order(as.Date(date), concepts,totalsum)),]
br <- within(br, cum <- cumsum(totalsum))
colnames(br)


#add a cumulative column for cell classes
cs1 <- aggregate(metadatacs[,1],by=list(metadatacs[,"date"]),FUN=length)
names(cs1)=c("date","concepts")
cs2 <- aggregate(metadatacs[,3],by=list(metadatacs[,"date"]),FUN=sum)
names(cs2)=c("date","totalsum")
cs <- merge(cs1,cs2,"date")
cs <- cs[with(cs, order(date, concepts,totalsum)),]
cs <- within(cs, cum <- cumsum(totalsum))
colnames(cs)

counts <- data.frame("species" = sp$cum, "brain regions" = br$cum, "cell class" = cs$cum)
matrix <- do.call(rbind, counts)
colnames(matrix) <- sp$date

concepts <- data.frame("species" = sp$concepts, "brain regions" = br$concepts, "cell class" = cs$concepts)
matrixc <- do.call(rbind,concepts)
colnames(matrixc) <- sp$date


par(new=TRUE)
#plot matrix on three groups
p <- barplot(matrix, main="growth rate of main axes",
        xlab="time", ylab = "#reconstructions", col=c("darkblue","red", "green"),
        legend = rownames(matrix), beside=TRUE,  args.legend = list(x = "topleft", bty="n"))
#add text on top of the bars
text(x = p, y = matrix[1,], label = matrixc[1,],pos=3,cex=0.8,col="purple")

library(dplyr)
library(wordcloud)

wordcloud(spcloud[,1],spcloud$x, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
wordcloud(brcloud[,1],brcloud$x, random.order=FALSE, colors=brewer.pal(8, "Dark2"))
wordcloud(cloud[,1],cloud$x, random.order=FALSE, colors=brewer.pal(8, "Dark2"))



#require(graphics)
neuron_name <- metricDataWhole[,1]
neuron_name
dim(metricDataWhole)
LMmetricsWhole <- metricDataWhole[,27:140]
#LMmetricsWhole <- cbind(neuron_name,LMmetricsWhole)
dim(LMmetricsWhole)

colnames(LMmetricsWhole)
LMmetricsWhole <- ClnMrgName(LMmetricsWhole)
rownames(LMmetricsWhole)
LMmetricsWhole$neuron_name

ncol(LMmetricsWhole)
#eliminate variables that were manually sorted out as similar to the existing ones
LMmetricsWhole <- LMmetricsWhole[, !(colnames(LMmetricsWhole) %in% c("Fragmentation_min","Fractal_Dim_min","Taper_1_min","Taper_2_min","Contraction_max","Daughter_Ratio_min"))]
ncol(LMmetricsWhole)
chkNegVals(LMmetricsWhole[2:115])
#if mean < 0, then remove that column from CV
chkMeanCols(LMmetricsWhole[2:115])
dim(LMmetricsWhole)
cvarr <- cvfun(LMmetricsWhole)

cvframe <- data.frame("metrics"=colnames(LMmetricsWhole), "coefficient of variation" = cvarr)
cvframe[order(cvframe$coefficient.of.variation,decreasing=TRUE),]
dim(cvframe)

#add neuron_name column to LMmetricsWhole
#LMmetricsWhole <- cbind(neuron_name,LMmetricsWhole)
#neuron_name_frm <- data.frame(neuron_name)
#neuron_name_frm$neuron_name

dim(LMmetricsWhole)
LMmetricsWhole <- ClnMrgName(LMmetricsWhole,neuron_name)
colnames(LMmetricsWhole)
LMmetricsWhole$neuron_name

#check for outliers in the data
(LMmetricsWhole$Taper_1_max)
(LMmetricsWhole$Diameter_min)
#diamoutlierpos <- with(LMmetricsWhole, which(LMmetricsWhole$Diameter_min <= 0.05))
#diameter_min <= 0.05
diamoutlierpos <- with(metricDataWhole, which(metricDataWhole$Diameter_min <= 0.05))
#DiamOutliers <-  LMmetricsWhole[diamoutlierpos,c("neuron_name","Diameter_min","archive_name")]
DiamOutliers <-  metricDataWhole[diamoutlierpos,c("neuron_name","Diameter_min","archive_name")]
str(DiamOutliers)
DiamOutliers <- DiamOutliers[order(DiamOutliers$Diameter_min),]

#diameter_max <= 0.05
diamoutlierpos_max <- with(metricDataWhole, which(metricDataWhole$Diameter_max <= 0.05))
DiamOutliers_max <-  metricDataWhole[diamoutlierpos_max,c("neuron_name","Diameter_max","archive_name")]
DiamOutliers_max <- DiamOutliers_max[order(DiamOutliers_max$Diameter_max),]
DiamOutliers_max
LMmetricsWhole <- LMmetricsWhole[-diamoutlierpos,]

#find row index where Taper_1_max has max value>=100
outlierpos <-  with(metricDataWhole, which(metricDataWhole$Taper_1_max>= 100))
Taper_1_max_outlier <- metricDataWhole[outlierpos,c("neuron_name","Taper_1_max","archive_name")]
Taper_1_max_outlier
#remove that from the LMmetricsWhole dataframe
LMmetricsWhole <- LMmetricsWhole[-outlierpos,]
dim(LMmetricsWhole)

#find row index where Pk_max>=60000
outlierpos <-  with(metricDataWhole, which(metricDataWhole$Pk_max>= 60000))
outlierpos

#find row index where N_bifs>=2000
outlierpos <-  with(metricDataWhole, which(metricDataWhole$N_bifs_total_sum>= 2000))
N_bifs_outliers <- metricDataWhole[outlierpos,c("neuron_name","N_bifs_total_sum","archive_name")]
N_bifs_outliers
testmetric <- LMmetricsWhole[-outlierpos,]
length(outlierpos)

plot(metricDataWhole$N_bifs_total_sum)
plot(metricDataWhole$Length_total_sum)
dim(metricDataWhole)
plot(metricDataWhole$N_bifs_total_sum, metricDataWhole$Length_total_sum, main="N_bifs vs total length",xlab="#bifs", ylab="total length")



outlierpos <-  with(testmetric, which(testmetric$Length_total_sum>= 100000))
testmetric <- testmetric[-outlierpos,]
dim(testmetric)
plot(testmetric$N_bifs_total_sum, testmetric$Length_total_sum, main="N_bifs vs total length",xlab="#bifs", ylab="total length")


hist(metricDataWhole$N_bifs_total_sum, width=0.33, offset=0.00, col="blue", main="Histogram of total length")
x <- logTransformData(metricDataWhole$Length_total_sum)
hist(x, width=0.33, offset=0.00, col="blue", main="Histogram of total length")
colnames(metricDataWhole[,14])
y <- logTransformData(metricDataWhole$N_bifs)
str(y)
hist(y, width=0.33, offset=0.00, col="red", main="Histogram of #bifs")
plot(y, x, main="N_bifs vs total length",xlab="#bifs", ylab="total length")
cor(y,x)

x <- logTransformData(metricDataWhole$Contraction_min)
hist(x, width=0.33, offset=0.00, col="blue", main="Histogram of Contraction min")

y <- logTransformData(metricDataWhole$HillmanThreshold_max)
hist(y, width=0.33, offset=0.00, col="red", main="Histogram of HillmanThreshold max")
plot(y, x, main="N_bifs vs total length",xlab="#contraction_min", ylab="HillmanThreshold_max")

plot(metricDataWhole$HillmanThreshold_max, metricDataWhole$Contraction_min, main="raw values",xlab="#contraction_min", ylab="HillmanThreshold_max")
hist(log(metricDataWhole$Contraction_min),breaks = 30)

hist(x, width=0.33, offset=0.00, col="blue", xlim=c(0,10),main="Histogram of #bifs, total length")
hist(y, width=0.33, offset=0.33, col="yellow", add=TRUE)

par(ask = TRUE)
plotTesting<- function(dataframe){
  numcols <- ncol(dataframe)
  for(i in 27:numcols){
    mycol <- dataframe[,i]
    if(sapply(mycol,is.numeric)){
      colname <- names(dataframe)[i] 
      print(paste("plotting..",colname)) 
      hist(mycol,breaks=30,xlab=colname,main=paste("Histogram of ",colname," before transformation"))
      t <- logTransformData(mycol)    
      hist(t, breaks=30, xlab=paste("log(",colname,")"),main=paste("Histogram of ",colname," after transformation"))    
    }  
  }
}

plotTesting(metricDataWhole)
colnames(metricDataWhole)
dim(metricDataWhole)
#eliminate few metrics that don't look useful
reducedMetricDataWhole <- metricDataWhole[, !(colnames(metricDataWhole) %in% c("Diameter_min","Diameter_max","Diameter_sd","Surface_min","Volume_min","Taper_1_avg","Taper_2_avg","Taper_2_max","Daughter_Ratio_avg","Parent_Daughter_Ratio_min","Parent_Daughter_Ratio_max","Rall_Power_max", "Rall_Power_sd","Pk_min","Bif_ampl_local_min","Bif_ampl_local_avg", "Bif_ampl_local_sd","Bif_ampl_remote_min","Bif_ampl_remote_avg","Bif_ampl_remote_sd","Bif_torque_local_avg","Helix_min","Helix_avg","Fractal_Dim_avg"))]
colnames(reducedMetricDataWhole[,27:116])
#check for most highly correlated pairs
leastcorrelated(reducedMetricDataWhole[,27:116],0.55,10)

cvframe <- data.frame("metrics"=colnames(LMmetricsWhole), "coefficient of variation" = cvarr)
cvframe[order(cvframe$coefficient.of.variation,decreasing=TRUE),]
dim(cvframe)
LMmetricsWhole <- ClnMrgName(LMmetricsWhole,neuron_name)

#tag the groups that are classified with the metrics data

#apply skewness function to all cols in the data frame


#text(LMmetricsWhole$neuron_name,cex=0.7, pos=4, col="red") # add labels
statval <- apply(LMmetricsWhole[,2:115],2,skewness)
summary(statval)
length(statval)

hist(LMmetricsWhole$Taper_1_max,breaks=30)
par(ask = TRUE)
numcols <- ncol(LMmetrics)
for(i in 1:numcols){
  mycol <- LMmetrics[,i]
  hist(mycol,breaks=30, xlab=skewarr[i])
}


#compute covariance and correlation matrix on metricsOfmasterMetaData
v <- cov(LMmetrics, y = NULL, use = "na.or.complete", method = "pearson",exact=FALSE)
r <- cov2cor(v)
dim(r)
r[r>0.85]
rownames(LMmetrics)[r>0.85]
colnames(LMmetrics)[r>0.85]
colnames(LMmetrics)
(v)
r <- cor(LMmetrics)
summary(r)
r[r>0.85]
(r)

#####

LMreduced <- mosthighlycorrelated(LMmetrics, 6441)
#get specific column name from a data frame
colnames(LMmetrics[,10])

dim(LMreduced)
(cor(LMmetrics$Taper_1_avg,LMmetrics$Taper_1_sd))

##############

#compute PCA

library(gclus)
lmcor <- cor(LMmetrics)
summary(lmcor)
previ = 1
for(i in seq(9, 114, 9 )){
  #pdf(paste("", i, ".pdf", sep = ""))
  
  lmm <- LMmetrics[,previ:i]
  previ = i
  lmcor <- abs(cor(lmm))
  lmcorColors <- dmat.color(lmcor)
  lmcorOrder <- order.single(cor(lmm))
  cpairs(lmm, lmcorOrder, panel.colors=lmcorColors,gap=.5,main="Variables Ordered and Colored by Correlation")
  
}
  

#dev.off()
#reduced metrics
lmprcomp <- prcomp(LMreduced, center = TRUE, scale = TRUE)
ls(lmprcomp)
summary(lmprcomp)
#choose eigenvalues >= 1
lmprcomp$sdev ^ 2
# pick upto 25 components as the change in variance reaches is pretty close to 0.
scree(lmprcomp, npcs = 97, main = "LM metrics",xlab="Components")
scree(lmprcomp, npcs = 97, type = "line", main = "LM metrics")

#complete metrics
lmprcomp <- prcomp(LMmetrics, center = TRUE, scale = TRUE)
ls(lmprcomp)
summary(lmprcomp)
#choose eigenvalues >= 1
lmprcomp$sdev ^ 2
# pick upto 25 components as the change in variance reaches is pretty close to 0.
scree(lmprcomp, npcs = 114, main = "LM metrics",xlab="Components")
scree(lmprcomp, npcs = 114, type = "line", main = "LM metrics")


# dot  of PC1
load = lmprcomp$rotation
sortLoad1 = load[order(load[,1]),1]
Main1 = "loadings  for PC1"
xlabs = "Variable loadings"
dotchart(sortLoad1, main = Main1, xlab = xlabs, cex = 0.5, col="red")
Main2 = "loadings  for PC2"
sortLoad2 = load[order(load[,2]),2]
dotchart(sortLoad2, main = Main2, xlab = xlabs, cex = 0.5, col="red")

#(prcomp(LMmetrics, scale = TRUE))
(lmprcomp$x[,1],lmprcomp$x[,2]) # make a scatter b/w PC1 & PC2
(lmprcomp$x[,2],lmprcomp$x[,3]) # make a scatter b/w PC2 & PC3
text(lmprcomp$x[,1],lmprcomp$x[,2], , cex=0.7, pos=4, col="red") # add labels

bi(lmprcomp, cex = c(0.5,0.5))

#apply varimax rotation
lmvar = varimax(lmprcomp$rotation)


LMreduced <- LMmetrics[,colnames(LMmetrics)%in%c("Diameter_avg","Height_total_sum","N_bifs_total_sum","Surface_min","Diameter_max","Volume_total_sum","Length_total_sum","PathDistance_max","EucDistance_avg","Volume_max","Volume_min","Taper_2_avg", "Taper_1_avg","Contraction_avg","Branch_pathlength_max","Parent_Daughter_Ratio_min","Fragmentation_avg","Daughter_Ratio_avg","Partition_asymmetry_sd","Rall_Power_min","Parent_Daughter_Ratio_max", "Pk_avg","Pk_classic_max","Pk_2_max", "Bif_ampl_local_min","Bif_ampl_remote_avg","Bif_tilt_local_sd","Bif_tilt_local_max","Bif_tilt_local_avg","Bif_tilt_remote_sd","Bif_torque_local_sd","Bif_torque_remote_avg", "Bif_torque_local_min", "HillmanThreshold_max", "HillmanThreshold_avg", "Diam_threshold_min", "Helix_min")]
v_r <- cov(LMreduced, y = NULL, use = "na.or.complete", method = "pearson")
r_r <- cov2cor(v_r)
dim(r_r)
(r_r)

lmcor <- cor(LMreduced)
summary(lmcor)
previ = 1
for(i in seq(9, 37, 9 )){
  #pdf(paste("", i, ".pdf", sep = ""))
  
  lmm <- LMreduced[,previ:i]
  previ = i
  lmcor <- abs(cor(lmm))
  lmcorColors <- dmat.color(lmcor)
  lmcorOrder <- order.single(cor(lmm))
  cpairs(lmm, lmcorOrder, panel.colors=lmcorColors,gap=.5,main="Variables Ordered and Colored by Correlation")
  
}

lmprcomp <- prcomp(LMreduced, center = TRUE, scale = TRUE)
ls(lmprcomp)
summary(lmprcomp)
#choose eigenvalues >= 1
lmprcomp$sdev ^ 2
# pick upto 25 components as the change in variance reaches is pretty close to 0.
scree(lmprcomp,npcs = 26, main = "LM metrics",xlab="Components")
scree(lmprcomp,npcs = 26, type = "line", main = "LM metrics")
# dot  of PC1
load = lmprcomp$rotation
sortLoad1 = load[order(load[,1]),1]
Main1 = "loadings  for PC1 after removing correlated variables"
xlabs = "Variable loadings"
dotchart(sortLoad1, main = Main1, xlab = xlabs, cex = 0.5, col="red")
Main2 = "loadings  for PC2 after removing correlated variables"
sortLoad2 = load[order(load[,2]),2]
dotchart(sortLoad2, main = Main2, xlab = xlabs, cex = 0.5, col="red")

#(prcomp(LMmetrics, scale = TRUE))

bi(lmprcomp, cex = c(0.5,0.5))

############Removed level2 part in threeHierarchyGrps##########

addgrps2MetricDF <- function(mydata,specifications=NULL){
  grpcol <- mydata$neuron_name
  tot_len <- 0
  print(length(grpcol))
  if (!is.null(specifications)){
    #print(str(specifications))
    print(length(specifications))
    GrpnameList <- list()
    FilteredSpecList <- list()
    neuronList <- list()
    #tmpGrpList <- list()
    for (i in length(specifications):1){
      #specType <- names(specifications)[i]
      #print(specType)
      if (is.list(specifications)){
        listelement <- specifications[[i]]
        GrpSpec <- listelement[1]$GroupSpec
        Grpname <- paste(GrpSpec,collapse="/")
        #Filtering redundant groups using string comparison on Grpname
        if(length(FilteredSpecList)==0){
          print("add first element to the list ")
          FilteredSpecList <- append(FilteredSpecList,list(Grpname))
          print(paste("lengh of FilteredSpecList:",length(FilteredSpecList)))
          print(paste("testing..",FilteredSpecList[1]))
          neuronList <- unlist(listelement[4])
          #tmpGrpList <- neuronList
          #mydata <- add2GrpCol(mydata,neuronList,length(FilteredSpecList))
          #print(unique(mydata$groupnum))
          tot_len <- length(neuronList)
        }
        #loop through FilteredSpecList and compare against the new Group names in Grpname
        else{
          print(paste("lengh of FilteredSpecList:",length(FilteredSpecList)))
          for(i in 1:length(FilteredSpecList)){
            existingGrp <- FilteredSpecList[[i]]
            #print(existingGrp)
            tmp <- substr(existingGrp,1,nchar(Grpname))
            #print(paste("substring:",tmp,"Grpname:",Grpname))
            if(length(tmp)>0 & Grpname==tmp){
              #print(paste("Skipping redundant groups:",Grpname))
              print(paste("skipping..",Grpname))
            }
            else{
              FilteredSpecList <- append(unique(FilteredSpecList),list(Grpname))
              print(paste("adding..",Grpname))
              neuronList <- unlist(listelement[4])
              #print(unique(neuronList))
              #tmpGrpList <- append(unique(tmpGrpList), unique(neuronList))
              print(paste("before..",length(FilteredSpecList)))
              #mydata <- add2GrpCol(mydata,neuronList,length(FilteredSpecList))
              #print(unique(mydata$groupnum))
              tot_len <- tot_len + length(neuronList)
              break
            }
          }
        }
        print(paste("neuronList:",Grpname, tot_len))
      }
    }
  }
  print(paste("tot_len",tot_len))
  print(FilteredSpecList)
  return (mydata)
}

colnames(plotDf)
dim(plotDf)
p <- subset(plotDf,!is.na(plotDf$groupnum))
dim(p)
diamoutlierpos <- with(p, which(p$Branch_Order_max >= 20 & p$Branch_Order_max <= 35))
DiamOutliers <-  p[diamoutlierpos,c("neuron_name","region2","groupnum","Branch_Order_max","archive_name")]
DiamOutliers
diamoutlierpos <- with(p, which(p$Bif_ampl_remote_avg >= 40 & p$Bif_ampl_remote_avg <= 50))
DiamOutliers <-  p[diamoutlierpos,c("neuron_name","region2","groupnum","Branch_Order_max","archive_name")]
DiamOutliers
unique(p[!is.na(p$groupnum),"groupnum"])
m1 <- "Branch_Order_max"
m2 <- "Bif_ampl_remote_avg"
x <- logTransformData(p[,m1])
y <- logTransformData(p[,m2])
x <- p[,m1]
length(x)
mean(x)/length(x)
y <- p[,m2]
#c(unique(plotDf[!is.na(plotDf$groupnum),"groupnum"]))
plot(x, y,col = p$groupnum, xlim=c(5,35), ylim = c(5,80), main=paste(m1," vs ",m2),xlab=m1, ylab=m2)
legend("bottomright",legend = c("CA1","CA3"), col = c(2,1),pch = 1, cex=.7)
#text(x, y, plotDf$groupnum, cex=0.7, pos=4, col="red")


library(ggplot2)
#plots crosshair of mean and SD on xy axes. Also computes if pair of groups are separated without overlapping.
crosshairSeparation<- function(s,group){
  ctr <- 0
  for(i in 1:(length(s)-1)){
    for(j in (i+1):length(s)){
      # x vs y
      m1 <- s[,i]
      m2 <- s[,j]
      print(paste("comparing",m1,m2))
      d <- sqrt((m1$x-m2$x) ^ 2 + (m1$y-m2$y) ^ 2)
      if(m1$xSD/d <= 1 && m2$ySD/d <= 1){
        ctr <- ctr + 1
      }
    }
  }
  print(ctr)
  #s[i,"SepGrps"] <- ctr
  return (ctr)
}

print(paste(m1,"vs",m2))
#xlabel = column1,
#ylabel = column2,
#first compute the summary metrics N, mean, and SD on X and y axes
z <- ddply(p,.(groupnum,group),function(dfp,column1,column2) {
  c(N = nrow(dfp),
    x = mean(dfp[,column1]),
    y = mean(dfp[,column2]),
    xSE = sqrt(var(dfp[,column1]))/nrow(dfp),
    ySE = sqrt(var(dfp[,column2]))/nrow(dfp),
    xSD = sd(dfp[,column1]),
    ySD = sd(dfp[,column2]))}, m1, m2)
print(z)
#compare how well the means b/w two groups are separated 
#print(paste("comparing ",z$group[1],"and",z$group[2]))
for(ii in 1:(nrow(z)-1)){
  #print(ii)
  row1 <- z[ii,]
  for(jj in (ii+1):nrow(z)){
    #print(jj)
    row2 <- z[jj,]
    #print(paste("comparing ",row1$group,"and",row2$group))
    #compute the difference of means on x-axis
    xmeanDistance = row1$x-row2$x
    #the distabce of means should be greater than the SDs, this will allow overlaps 
    xmaxdis <- xmeanDistance > row1$xSD && xmeanDistance > row2$xSD
    #compute the difference of means on y-axis
    ymeanDistance = row1$y-row2$y
    ymaxdis <- ymeanDistance > row1$ySD && ymeanDistance > row2$ySD
    if(xmaxdis && ymaxdis){
      print(paste(row1$group, "and", row2$group, "is a non-overlapping pair"))  
    }
  }
}
#geom_text(aes(label = paste("ln(",m1,")"))) +
#geom_text(aes(label = paste("ln(",m2,")"))) +
#xlab(paste("ln(",m1,")")) +
#ylab(paste("ln(",m2,")")) +
#,xlab=xlabel,ylab=ylabel
#plot the crosshairs
ggplot(data=z,aes(x=x,y=y,size=0.1,label=FALSE)) +
  geom_point(aes(colour = group)) +
  xlab(m1) +
  ylab(m2) +
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD, colour = group)) + 
  geom_errorbarh(aes(xmin = x - xSD, xmax = x + xSD, colour = group)) +
  geom_errorbar(aes(ymin = y - ySE, ymax = y + ySE, colour = group)) + 
  geom_errorbarh(aes(xmin = x - xSE, xmax = x + xSE, colour = group)) +
  scale_colour_hue(name="cell types", # Legend label, use darker colors
                   breaks=unique(group),
                   labels=unique(group),
                   l=50,       # Use darker colors, lightness=50
                   c=100)      #chroma (intensity of color)
}

n <- length(unique(group))
print(n*(n-1)/2-n)
print((n^2)/2-n)
bestpairs <- list()
#loop through metrics for all cominations n(n-1)/2 and pick the best uncorrelated pair
for(i in 33:50){
  for(j in (i+1):51){
    npairs = crosshairscatter(p,names(p)[i],names(p)[j],group)
    #bestpairs <- append(bestpairs,list(npairs))
    if(npairs >= n){
      print(paste(names(p)[i],names(p)[j],npairs))
      bestpairs <- append(bestpairs,paste(names(p)[i],names(p)[j],npairs))
    }
  }
}
length(bestpairs)
bestpairs

#compare how well the means b/w two groups are separated 
#print(paste("comparing ",z$group[1],"and",z$group[2]))
for(ii in 1:(nrow(z)-1)){
  #print(ii)
  row1 <- z[ii,]
  for(jj in (ii+1):nrow(z)){
    #print(jj)
    row2 <- z[jj,]
    #print(paste("comparing ",row1$group,"and",row2$group))
    #compute the difference of means on x-axis
    xmeanDistance = row1$x-row2$x
    #the distabce of means should be greater than the SDs, for good separation 
    xmaxdis <- xmeanDistance >= row1$xSD && xmeanDistance >= row2$xSD
    #compute the difference of means on y-axis
    ymeanDistance = row1$y-row2$y
    ymaxdis <- ymeanDistance >= row1$ySD && ymeanDistance >= row2$ySD
    if(xmaxdis && ymaxdis){
      print(paste(row1$group, "and", row2$group, "is a non-overlapping pair"))  
      ctr = ctr+ 1
    }
  }
}

#overlay the dotted crosshairs for 'Others' group
ggplotobj <- chp + 
  #geom_point(aes(labels = unique(group))) +
  geom_errorbar(data=z, aes(x=x, y=y, ymin = y - ySD, ymax = y + ySD, colour = group),linetype="dashed") + 
  geom_errorbarh(data=z, aes(x=x, y=y, xmin = x - xSD, xmax = x + xSD, colour = group), linetype="dashed") +
  scale_colour_manual(values = unique(group), labels = unique(group))
#geom_text(aes(label="Others"))

#testing this part...
if(length(typSpec)>1){
  #get unique values from colnmvar 
  nwparamlst <- unique(usedMetaData[,colnmvar])
  for(i in 1:length(nwparamlst)){
    typSpec <- append(typSpec, list(nwparamlst[i]))
    #Assign name to newly added value to the specification list
    names(typSpec)[length(typSpec)] = colnmvar
    
  }
  
  printGrps(ctxLobesOnlyHierarchy)
  
  primary <- list(expercond="Control", region1="Hippocampus")
  level1 <- "region2"
  level2 <- "region3"
  HipOnlyHierarchy <- threeLevelHierarchyGrping(primary, metricwhole, typeMetaData= usedMetaData, level1, level2,minSize=70)
  printGrps(HipOnlyHierarchy)
  #HipOnlyHierarchy <- HipOnlyHierarchy[-6]#six times
  HipOnlyHierarchy <- HipOnlyHierarchy[-6]#five times
  
  primary <- list(expercond="Control", region1="Olfactory bulb")
  level1 <- "region2"
  level2 <- "region3"
  OlfOnlyHierarchy <- threeLevelHierarchyGrping(primary, metricwhole, typeMetaData= usedMetaData, level1, level2,minSize=70)
  printGrps(OlfOnlyHierarchy)
  #OlfOnlyHierarchy <- OlfOnlyHierarchy[-5]#four times
  OlfOnlyHierarchy <- OlfOnlyHierarchy[-5]#two times
  
  #brMetricHierarchy <- ctxLobesOnlyHierarchy
  #brMetricHierarchy <- append(brMetricHierarchy,HipOnlyHierarchy)
  #brMetricHierarchy <- append(brMetricHierarchy,OlfOnlyHierarchy)
  #printGrps(brMetricHierarchy)
  
  brMetricHierarchy_1 <- makeHierarchyGroups(primary, metricwhole, nctxData, level1, level2,level3,minSize=300)
  #x <- getSubset2(metricwhole, primary, usedMetaData)
  #length(x)
  #ctxallHierarchy <- makeHierarchyGroups(primary, metricwhole, usedMetaData, level1, level2,level3,minSize=300)
  #printGrps(ctxallHierarchy)
  #regions not Neocortex
  nctxData <- subset(usedMetaData,region1=="Neocortex")
  dim(nctxData)
  nonctxData <- subset(usedMetaData,region1!="Neocortex")
  dim(nonctxData)
  brMetricHierarchy_1 <- makeHierarchyGroups(primary, metricwhole, nctxData, level1, level2,level3,minSize=300)
  printGrps(brMetricHierarchy_1)
  
  brMetricHierarchy_2 <- makeHierarchyGroups(primary, metricwhole, nonctxData, level1, level2,level3,minSize=300)
  printGrps(brMetricHierarchy_2)
  
  brMetricHierarchy <- brMetricHierarchy_1
  brMetricHierarchy <- append(brMetricHierarchy,brMetricHierarchy_2)
  
  str(pcawhole_matrix)
  pcasumm <- summary(pcawhole_matrix)
  str(pcasumm)
  names(pcasumm)
  
  dim(pcawhole_matrix$x)
  str(pcawhole_matrix)
  plot(pcawhole_matrix$x)
  plot(pcawhole_matrix$x[,1],pcawhole_matrix$x[,2],col=c(1,2)) # make a scatterplot
  
  
  pcawhole_matrix$x
  
  screeplot(pcawhole_matrix, npcs = 32, type = "lines", main = "whole arbor metrics")#,xlab="Components")
  
  screeplot(pcawhole_matrix,npcs = 7, type = "lines", main = "whole arbor metrics")
  loadings <- pcawhole_matrix$rotation
  dim(loadings)
  #square loadings for easy comparison
  loadings2 <- loadings[,1:2]^2
  
  #the sum of square of loadings for each PC is equal to 1
  sum(loadings2[,1])
  #the most contributing variables for PC1
  sortLoad1 <- loadings2[order(loadings2[,1]),1]
  Main1 = "loadings for PC1"
  xlabs = "Variable loadings^2"
  dotchart(sortLoad1, main = Main1, xlab = xlabs, cex = 0.5, col="red")
  #the most contributing vairables for PC2
  sortLoad2 <- loadings2[order(loadings2[,2]),2]
  Main1 = "loadings for PC2"
  xlabs = "Variable loadings^2"
  dotchart(sortLoad2, main = Main1, xlab = xlabs, cex = 0.5, col="red")
  # Determine number of clusters
  pcawhole[,2:33]
  dim(pcawhole)
  str(pcawhole_matrix)
  dim(pcawhole_matrix$x)
  pcawhole_matrix$x[,1]
  kmeansObj <- kmeans(pcawhole_matrix$x[,1:32],centers=4)
  
  names(kmeansObj)
  kmeansObj$withinss
  plot(pcawhole[,3],pcawhole[,9],col=c(1,2))
  plot
  plot(pcawhole_matrix$x[,3],pcawhole_matrix$x[,4],col=c(1,2))
  #wss <- (nrow(pcawhole_matrix$x[,1:32])-1)*sum(apply(pcawhole[,2:33],2,var))
  #length(wss)
  
  
  unique(kmeansObj$cluster)
  
  pairs(~PC1+PC2+PC3+PC4,data=pcawhole_matrix$x,col=kmeansObj$cluster,main="k-means scatterplot Matrix")
  
  
  plot(pcawhole_matrix$x[,1:2],col=kmeansObj$cluster)
  plot(pcawhole_matrix$x[,2:3],col=kmeansObj$cluster)
  plot(pcawhole_matrix$x[,3:4],col=kmeansObj$cluster)
  plot(pcawhole_matrix$x[,4:5],col=kmeansObj$cluster)
  t <- subset(pcawhole_matrix,kmeansObj$cluster==1)
  dim(t)