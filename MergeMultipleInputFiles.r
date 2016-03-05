#open metadata file of v5.4 
# Load meta data
# Set metadata filename
#metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/all_7986_v5.4_metadata_1.csv"
#masterMetaData <- read.delim(metadataFilename,header = TRUE, sep="}",na.strings=c("NA",""),quote="\"",stringsAsFactors=FALSE)
#print("finished reading...")
#removing the extra column named 'X'
#masterMetaData <- masterMetaData[, !(colnames(masterMetaData) %in% c("X"))]

#metricAxonFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/Report_axon.txt"
#metricDataAxon <- read.delim(metricAxonFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#colnames(metricDataAxon)[1] <- 'neuron_name'
#print("finished reading Axon...")

#metricDendFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/Report_dendrite.txt"
#metricDataDend <- read.delim(metricDendFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#colnames(metricDataDend)[1] <- 'neuron_name'
#print("finished reading Dend...")

#metricApicalFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/Report_apical.txt"
#metricDataApical <- read.delim(metricApicalFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#colnames(metricDataApical)[1] <- 'Neuron_Name'
#print("finished reading Apical...")

#metricApiBasFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/Report_apicalbasal.txt"
#metricDataApiBas <- read.delim(metricApiBasFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#colnames(metricDataApiBas)[1] <- 'Neuron_Name'
#print("finished reading apical & basal...")

#check for -ve values in metric data frame
chkNegVals <- function(metricDataframe){
  x <- metricDataframe
  negcols <- apply(metricDataframe[,sapply(x,is.numeric)], 2, function(col) any(col < 0)) 
  #print(length(negcols))
  print(negcols[negcols==TRUE])
  return(negcols[negcols==TRUE])
}

rmvNegVals <- function(metricDataframe, fnlist){
  flist <- unlist(fnlist)
  print(str(flist))
  #if(is.list(flist)){
    for( i in 1:length(flist)){
      tmp <- flist[i]
      print(tmp)
      negpos <- with(metricDataframe, which(metricDataframe$tmp < 0))
      print(negpos)
      cols <-c("neuron_name", "archive_name")
      cols <- append(cols, list(as.character(tmp)))
      if(length(negpos)>0){
        negrows <- metricDataframe[negpos,cols]
        print(negrows)
        metricDataframe <- metricDataframe[-negpos]
      }else
        print("no negative data found!")
    }
  #}
  return (metricDataframe)
}

#check outlier values that are greater than a given val in the metrics
chkOutlierVals <- function(metricDataframe,lowlim, uplim){
  outcols <- apply(metricDataframe, 2, function(col) any(col > lowlim & col<uplim)) 
  print(outcols[outcols>val])
}

#count NAs per each row in the data frame
cntNRmvNAs <- function(metricDataframe,NAtotal){
  nasperrow <- apply(is.na(metricDataframe),1,sum)
  #print(nasperrow)
  #print values that have NAs same as NAtotal
  cntNArows <- nasperrow[nasperrow==NAtotal]
  #count the no.of rows that are NAs
  print(paste("Found ", length(cntNArows), NAtotal, "NA rows from the table"))
  #check the rownames and colnames of the NA values
  #print(rownames(metricDataframe)[nasperrow>0])
  #print(nasperrow[nasperrow>0])
  #print(paste("no more columns with NAS", colnames(metricDataframe)[naspercol>0]))
  #list neuron_name that are NAs 
  #print(metricDataframe[nasperrow>0,"neuron_name"])
  #remove the NA row from metricDataframe
  retDataframe <- metricDataframe[nasperrow!=NAtotal,]
  
  #print(dim(retDataframe))
  return (retDataframe)

 
}

#Remove the extension
rmvExtnsn <- function(mydataframe){
  mydataframe$neuron_name <- sapply(strsplit(mydataframe$neuron_name, "\\."), "[[", 1)
  mydataframe <- mydataframe[, !(colnames(mydataframe) %in% c("X"))]
  print(dim(mydataframe))
  return (mydataframe)
}

#function to remove .CNG.swc from Neuron_Name column, and remove the extra "X" column and merge metadata with metrics
clnData <- function(mydataframe){
  mydataframe$neuron_name <- sapply(strsplit(mydataframe$neuron_name, "\\."), "[[", 1)
  mydataframe <- mydataframe[, !(colnames(mydataframe) %in% c("X"))]
  print(dim(mydataframe))
  #redDataframe <- cntNRmvNAs(mydataframe,113)
  #print("Merging metadata with metrics..")
  return (mydataframe)
}

#function to check for NA values and merge metrics with neuron_name column
ClnMrgName <- function(mydataframe){
  mydataframe <- mydataframe[, !(colnames(mydataframe) %in% c("X"))]
  print(dim(mydataframe))
  redDataframe <- cntNRmvNAs(mydataframe)
  print(paste("Merging metadata with metrics..",colnames(mydataframe$neuron_name)))
  mrgddataframe <- cbind(neuron_name,redDataframe)
  print(dim(mrgddataframe))
  return (mrgddataframe)
}

#read only the first column after the split at '.'

dim(metricDataWhole)
dim(masterMetaData)
#metricDataWhole$Neuron_Name <- sapply(strsplit(metricDataWhole$Neuron_Name, "\\."), "[[", 1)
#remove extenstions in neuron_name column and merge metadata with metrics

#metricDataAxon <- clnData(metricDataAxon)
#dim(metricDataAxon)
#chkNegVals(metricDataAxon)
#names(metricDataAxon)

#dim(metricDataDend)

#metricDataDend <- clnData(metricDataDend)
#chkNegVals(metricDataDend)
#names(metricDataDend)
#This group contains cells that have any of apical or basal dendrites
#dim(metricDataApical)
#metricDataApical <- clnData(metricDataApical)
#chkNegVals(metricDataApical)
#names(metricDataApical)

#dim(metricDataApiBas)
#metricDataApiBas <- clnData(metricDataApiBas)
#chkNegVals(metricDataApiBas)
#names(metricDataApiBas)

#metricDataWhole$Neuron_Name <- sapply(strsplit(metricDataWhole$Neuron_Name, "\\."), "[[", 1)
#metricDataWhole <- metricDataWhole[, !(colnames(metricDataWhole) %in% c("X"))]

#metricDataAxon$Neuron_Name <- sapply(strsplit(metricDataAxon$Neuron_Name, "\\."), "[[", 1)
#metricDataAxon <- metricDataAxon[, !(colnames(metricDataAxon) %in% c("X"))]

#metricDataDend$Neuron_Name <- sapply(strsplit(metricDataDend$Neuron_Name, "\\."), "[[", 1)
#metricDataDend <- metricDataDend[,!(colnames(metricDataDend) %in% c("X"))]


#checking for mismatch is neuron names between masterMetaData and metricDataWhole
#testoverlap <- intersect(masterMetaData$neuron_name, metricDataWhole$Neuron_Name)
#testdiff1 <- setdiff(masterMetaData$neuron_name, metricDataWhole$Neuron_Name)
#testdiff2 <- setdiff(masterMetaData$neuron_name, metricDataWhole$Neuron_Name)

#masterMetaData <- merge(masterMetaData,metricDataWhole,by.x="neuron_name",by.y="Neuron_Name")
#names(masterMetaData)
#dim(masterMetaData)

#make a copy of just the metrics cols
#neuron_name <- masterMetaData[,1]
#metricsOfmasterMetaData <- masterMetaData[,27:146]
#metricsOfmasterMetaData <- cbind(neuron_name,metricsOfmasterMetaData)

#make a copy of just metadata
#metadataOfMaster <- masterMetaData[,1:26]


#nasperrow <- apply(is.na(metricsOfmasterMetaData),1,sum)
#print the values that are >0
#nasperrow[nasperrow>0]
#count NAs per each column
#naspercol <- apply(is.na(metricsOfmasterMetaData),2,sum)
#naspercol[naspercol>0]
#check the rownames and colnames of the NA values
#rownames(metricsOfmasterMetaData)[nasperrow>0]
#colnames(metricsOfmasterMetaData)[naspercol>0]

#colnames(metricsOfmasterMetaData)
#list neuron_name that are NAs 
#metricsOfmasterMetaData[nasperrow>0,"neuron_name"]
#remove the NA row from metrics part
#metricsOfmasterMetaData <- metricsOfmasterMetaData[nasperrow==0,]
#remove the NA row from metadata part
#metadataOfMaster <- metadataOfMaster[nasperrow==0,]

#combine the metrics and metadata after removing NA values.all.y=F will drop rows in y that are not matched with x.
#masterMetaData <- merge(metadataOfMaster,metricsOfmasterMetaData,by.x="neuron_name",by.y="neuron_name",all.x=T,all.y=F)

#print coefficent of variation for each metric 
printCV <- function(metricdataframe){
  metricmean <- sapply(metricdataframe,mean)
  #zeromean <- metricmean[metricmean > -1 & metricmean < 1]
  #print(paste("zero mean", zeromean))
  #metricmeanfrmt <- format(metricmean,digits =5,nsmall=3,justify="left",trim = FALSE)
  metricsd <- sapply(metricdataframe,sd)
  #metricsdfrmt <- format(metricsd,digits =5,nsmall=3,justify="left",trim = FALSE)
  #cv <- sapply(metricdataframe,sd/mean)
  print(paste("mean of", length(metricmean), "metrics:"))
  print(metricmean)
  print(paste("sd of", length(metricsd), "metrics:"))
  print(metricsd)
  cv <- (metricsd/metricmean)
  print(paste("cv of", length(cv), "metrics:"))
  print(cv)
}

#method 2 for calculating CV only for the numeric columns
cvfun <- function(x) (apply(x[,sapply(x,is.numeric)],2,sd)/colMeans(x[,sapply(x,is.numeric)]) )

#function to remove columns with mean <= 0
chkMeanCols <- function(metricdataframe){
  metricmean <- sapply(metricdataframe,mean)
  zeromeanpos <-  with(metricdataframe,which(metricmean<0))
  print(paste(zeromeanpos,colnames(metricdataframe[,zeromeanpos])))
  #redframe <- metricdataframe[,-zeromeanpos]
  #return (redframe)
}


#compute p-values for least correlated variables
# apply cor.test() with extract to each element of the list
mostcorrelated <- function(mydataframe,threshold=0.8,numreport=50)
{
  ctr = 0   
  # find the correlations
  cormatrix <- cor( mydataframe, y = NULL, use = "na.or.complete", method = "pearson")
  print(dim(cormatrix))
  #print(summary(cormatrix))
  # set the correlations on the diagonal or lower triangle to zero,
  # so they will not be reported as the highest ones:
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  # find the dimensions of the matrix, and the row names:
  numrows <- nrow(cormatrix)
  therownames <- rownames(cormatrix)
  #print(ls(cormatrix))
  # find the highest correlations
  sorted <- sort(abs(cormatrix),decreasing=TRUE)
  #nonzerocor <- sorted[1:length(cormatrix)/2+1]
  nonzerocor <- sorted[sorted>0]
  print(length(nonzerocor))
  #print(summary(nonzerocor))
  print(paste(ctr,"<=",numreport))
  #dataframe to return with features and corresponding R^2 values sorted by least correlated
  featurecorr <- data.frame("f1"='f1',"f2"='f2',"r2"=0.0,"pval"=0.0)
  for (i in length(nonzerocor):1)
  {
    if(ctr<=numreport){
      corri <- nonzerocor[i]
      #print(corri)
      # find the pair of variables with this correlation
      for (j in 1:(numrows-1))
      {
        for (k in (j+1):numrows)
        {
          corrjk <- cormatrix[j,k]
          #cor <- cor(cormatrix[,j],cormatrix[,k])
          #print(paste("cor= ",corri," corrjk= ",corrjk))
          if(corri == abs(corrjk) && corri^2>threshold){
            #ctr = ctr+1
            rowname <- therownames[j]
            colname <- therownames[k]
            #print(paste("**",length(mydataframe[,rowname]),length(mydataframe[,colname])))
            pval <- cor.test(x = mydataframe[,rowname],y = mydataframe[,colname],method="pearson")$p.value
            #add the least correlated feature pairs
            featurecorr <- rbind(featurecorr,data.frame("f1"=rowname,"f2"=colname,"r2"= round(corri^2,digits=5),"pval"=round(pval,digits=5)))
            print(paste("i= ",i,"variables ",rowname, "and ", colname,"correlation= ",corrjk,"r^2= ",round(corri^2,digits=5), "p.value=",round(pval,digits=5)))
            
          }
        }
      }
    }else
      break
  } 
  print(paste("#feature pairs compared=",nrow(featurecorr)))
  #remove the first dummy row from featurecorr
  featurecorr <- featurecorr[-1,]
  
  return (featurecorr)
  
}

#removes one of the correlated pairs that has 1 or more correlations
eliminateCorrFeat <- function(corrDf){
  feliminate <- data.frame(feature=character(0),occurence=integer(0),stringsAsFactors=FALSE)
  f1 <- as.character(corrDf[,'f1'])
  f2 <- as.character(corrDf[,'f2'])
  for(i in 1:nrow(corrDf)){
    #read the pair
    a <- as.character(corrDf[i,1])
    b <- as.character(corrDf[i,2])
    #count occurence of a
    repA <- length(f2[f2==a]) + length(f1[f1==a])
    #count occurence of b
    repB <- length(f2[f2==b]) + length(f1[f1==b])
    #add the most occured feature for eliminating from the matrix
    if(repA >=repB && length(agrep(a,feliminate$feature))==0){
      feliminate <- rbind(feliminate, data.frame(feature=a, occurence=repA))
    }else if(repB >= repA && length(agrep(b,feliminate$feature))==0){
      feliminate <- rbind(feliminate, data.frame(feature=b, occurence=repB)) 
    }   
  }
  return (feliminate)
}

adjustpval <- function(col,adjustby){
  col <- col*adjustby
  names(col) <- "pval_bonferroni"
  return(col)
}


#find skewness,mean,median of each column passed as array parameter
skewness <- function(x){
  m3 <- mean((x-mean(x))^3)
  skew <- m3/(sd(x)^3)
  #mean <- mean(x)
  #median <- median(x)
  #stat3 <- c(skew,mean,median)
  return(skew)
}

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#natural logarithm
ln <- function(DataCol){
  return (log(DataCol))
}

#log transform the skewed distribution
logTransformData <- function(myData,threshold=1,maxIterations=2){
  dataSkewness <- skewness(myData)
  print(paste("skewness:",dataSkewness))
  #only distributions with skewness > 1 will be transformed
  if (abs(dataSkewness) > threshold){
    if (dataSkewness < 0){
      #making sure the raw data > 1 for negatively skewed observations
      myData <- -myData + 1.5 + max(myData)
     
    }
    iterations <- 0
    repeat {
      if (min(myData) <= 1){
        #making sure the raw data >=1 for positiviely skewed observations
        print(min(myData))
        myData <- myData + (1-min(myData))
      }
      myData <- ln(myData)
      iterations <- iterations + 1
      dataSkewness <- skewness(myData)
      if (dataSkewness <= threshold || (maxIterations != -1 && iterations >= maxIterations)){
        break;
      }
    } 
  }
  return(myData)
}

  

