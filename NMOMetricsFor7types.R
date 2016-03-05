
source("C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/LoadMetricFiles.r")
dim(masterMetaData)

#metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/all_8858_v5.6_metadata.csv"
#masterMetaData <- read.csv(metadataFilename,header = TRUE, sep="}",na.strings=c("NA",""),quote="\"",stringsAsFactors=FALSE)
#print("finished reading...")
unique(masterMetaData$species)
unique(masterMetaData$region1)
unique(masterMetaData$region2)
unique(masterMetaData$region3)
colnames(masterMetaData)
masterMetaData$neuron_name

#removing duplicate rows
dpos <- duplicated(masterMetaData)
length(dpos[dpos==TRUE])
dim(masterMetaData)
#combining soma_surface, N_stems with other metrics
metricWholeFile_1 <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_wholearbor_10022.txt"
#read the part with soma_surface
metricWholeFile_2 <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_somaonly_10022.txt"

neuronByDateFile <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/neuronnamebydate.csv"

neuronByDate <- read.delim(neuronByDateFile, header = TRUE, sep="\t", stringsAsFactors=FALSE)
dim(neuronByDate)
colnames(neuronByDate)

Whole_1 <- read.delim(metricWholeFile_1, header = TRUE, sep="\t",dec=".",stringsAsFactors=FALSE)

#rename the column name 
colnames(Whole_1)[1] <- 'neuron_name'
colnames(Whole_1)
dim(Whole_1)
Whole_1 <- rmvExtnsn(Whole_1)
Whole_1$neuron_name
colnames(Whole_1)
#making sure no extraneous columns like "x" and "Soma_Surface_total_sum" are not in the first part
Whole_1 <- Whole_1[, !(colnames(Whole_1) %in% c("X"))]

Whole_2 <- read.delim(metricWholeFile_2, header = TRUE, sep="\t",dec=".",stringsAsFactors=FALSE)
#Whole_2 <- read.table(metricWholeFile_2, header = TRUE, sep="\t", dec=".")
colnames(Whole_2)[1] <- 'neuron_name'
colnames(Whole_2)
Whole_2 <- rmvExtnsn(Whole_2)
#largeSomaVals <- Whole_2[with(Whole_2, which(Whole_2$Soma_Surface_total_sum>10000000)),c("neuron_name")]
#Whole_1[with(Whole_1,which(Whole_1$neuron_name%in%largeSomaVals)),]
#masterMetaData[with(masterMetaData,which(masterMetaData$neuron_name%in%largeSomaVals)),c("archive_name","neuron_name")]
#Whole_2$Soma_Surface_total_sum[with(Whole_2,which(Whole_2$neuron_name=='0-L5int-na'))]
#remove extraneous columns
Whole_2 <- Whole_2[, !(colnames(Whole_2) %in% c("X"))]
colnames(Whole_2)
dim(Whole_1)
dim(Whole_2)
dim(masterMetaData)
Whole_2[with(Whole_2, which(Whole_2$neuron_name=='0-2b')),]
#Whole_2[with(Whole_2, which(Whole_2$neuron_name=='TM-pyr-mouse-Oct-7-2009-B-b-R-ASC')),]

str(Whole_2)
#merge all metrics over neuron_name
t <- setdiff(Whole_1$neuron_name,Whole_2$neuron_name)
#difference between database and dableFiles
t <- setdiff(masterMetaData$neuron_name,Whole_1$neuron_name)
t <- setdiff(Whole_1$neuron_name, masterMetaData$neuron_name)
masterMetaData[with(masterMetaData,which(masterMetaData$neuron_name%in%t)),c("archive_name","neuron_name")]
Whole_2[with(Whole_2,which(Whole_2$neuron_name%in%t)),"neuron_name"]
length(t)
#tt <- aggregate(Whole_1$neuron_name,by=list(Whole_1$neuron_name),length)
#tt[order(tt$x,decreasing=T),]
#colnames(tt)
#head(tt,n=15)
#nrow(tt)
nrow(Whole_1)
#removing duplicate rows
#dpos <- duplicated(Whole_1)
#unique(dpos)
#dups_vals <- subset(Whole_1, dpos==TRUE)
#dups_vals$neuron_name
dpos <- duplicated(Whole_1)
length(dpos[dpos==TRUE])
dup_neuron_names <- Whole_1[dpos,"neuron_name"]
dup_neuron_names
Whole_1 <- subset(Whole_1, dpos==FALSE)
dim(Whole_1)
dpos <- duplicated(Whole_2)
length(dpos[dpos==TRUE])
dup_neuron_names <- Whole_2[dpos,"neuron_name"]
dup_neuron_names
Whole_2 <- subset(Whole_2, dpos==FALSE)
dim(Whole_2)
#colnames(Whole_2)
setdiff(Whole_1$neuron_name,masterMetaData$neuron_name)
setdiff(masterMetaData$neuron_name,Whole_1$neuron_name)
Whole_1$neuron_name[10009]
Whole_2$neuron_name[10009]

#NMOmetricData <- cbind(Whole_2,Whole_1[2:114])
NMOmetricData <- merge(x=Whole_2,y=Whole_1, by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(NMOmetricData)
NMOmetricData$neuron_name[10009]
NMOmetricData <- cntNRmvNAs(NMOmetricData,115)
print("finished reading whole...")
dim(NMOmetricData)
colnames(NMOmetricData)
#attach neuron_ID column to the report
colnames(masterMetaData[,1:2])

dim(masterMetaData)
t <- intersect(NMOmetricData$neuron_name,masterMetaData$neuron_name)
t <- setdiff(NMOmetricData$neuron_name,masterMetaData$neuron_name)
t
length(t)

addNeuronID <- function(idDf, metricDf){
  metricwneuronID <- merge(x=idDf[,c('neuron_name','neuron_id')],y=metricDf,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
  print(dim(metricwneuronID))
  print(colnames(metricwneuronID))
  #metricwneuronID[,c('neuron_name','neuron_id')]  
  return (metricwneuronID)
}

masterMetaforneuronID <- merge(x=masterMetaData[,c('neuron_name','neuron_id')],y=NMOmetricData,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(masterMetaforneuronID)
colnames(masterMetaforneuronID)
masterMetaforneuronID[,c('neuron_name','neuron_id')]


#write the combined wholearbor
write.table(masterMetaforneuronID, file="wholeArborCombined_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)
dim(neuronByDate)
dim(masterMetaforneuronID)
setdiff(neuronByDate$neuron_name,masterMetaforneuronID$neuron_name)
imagecsv <- merge(x=neuronByDate,y=masterMetaforneuronID,by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(imagecsv)
colnames(imagecsv)
write.table(imagecsv, file="tilingimagelookup.csv",quote=FALSE, sep="\t",row.names=FALSE)

#write the basal only
basalFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_dendriteonly_10022.txt"

basal_1 <- read.delim(basalFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(basal_1)[1] <- 'neuron_name'

colnames(basal_1)
dim(basal_1)
#remove extension and remove 'X' column 
basal_1 <- rmvExtnsn(basal_1)
basal_1$neuron_name
colnames(basal_1)
basal_1 <- subset(basal_1, dpos==FALSE)
basalreport <- addNeuronID(masterMetaData[,1:2], basal_1)
dim(basalreport)
colnames(basalreport)
#write the basal Only
write.table(basalreport, file="basalOnly_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)

#write the axon only
axonFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_axononly_10022.txt"

axon_1 <- read.delim(axonFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(axon_1)[1] <- 'neuron_name'
colnames(axon_1)
dim(axon_1)
axon_1 <- rmvExtnsn(axon_1)
axon_1 <- subset(axon_1, dpos==FALSE)
axon_1$neuron_name
colnames(axon_1)
axonreport <- addNeuronID(masterMetaData[,1:2], axon_1)
dim(axonreport)
#write the axon Only 
write.table(axonreport, file="axonOnly_v5.6.txt",quote=FALSE,sep="\t",row.names=FALSE)

#write the apical only
apicalFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_apicalOnly_10022.txt"
apical_1 <- read.delim(apicalFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(apical_1)[1] <- 'neuron_name'
colnames(apical_1)
dim(apical_1)
apical_1 <- rmvExtnsn(apical_1)
apical_1 <- subset(apical_1, dpos==FALSE)
apical_1$neuron_name
colnames(apical_1)
apicalreport <- addNeuronID(masterMetaData[,1:2], apical_1)
dim(apicalreport)
#write the apical Only 
write.table(apicalreport, file="apicalOnly_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)

#write the apical & basal
apicalbasalFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_apicalBasal_10022.txt"

apicalbasal_1 <- read.delim(apicalbasalFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(apicalbasal_1)[1] <- 'neuron_name'
setdiff(Whole_1$neuron_name,Whole_2$neuron_name)
colnames(apicalbasal_1)
dim(apicalbasal_1)
apicalbasal_1 <- rmvExtnsn(apicalbasal_1)
apicalbasal_1 <- subset(apicalbasal_1, dpos==FALSE)
apicalbasal_1$neuron_name
colnames(apicalbasal_1)
apicalbasalreport <- addNeuronID(masterMetaData[,1:2], apicalbasal_1)
dim(apicalbasalreport)
#write the apical basal 
write.table(apicalbasalreport, file="apicalbasal_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)

#write apical and axon 
apicalaxonFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_apicalAxon_10022.txt"

apicalaxon_1 <- read.delim(apicalaxonFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(apicalaxon_1)[1] <- 'neuron_name'
colnames(apicalaxon_1)
dim(apicalaxon_1)
apicalaxon_1 <- rmvExtnsn(apicalaxon_1)
apicalaxon_1 <- subset(apicalaxon_1, dpos==FALSE)
apicalaxon_1$neuron_name
colnames(apicalaxon_1)
apicalaxonreport <- addNeuronID(masterMetaData[,1:2], apicalaxon_1)
dim(apicalaxonreport)
#write the axon Only 
write.table(apicalaxonreport, file="apicalaxon_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)

#write the basal and axon
basalaxonFile<- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/ReportsForNMO_10022/Report_basalAxon_1022.txt"

basalaxon_1 <- read.delim(basalaxonFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(basalaxon_1)[1] <- 'neuron_name'

colnames(basalaxon_1)
dim(basalaxon_1)
basalaxon_1 <- rmvExtnsn(basalaxon_1)
basalaxon_1 <- subset(basalaxon_1, dpos==FALSE)

basalaxon_1$neuron_name
colnames(basalaxon_1)
basalaxonreport <- addNeuronID(masterMetaData[,1:2], basalaxon_1)
dim(basalaxonreport)
#write the axon Only 
write.table(basalaxonreport, file="basalaxon_v5.6.txt",quote=FALSE, sep="\t",row.names=FALSE)