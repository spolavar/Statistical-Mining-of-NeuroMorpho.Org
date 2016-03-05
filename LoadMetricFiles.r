#open metadata file of v5.4 
# Load meta data
# Set metadata filename
#metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/all_7986_v5.4_metadata_2.csv"
#metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/all_7986_v5.4_metadata_3.csv"
#metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/all_8858_v5.5_metadata.csv"
metadataFilename <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/NMO_DB_metadata.csv"
masterMetaData <- read.csv(metadataFilename,header = TRUE, sep="}",na.strings=c("NA",""),quote="\"",stringsAsFactors=FALSE)
print("finished reading...")
colnames(masterMetaData)
unique(masterMetaData$order)
#removing duplicate rows
dpos <- duplicated(masterMetaData)
length(dpos[dpos==TRUE])
masterMetaData <- masterMetaData[, !(colnames(masterMetaData) %in% c("X"))]

#read the part with only metrics for the whole arbor (with specificity type > 1)
#metricWholeFile_1 <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/whole.txt"
metricWholeFile_1 <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_wholeArbor_v5.6.txt"
#read the part with soma_surface
#metricWholeFile_2 <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/somawhole.txt"
metricWholeFile_2 <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_somaOnly_v5.6.txt"

metricDataWhole_1 <- read.delim(metricWholeFile_1, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
#rename the column name 
colnames(metricDataWhole_1)[1] <- 'neuron_name'

colnames(metricDataWhole_1)
dim(metricDataWhole_1)
metricDataWhole_1 <- rmvExtnsn(metricDataWhole_1)
metricDataWhole_1$neuron_name
colnames(metricDataWhole_1)
#making sure no extraneous columns like "x" and "Soma_Surface_total_sum" are not in the first part
metricDataWhole_1 <- metricDataWhole_1[, !(colnames(metricDataWhole_1) %in% c("X","Soma_Surface_total_sum"))]


metricDataWhole_2 <- read.delim(metricWholeFile_2, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(metricDataWhole_2)[1] <- 'neuron_name'
colnames(metricDataWhole_2)
metricDataWhole_2 <- rmvExtnsn(metricDataWhole_2)

#remove extraneous columns
metricDataWhole_2 <- metricDataWhole_2[, !(colnames(metricDataWhole_2) %in% c("X","N_stems_total_sum"))]
colnames(metricDataWhole_2)
dim(metricDataWhole_1)
dim(metricDataWhole_2)
#merge all metrics over neuron_name
t <- setdiff(metricDataWhole_1$neuron_name,metricDataWhole_2$neuron_name)
length(t)
#tt <- aggregate(metricDataWhole_1$neuron_name,by=list(metricDataWhole_1$neuron_name),length)
#tt[order(tt$x,decreasing=T),]
#colnames(tt)
#head(tt,n=15)
#nrow(tt)
nrow(metricDataWhole_1)
#removing duplicate rows
dpos <- duplicated(metricDataWhole_1)
metricDataWhole_1 <- subset(metricDataWhole_1, dpos==FALSE)
dim(metricDataWhole_1)
dpos <- duplicated(metricDataWhole_2)
metricDataWhole_2 <- subset(metricDataWhole_2, dpos==FALSE)
dim(metricDataWhole_2)
metricDataWhole <- merge(x=metricDataWhole_2,y=metricDataWhole_1, by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(metricDataWhole)
print("finished reading whole...")
#checking the neurons that are present in masterMetaData but not in metricDataWhole
setdiff(masterMetaData$neuron_name,metricDataWhole$neuron_name)
#files that are flat (<1um) should have depth as NA
depthpos <- with(metricDataWhole,which(Depth_total_sum<1))
metricDataWhole[depthpos,"Depth_total_sum"] <- NA
colnames(metricDataWhole)
metricDataWhole$Depth_total_sum
metricDataWhole$Soma_Surface_total_sum
dim(metricDataWhole)

#v5.5 neurons without soma
somaNA <- subset(metricDataWhole,is.na(metricDataWhole$Soma_Surface_total_sum))
dim(somaNA)
#open apical dendrites only
#metricApicalFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/apical.txt"
metricApicalFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_apicalOnly.txt"
metricDataApical <- read.delim(metricApicalFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(metricDataApical)[1] <- 'neuron_name'
print("finished reading Apical...")
metricDataApical <- metricDataApical[, !(colnames(metricDataApical) %in% c("X","Soma_Surface_total_sum"))]
colnames(metricDataApical)
#removing duplicate rows
dpos <- duplicated(metricDataApical)
metricDataApical <- subset(metricDataApical, dpos==FALSE)
dim(metricDataApical)
#files that are flat (<1um) should have depth as NA
depthpos <- with(metricDataApical,which(Depth_total_sum<1))
metricDataApical[depthpos,"Depth_total_sum"] <- NA


#open apical & basal dendrites only
#metricApiBasFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/apicalbasal.txt"
metricApiBasFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_apicalBasal_v5.6.txt"
metricDataApiBas <- read.delim(metricApiBasFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(metricDataApiBas)[1] <- 'neuron_name'
print("finished reading apical & basal...")
metricDataApiBas <- metricDataApiBas[, !(colnames(metricDataApiBas) %in% c("X","Soma_Surface_total_sum"))]
colnames(metricDataApiBas)
#removing duplicate rows
dpos <- duplicated(metricDataApiBas)
metricDataApiBas <- subset(metricDataApiBas, dpos==FALSE)
dim(metricDataApiBas)
#files that are flat (<1um) should have depth as NA
depthpos <- with(metricDataApiBas,which(Depth_total_sum<1))
metricDataApiBas[depthpos,"Depth_total_sum"] <- NA

#open basal dendrites only
#metricDendFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/basal.txt"
metricDendFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_basalOnly.txt"
metricDataDend <- read.delim(metricDendFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(metricDataDend)[1] <- 'neuron_name'
print("finished reading Dend...")
metricDataDend <- metricDataDend[, !(colnames(metricDataDend) %in% c("X","Soma_Surface_total_sum"))]
colnames(metricDataDend)
#removing duplicate rows
dpos <- duplicated(metricDataDend)
metricDataDend <- subset(metricDataDend, dpos==FALSE)
dim(metricDataDend)
#files that are flat (<1um) should have depth as NA
depthpos <- with(metricDataDend,which(Depth_total_sum<1))
metricDataDend[depthpos,"Depth_total_sum"] <- NA
dim(metricDataDend)

#open axons only
#metricAxonFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/metricData/axon.txt"
metricAxonFile <- "C:/Users/sridevi/thesismaterial/proposal/LMpaper/Rcode/Report_axonOnly.txt"
metricDataAxon <- read.delim(metricAxonFile, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(metricDataAxon)[1] <- 'neuron_name'
print("finished reading Axon...")
metricDataAxon <- metricDataAxon[, !(colnames(metricDataAxon) %in% c("X","Soma_Surface_total_sum"))]
colnames(metricDataAxon)
#removing duplicate rows
dpos <- duplicated(metricDataAxon)
metricDataAxon <- subset(metricDataAxon, dpos==FALSE)
dim(metricDataAxon)
#files that are flat (<1um) should have depth as NA
depthpos <- with(metricDataAxon,which(Depth_total_sum<1))
metricDataAxon[depthpos,"Depth_total_sum"] <- NA
dim(metricDataAxon)