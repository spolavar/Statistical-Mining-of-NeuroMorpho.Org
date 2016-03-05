#fly data
f1 <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_SomaOnly.txt"
f2 <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_wholeArbor.txt"
f1_df <- read.delim(f1, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
f2_df <- read.delim(f2, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(f1_df)[1] <- 'neuron_name'
colnames(f2_df)[1] <- 'neuron_name'
f1_df <- rmvExtnsn(f1_df)
f2_df <- rmvExtnsn(f2_df)
f1_df <- f1_df[, !(colnames(f1_df) %in% c("X"))]
f2_df <- f2_df[, !(colnames(f2_df) %in% c("X"))]

dim(f1_df)
dim(f2_df)

flyMetricData <- merge(x=f1_df,y=f2_df, by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
dim(flyMetricData)
colnames(flyMetricData)
write.table(flyMetricData, file="flyWholearbor.txt",quote=FALSE, sep="\t",row.names=F)

#Basal only file
fb <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_BasalOnly.txt"
fb <- read.delim(fb, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(fb)[1] <- 'neuron_name'
fb <- rmvExtnsn(fb)
fb <- fb[, !(colnames(fb) %in% c("X"))]
dim(fb)
fb$neuron_name
colnames(fb)
write.table(fb, file="flyBasalOnly.txt",quote=FALSE, sep="\t",row.names=F)

#Basal axon 
f <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_BasalAxon.txt"
fba <- read.delim(f, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(fba)[1] <- 'neuron_name'
fba <- rmvExtnsn(fba)
fba <- fba[, !(colnames(fba) %in% c("X"))]
dim(fba)
fba$neuron_name
colnames(fba)
write.table(fba, file="flyBasalAxon.txt",quote=FALSE, sep="\t",row.names=F)

#Axon only
f <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_AxonOnly.txt"
fa <- read.delim(f, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(fa)[1] <- 'neuron_name'
fa <- rmvExtnsn(fa)
fa <- fa[, !(colnames(fa) %in% c("X"))]
dim(fa)
fa$neuron_name
colnames(fa)
write.table(fa, file="flyAxonOnly.txt",quote=FALSE, sep="\t",row.names=F)

#apical only
f <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_ApicalOnly.txt"
fap <- read.delim(f, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(fap)[1] <- 'neuron_name'
fap <- rmvExtnsn(fap)
fap <- fap[, !(colnames(fap) %in% c("X"))]
dim(fap)
fap$neuron_name
colnames(fap)
write.table(fap, file="flyApicalOnly.txt",quote=FALSE, sep="\t",row.names=F)

#apical axon 
f <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_ApicalAxon.txt"
faa <- read.delim(f, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(faa)[1] <- 'neuron_name'
faa <- rmvExtnsn(faa)
faa <- faa[, !(colnames(faa) %in% c("X"))]
dim(faa)
faa$neuron_name
colnames(faa)
write.table(faa, file="flyApicalAxon.txt",quote=FALSE, sep="\t",row.names=F)

#apical basal
f <- "Z:/LMAutoTesting/MetricExtraction/MetricExtractionReports/Fly_Reports/Report_ApicalBasal.txt"
fab <- read.delim(f, header = TRUE, sep="\t", dec=".",stringsAsFactors=FALSE)
colnames(fab)[1] <- 'neuron_name'
fab <- rmvExtnsn(fab)
fab <- fab[, !(colnames(fab) %in% c("X"))]
dim(fab)
fab$neuron_name
colnames(fab)
write.table(fab, file="flyApicalBasal.txt",quote=FALSE, sep="\t",row.names=F)


