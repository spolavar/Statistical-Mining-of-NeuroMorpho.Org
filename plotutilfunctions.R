#the initial set of metrics

#whole arbor pca metrics
pcawholeMetrics <- c("neuron_name","Soma_Surface_total_sum","N_bifs_total_sum","N_tips_total_sum","Width_total_sum","Height_total_sum","Depth_total_sum","Length_total_sum","EucDistance_avg","EucDistance_max","PathDistance_avg","PathDistance_max","Branch_Order_max","Terminal_degree_avg","Branch_pathlength_avg","Branch_pathlength_max","Contraction_avg","Partition_asymmetry_avg","Bif_ampl_local_avg","Bif_ampl_local_max","Bif_ampl_remote_avg","Bif_ampl_remote_max","Bif_tilt_local_avg","Bif_tilt_local_max","Bif_tilt_remote_avg","Bif_tilt_remote_max","Bif_torque_local_avg", "Bif_torque_local_max","Bif_torque_remote_avg","Bif_torque_remote_max","Helix_max","Fractal_Dim_avg","Fractal_Dim_max")
length(pcawholeMetrics)

#mypalette <- c("magenta","orange","cyan","green","blue","seagreen4","darkorchid","darkcyan","darkmagenta","yellow","black","darkred","violetred4", "pink","darkgreen","darkgrey")
mypalette <- c("red","cyan","orange","black","blue","magenta","green","seagreen4","darkorchid","darkcyan","darkmagenta","yellow","darkred","violetred4", "pink","darkgreen","darkgrey")
alphabets <- c('a','b','c','d','e','f','g','h','i')
zscralphconvers <- c("f","e","a","c","d","b")
rgblim <- 200
rgbrange <- c(0:255)
colors()[c(552,254,26)]
Rcolorchart_100 <- colors()[1:100]
length(Rcolorchart_100)
###test###
test <- c("zero")
test <- append(test, mypalette)
test
mypalette
any(mypalette=="darkred")
#choose 50 random colors after eliminating while shaded ones
mypalette_50 <- c()
cnt <- 50
while(cnt>=1){
  rcol <- sample(0:255, 1)
  gcol <- sample(0:255, 1)
  bcol <- sample(0:255, 1)
  rgbcol <- rgb(rcol, gcol, bcol,maxColorValue=255)
  #make sure that all three rgb values are not >200 at the same time. this will avoid pale white like colors.
  if(rcol>200 & gcol>200 & bcol>200){
    print("skipping color..")
  }else if(!any(mypalette_50==rgbcol)){
    cnt <- cnt - 1
    mypalette_50 <- append(mypalette_50,rgbcol)
  }
}
mypalette_50
plot(x=1:length(mypalette_50),y=rep(1,length(mypalette_50)),pch=19,col=mypalette_50)

mypalette_50 <- colors()[c(76:125, 451:550)]
whitecolgrp <- grep("white",colors())
whitecolgrp
lightcolgrp
probofsampling <- seq(1,50)

sample(76:125, 50, replace=F)
sample(76:125, 50, replace=F, prob=)

colorcode <- function(groupnum, order){
  if(!missing(order))
    colcode <- order
  else
    colcode <- unique(groupnum)
    
  colarr <- c()
  #loop and assign colors to each data point 
  for(i in 1:length(colcode)){
    colarr[groupnum==colcode[i]] <- i
  }
  
  #plot(x=1:16,y=rep(1,16),pch=19,col=mypalette)
  
  unique(colarr)
  return (colarr)
}

alphabetcode <- function(groupnum){
  groupnum <- zscralphconvers[groupnum]
  return(groupnum)
}

  #DendReduced[,"alphclusters"] <- zscralphconvers[DendReduced$classification]
  #alpharr <- c()
  #if(!missing(order)){
    
  #}else
   # for(i in 1:length(groupnum)){
    #  alpharr[groupnum==order[i]] <- i
    #}
  #loop and assign indices to each data point 
  #unique(alpharr)
  #return (alpharr)
#}

#Density function
kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=50,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,       # see option nlevels in contour
                      theta=30,         # see option theta in persp
                      phi=30)           # see option phi in persp
{
  z   <- d$z
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
  fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
  dim(fcol) <- c(nrz,ncz)
  fcol      <- fcol[-nrz,-ncz]
  
  par(mfrow=c(1,2),mar=c(0.5,0.5,0.5,0.5))
  persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,zlab="density")
  
  par(mar=c(2,2,2,2))
  image(d,col=couleurs)
  contour(d,add=T,nlevels=nlevels)
  box()
  par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
}

#compute distance between two multidimensional vectors
pairwiseMDdist <- function(p,q){
  S <- 0
  for(i in 1:length(p)){
    S <- S +(p[i]-q[i])^2
  }
  dist <- sqrt(S)
  return (dist)
}


editBrainRegionLabels <- function(tmpDf,colname){
  tmpDf[,colname] <- as.character(tmpDf[,colname])
  #replace the groupttl column for brain regions with legend titles
  tmpDf[tmpDf[,colname]== "Control/Hippocampus/Hippocampus" ,colname] <-  "Hippocampus"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory" ,colname] <-  "Neocortex/Parietal/S1"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory/Barrel" ,colname] <-  "Neocortex/Parietal/S1/Barrel"
  tmpDf[tmpDf[,colname]== "Control/Hippocampus/Hippocampus/CA1" ,colname] <-  "Hippocampus/CA1"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory/Others"  ,colname] <-  "Other Neocortex/Parietal/S1"
  tmpDf[tmpDf[,colname]== "Control/Retina/Retina",colname] <-  "Retina"
  tmpDf[tmpDf[,colname]== "Control/Retina/Retina/Ganglion cell layer",colname] <-  "Retina/Ganglion cell layer"
  tmpDf[tmpDf[,colname]== "Control/Hippocampus/Hippocampus/Others",colname] <-  "Other Hippocampus"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory/Layer 4",colname] <- "Neocortex/Parietal/S1/L4"
  tmpDf[tmpDf[,colname]== "Control/Olfactory bulb/Olfactory bulb",colname] <-  "Olfactory bulb"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Others",colname] <-"Other Neocortex/Parietal"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Occipital lobe",colname] <-  "Neocortex/Occipital"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Prefrontal",colname] <-   "Neocortex/Frontal/Prefrontal"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Motor",colname] <- "Neocortex/Frontal/Motor"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Others",colname] <- "Other Neocortex/Frontal"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Insula",colname] <-  "Neocortex/Insula"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Others",colname] <-   "Other Neocortex"
  tmpDf[tmpDf[,colname]== "Control/Others",colname] <-  "Other Brain regions"
  print(unique(tmpDf[,colname]))
  return (tmpDf)
}

editCelltypeLabels <- function(tmpDf,colname){
  tmpDf[,colname] <- as.character(tmpDf[,colname])
  #replace the groupttl column for cell types with legend titles
 
  tmpDf[tmpDf[,colname]== "Control/Interneuron/Martinotti cell" ,colname] <-  "Martinotti"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Pyramidal cell/Others" ,colname] <-  "Other Pyramidal"
  tmpDf[tmpDf[,colname]== "Control/Interneuron/Others"  ,colname] <-  "Other Interneurons"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Ganglion cell/Others",colname] <-  "Other Ganglion"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Others",colname] <-  "Other Principal"
  tmpDf[tmpDf[,colname]== "Control/Interneuron/Basket cell",colname] <-  "Basket"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Granule cell",colname] <-  "Granule"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Medium spiny cell",colname] <-  "Medium spiny"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Pyramidal cell/Magnopyramidal",colname] <-  "Magnopyramidal"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Motoneuron",colname] <-  "Motoneuron"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Ganglion cell/Monostratified",colname] <-  "Monostratified"
  tmpDf[tmpDf[,colname]== "Control/Axonal terminal/Uniglomerular projection neuron/iPN",colname] <-  "AT/Uniglomerular/iPN"
  tmpDf[tmpDf[,colname]== "Control/Principal cell/Pyramidal cell",colname] <-  "Pyramidal"
  tmpDf[tmpDf[,colname]== "Control/Axonal terminal/Others",colname] <-  "Other AT"
  tmpDf[tmpDf[,colname]== "Control/Axonal terminal/Uniglomerular projection neuron/Others",colname] <-  "Other Uniglomerular PNs"
  print(unique(tmpDf[,colname]))
  return (tmpDf)
}

editspeciesLabels <- function(tmpDf,colname){
  #eliminate NA rows
  tmpDf <- subset(tmpDf,!is.na(tmpDf[,colname]))
  tmpDf[,colname] <- as.character(tmpDf[,colname])
  #replace the groupttl column for brain regions with legend titles
  tmpDf[tmpDf[,colname]== "Control/Rodent/Rat" ,colname] <-  "Rodent/Rat"
  tmpDf[tmpDf[,colname]== "Control/Rodent/Mouse" ,colname] <-  "Rodent/Mouse"
  tmpDf[tmpDf[,colname]== "Control/Primate/Human",colname] <-   "Primate/Human"
  tmpDf[tmpDf[,colname]== "Control/Primate/Monkey",colname] <-  "Primate/Monkey"
  tmpDf[tmpDf[,colname]== "Control/Carnivora/Cat",colname] <-  "Carnivora/Cat"
  tmpDf[tmpDf[,colname]== "Control/Proboscidae/Elephant",colname] <-  "Elephant"
  tmpDf[tmpDf[,colname]== "Control/Bony fishes/Goldfish",colname] <-  "Goldfish"
  tmpDf[tmpDf[,colname]== "Control/Others",colname] <-  "Other species"
  tmpDf[tmpDf[,colname]== "Control/Rodent/Others",colname] <-  "Other Rodents"
  tmpDf[tmpDf[,colname]== "Control/Insects/Blowfly" ,colname] <-  "Insects/Blowfly"  
  tmpDf[tmpDf[,colname]== "Control/Amphibians/Salamander" ,colname] <-  "Amphibians/Salamander"
  tmpDf[tmpDf[,colname]== "Control/Insects/Drosophila" ,colname] <-  "Insects/Drosophila"
  tmpDf[tmpDf[,colname]== "Control/Insects/Cricket" ,colname] <-  "Insects/Cricket"
  
  print(unique(tmpDf[,colname]))
  return (tmpDf)
}

addmocklabels <- function(tmpDf,colname){
  tmpDf[tmpDf[,colname]== "Control/Hippocampus/Hippocampus/CA1" ,colname] <-  "CA1"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory/Others"  ,colname] <-  "O S1"
  tmpDf[tmpDf[,colname]== "Control/Retina/Retina",colname] <-  "R"
  tmpDf[tmpDf[,colname]== "Control/Hippocampus/Hippocampus/Others",colname] <-  "O H"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Somatosensory/Layer 4",colname] <- "L4"
  tmpDf[tmpDf[,colname]== "Control/Olfactory bulb/Olfactory bulb",colname] <-  "OB"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Parietal lobe/Others",colname] <-"O P"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Occipital lobe",colname] <-  "Occ"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Prefrontal",colname] <-   "Pr"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Motor",colname] <-   "Mr"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Frontal lobe/Others",colname] <-   "O Fr"
  tmpDf[tmpDf[,colname]== "Control/Neocortex/Others",colname] <-   "O N"
  tmpDf[tmpDf[,colname]== "Control/Others",colname] <-  "O Br"
  print(unique(tmpDf[,colname]))
  return (tmpDf)
}

cores <- function() {
  
  par(mar=c(0,0,0,0),mgp=c(0,0,0))
  plot(c(0:24),type='n') 
  c <- 0
  
  mouse <- function(b, x, y) {
    x <- as.integer(x*26)
    y <- as.integer(y*26)
    print(colors()[(x+26*y) %% 657 + 1])
    return()
  } 
  k <- colours()[(1:26^2 - 1) %% 657 + 1]
  for (i in 1:26) {
    for (j in 1:26) {
      c <- c+1
      polygon(c(j,j,j-1,j-1),c(i,i-1,i-1,i)-1,col=k[(c-1) %% 657 + 1])
    }
  }
  getGraphicsEvent('Click on a colour!',onMouseDown=mouse)
  
}


#plot cluster
clustercolorcode <- function(subgrpDf,grpnum){
  colcode <- unique(subgrpDf[,grpnum])
  colarr <- c()
  #loop and assign colors to each data point for max 15 colors
  for(i in 1:15){
    colarr[subgrpDf[,grpnum]==colcode[i]] <- i
  }
  print(unique(colarr))
  return(colarr)
}
 
#make 4D groupnum
addmixedgrpcol <- function(mixedgrpsDf, hybrid3pcaDf){
  #the data frame to be returned is first initialized tohybrid3pcaDf
  resDf <- data.frame()
  #assign groupnums to 
  for(i in 1:nrow(mixedgrpsDf)){
    print(paste(mixedgrpsDf[i,"archive_name"], mixedgrpsDf[i,"ct_groupttl"],mixedgrpsDf[i,"br_groupttl"],mixedgrpsDf[i,"sp_groupttl"]))
    #condition to select a specific subgroup of hybrid3pcaDf
    mxdcond <-hybrid3pcaDf$archive_name == mixedgrpsDf[i,"archive_name"] & 
     hybrid3pcaDf$ct_groupttl == mixedgrpsDf[i,"ct_groupttl"] &
     hybrid3pcaDf$br_groupttl == mixedgrpsDf[i,"br_groupttl"] &
     hybrid3pcaDf$sp_groupttl == mixedgrpsDf[i,"sp_groupttl"]
    #print(mxdcond)
    tmp <- data.frame()  
    tmp <- subset(hybrid3pcaDf, rownames(hybrid3pcaDf)%in%rownames(subset(hybrid3pcaDf, mxdcond)))
    #tmp2 <- subset(hybrid3pcaDf[,c(1:32,65:70)], rownames(dendriteReduced)%in%rownames(subset(hybrid3pcaDf, mxdcond)))
    if(length(tmp)>0){
      tmp["mxdgroupnum"] <- i #add the mxdgroupnum
      #tmp <- cbind(tmp,tmp2) #add metadata to PCA matrix
    }
    resDf <- rbind(resDf,tmp)
    
  }
  print(dim(resDf))
  print(colnames(resDf))
  print(unique(resDf$mxdgroupnum))
  
  return(resDf)
}