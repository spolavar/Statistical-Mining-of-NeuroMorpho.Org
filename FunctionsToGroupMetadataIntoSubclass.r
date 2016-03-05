# Meta-data exploration: explore groups and their subgroups while size is greater than some minimum
# Create order column from species
masterMetaData$order <- rep("",nrow(masterMetaData))
masterMetaData[masterMetaData$species%in%c("Human","Monkey"),"order"] <- "Primate"
masterMetaData[masterMetaData$species%in%c("Rat","Mouse","Guinea-pig", "Agouti","Proechimys"),"order"] <- "Rodent"
masterMetaData[masterMetaData$species=="Elephant","order"] <- "Proboscidae"
masterMetaData[masterMetaData$species%in%c("Blowfly","Drosophila","Cricket"),"order"] <- "Insects"
masterMetaData[masterMetaData$species%in%c("Goldfish","Zebrafish"),"order"] <- "Bony fishes"
masterMetaData[masterMetaData$species=="Rabbit","order"] <- "Lagomorpha"
masterMetaData[masterMetaData$species=="Cat","order"] <- "Carnivora"
#masterMetaData[masterMetaData$species=="Cricket","order"] <- "Orthoptera"
masterMetaData[masterMetaData$species%in%c("Salamander","Frog"),"order"] <- "Amphibians"
masterMetaData[masterMetaData$species=="Spiny lobster","order"] <- "Crustaceans"
#masterMetaData[masterMetaData$species=="Frog","order"] <- "Anura"
masterMetaData[masterMetaData$species=="Turtle","order"] <- "Reptiles"
masterMetaData[masterMetaData$species=="Chicken","order"] <- "Birds"
masterMetaData[masterMetaData$species=="C. elegans","order"] <- "Nematodes"

aggregate(masterMetaData$order,by=list(masterMetaData$order),length)
#cortical data
cortical <- masterMetaData$region1 == "Neocortex"

# Create Layer column based on Region3
masterMetaData$layer <- masterMetaData$region3
#move the layer specific metadata from region3 to layer
layerInds <- with(masterMetaData,grepl("layer",region3,ig=T) | grepl("Stratum",region3,ig=T) |
  grepl("granul",region3,ig=T) | grepl ("Hilus",region3,ig=T) | grepl("Superficial",region3,ig=T) | grepl("Deep",region3,ig=T) | grepl("Middle",region3,ig=T))

unique(masterMetaData$layer[layerInds])
unique(masterMetaData$layer[!layerInds])
masterMetaData[!layerInds,"layer"] <- "Not reported"
unique(masterMetaData$layer)

#move layer specific values from region2 to layer
layerInds <- with(masterMetaData,grepl("layer",region2,ig=T))
unique(masterMetaData$layer[layerInds])
unique(masterMetaData$region2[layerInds])
masterMetaData[layerInds,"layer"] <- masterMetaData[layerInds,"region2"]
unique(masterMetaData$layer)

#move layer specific values from region3b to layer
layerInds <-  with(masterMetaData,grepl("layer",region3b,ig=T) | grepl("Stratum",region3b,ig=T) |
  grepl("granul",region3b,ig=T) | grepl ("Hilus",region3b,ig=T) | grepl("Superficial",region3b,ig=T) | grepl("Deep",region3b,ig=T) | grepl("Middle",region3b,ig=T))

unique(masterMetaData$layer[layerInds])
#unique(masterMetaData[layerInds,c("region2","region3","layer","region3b")])
masterMetaData[layerInds,"layer"] <- masterMetaData[layerInds,"region3b"]
unique(masterMetaData$layer)

aggregate(masterMetaData$layer,by=list(masterMetaData$layer),length)

#merge Drosophila 'Not reported' strain with 'melanogaster'
pos <- with(masterMetaData,species =="Drosophila")
masterMetaData[pos, "strain"] <- "melanogaster"

# group region1 areas under a superclass
#follow the metadata hierarchy to group region1,2 and 3 as single hierarachy
#masterMetaData[masterMetaData$region1%in%c("Dorsal thalamus","Ventral thalamus"),"region1"] <- "Brainstem"
masterMetaData[masterMetaData$region1%in%c("Midbrain","Cerebellum","Medulla"),"region1"] <- "Brainstem"

masterMetaData[masterMetaData$region2%in%c("Cochlear nucleus","Hypoglossal nucleus","Nucleus laminaris (NL)"),"region1"] <- "Brainstem"

# replace region1 "olfactory bulb" of antennal lobe with "deuterocerebrum"
pos <- with(masterMetaData,region2 =="Antennal lobe")
masterMetaData[pos,"region1"] <- "Deuterocerebrum"

# join Anterior olfactory nucleus with olfactory bulb in region1
masterMetaData[masterMetaData$region1 == "Anterior olfactory nucleus","region1"] <- "Olfactory bulb"

#move all 'Entorhinal' in region2 to region1 exclusively
masterMetaData[masterMetaData$region2%in%c("Entorhinal"),"region1"] <- "Hippocampus"

# Create lobes column based on region2 of Neocortex
pos <- cortical & grepl('lobe',masterMetaData$region2)
length(pos[pos==TRUE])
unique(masterMetaData$region2[pos])
#( || grep('Insula',masterMetaData$region2))
unique(masterMetaData$region2[pos])
masterMetaData[pos, "lobes"] <- masterMetaData$region2[pos]
masterMetaData[!pos,"lobes"] <- "Not reported"
#add non-neocortical region to lobes column
pos <- !cortical#masterMetaData$region1 != "Neocortex"  
masterMetaData[pos, "lobes"] <- masterMetaData$region1[pos]
unique(masterMetaData$lobes)
aggregate(masterMetaData$lobes,by=list(masterMetaData$lobes),length)

# merge Middle short insular gyrus,Anterior long insular gyrus,Posterior short insular gyrus and fronto-insula in region2 with insula in region2
masterMetaData[masterMetaData$region2%in%c("Middle short insular gyrus","Anterior long insular gyrus","Posterior short insular gyrus","Fronto-insula","Insular cortex"),"lobes"] <- "Insula"

# merge anterior cingulate with limbic lobe in region2
masterMetaData[masterMetaData$region2%in%c("Anterior cingulate","Perirhinal"),"lobes"] <- "Limbic lobe"

#merge visual with occipital lobe in region2
masterMetaData[masterMetaData$region2%in%c("Visual","Suprasylvian gyrus"),"lobes"] <- "Occipital lobe"
masterMetaData[masterMetaData$region3%in%c("Area 17","Brodmann area 18"),"lobes"] <- "Occipital lobe"


#merge Primary auditory with temporal lobe in region2
masterMetaData[masterMetaData$region2%in%c("Primary auditory"),"lobes"] <- "Temporal lobe"
masterMetaData[masterMetaData$region3%in%c("Brodmann area 22"),"lobes"] <- "Temporal lobe"

#merge Postcentral gyri, Brodman area 3,1,2 with Somatosensory and Somatosensory with Temporal lobe 
masterMetaData[masterMetaData$region2%in%c("Postcentral gyri"),"region2"] <- "Somatosensory"
masterMetaData[masterMetaData$region2%in%c("Precentral gyri"),"region2"] <- "Somatosensory"
masterMetaData[masterMetaData$region3%in%c("Brodmann area 3,1,2"),"region2"] <- "Somatosensory"
masterMetaData[masterMetaData$region2%in%c("Somatosensory"),"lobes"] <- "Parietal lobe"
masterMetaData[masterMetaData$region3%in%c("Brodmann area 39"),"lobes"] <- "Parietal lobe"


#merge Temporal sulcus with Parietal lobe
masterMetaData[masterMetaData$region2%in%c("Temporal sulcus"),"lobes"] <- "Parietal lobe"

#testing...
#pos <- masterMetaData$lobes == "Parietal lobe"
#aggregate(pos,by=list(pos),length)

#merge region2 and region3 areas into Motor, which in turn is merged to Frontal lobe
masterMetaData[masterMetaData$region3%in%c("Brodmann area 6","Brodmann area 4","M1/M2"),"region2"] <- "Motor"
masterMetaData[masterMetaData$region2%in%c("Precentral gyri"),"region2"]<- "Motor"
masterMetaData[masterMetaData$region2%in%c("Orbital cortex","Motor"),"lobes"] <- "Frontal lobe"

#merge region2 and region3 areas with prefrontal
masterMetaData[masterMetaData$region3%in%c("Brodmann area 44","Brodmann area 10","Brodmann area 11","Inferior frontal gyrus","Superior frontal gyrus","Prelimbic","Infralimbic"),"region2"] <- "Prefrontal"
masterMetaData[masterMetaData$region2%in%c("Frontopolar","Medial prefrontal cortex"),"region2"] <- "Prefrontal"
masterMetaData[masterMetaData$region2%in%c("Prefrontal"),"lobes"] <- "Frontal lobe"


unique(masterMetaData$cellclass2)
#combine parvalbumin and basket cells; combine Martinotti and SOM containing cell
masterMetaData[cortical&masterMetaData$cellclass2=="Parvalbumin (PV) containing cell","cellclass2"] <- "Basket cell"
masterMetaData[cortical&masterMetaData$cellclass2=="Somatostatin (SOM) containing cell","cellclass2"] <- "Martinotti cell"
#combine long projecting pyramidal cells with Ipsilateral-projecting (original category has metadata error)
unique(masterMetaData$cellclass3)
masterMetaData[masterMetaData$cellclass3=="Local projecting","cellclass3"] <- "Ipsilateral-projecting"

#combine C57BL/6, C57BL/6J into C57BL6 Mouse in strain
masterMetaData[masterMetaData$strain%in%c("C57BL/6J"," C57BL/6J"),"strain"] <- "C57BL6 Mouse"

#combine C57BL6/Thy1-GFP Mouse;GIN;thy-1-YFP-16;Dbx1cre;ROSA26 YFP;Nkx2.1Cre MADM;G42;Somatostatin-GFP Mouse;Calretinin-EGFP;5HT3-EGFP into "Transgenic mice"
masterMetaData[masterMetaData$strain%in%c("C57BL6/Thy1-GFP Mouse","GIN","thy-1-YFP-16","Dbx1cre;ROSA26 YFP","Nkx2.1Cre MADM","G42","Somatostatin-GFP Mouse","Calretinin-EGFP","5HT3-EGFP","rd10/Thy1-GFP","CXCR4-EGFP"),"strain"] <- "Transgenic mice"

# combine C57BL6/TrkB.T1 deficient mouse with "Knockout mice"
masterMetaData[masterMetaData$strain%in%c("C57BL6/TrkB.T1 deficient mouse"),"strain"] <- "Knockout mice"

# group region1 areas under a superclass
#brainstemInds <- with(masterMetaData,grepl("rodent",Order,ig=T) | grepl("primate",Order,ig=T))| grepl("carnivora",Order,ig=T))

#masterMetaData[masterMetaData$region1%in%c("Dorsal thalamus","Ventral thalamus"),"region1"] <- "Brainstem"

#masterMetaData[masterMetaData$region1 == "Anterior olfactory nucleus","region1"] <- "Olfactory bulb"

#region2 grouping
#need to ask how to achive this grouping alternatively with grepl front
#masterMetaData[masterMetaData$region2%in%c("Frontopolar","Prefrontal","Medial prefrontal cortex","Fronto-insula", "Precentral gyri", "Inferior frontal gyrus"),"region2"] <- "Frontal lobe"

#masterMetaData[masterMetaData$region2 == "Temporal sulcus","region2"] <- "Parietal lobe"

#masterMetaData[masterMetaData$region2 == "Visual","region2"] <- "Occipital lobe"

#masterMetaData[masterMetaData$region2%in%c("Anterior cingulate","Perirhinal"),"region2"] <- "Limbic lobe"


# Returns the subset of a data frame given a subtype specification
getSubset2 <- function(mydata,specifications=NULL,typeMetaData=masterMetaData,minLength=NaN,maxLength=NaN){
  grpcol <- typeMetaData$neuron_name
  #print(length(grpcol))
  if (!is.null(specifications)){
    for (i in 1:length(specifications)){
      specType <- names(specifications)[i]
      #print(specType)
      if (is.list(specifications)){
        #print(specType)
        specValue <- specifications[[i]]
      }
      else{
        specValue <- specifications[i]
      }
      #cat(specType,":",specValue,":",length(specValue))
      if (length(specValue) > 1){
        print("should'nt enter here...")
        thegrpcol <- c()
        for (j in specValue){
          thegrpcol <- union(thegrpcol,typeMetaData[typeMetaData[,specType]==j,"neuron_name"])
        }
        grpcol <- intersect(grpcol,thegrpcol)
      }
      else{
        #flush.console()
        grpcol <- intersect(grpcol,typeMetaData[typeMetaData[,as.character(specType)]==as.character(specValue),"neuron_name"])
        #print(length(grpcol))
      }
    }
  }
  #print(minLength)
  #if (!is.nan(minLength)){
  #  size <- intersect(size,as.character(subset(typeMetaData,SeqLen>=minLength)$n_bifs))
  #}
  #print(maxLength)
  #if (!is.nan(maxLength)){
  #  size <- intersect(size,as.character(subset(typeMetaData,SeqLen<=maxLength)$n_bifs))
  #}
  #if ("neuron_name"%in%colnames(mydata)){
  #print("returning..")
  #ret <- subset(mydata,neuron_name%in%grpcol)
  #dim(ret)
  #subset of rows that satisfy specification list
  returndf <- subset(mydata,neuron_name%in%grpcol)
  #add other column for small groups
  #if(length(grpcol) < minLength){
    #create a single col dataframe named "OtherCol"
  #  otherdf <- data.frame("otherCol"='Other')
    #add "other" as value for each neuron_name
  #  otherdf <- do.call(rbind, rep(list(otherdf), length(grpcol)))
  #  returndf <- cbind(returndf,otherdf) 
  #  print("added otherCol..")
  #}
  
  
  
  return(returndf)
  #return(mydata[grpcol,])  
  #}
  #else{
  #  size <- intersect(size,rownames(mydata))
  #  return(mydata[size,])
  #}
}

#Testing function call
#primary <- list(protocol = c("In vivo","In vitro","Not reported","Culture"))
#par(ask=TRUE)
#x <- getSubset2(masterMetaData, primary, typeMetaData=masterMetaData)
#dim(x)
#colnames(x)