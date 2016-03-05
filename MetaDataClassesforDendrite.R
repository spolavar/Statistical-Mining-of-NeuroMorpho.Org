
colnames(metricDataDendrite)
reducedMetrics <- c("neuron_name","Soma_Surface_total_sum","N_bifs_total_sum", "N_tips_total_sum","Width_total_sum","Height_total_sum","Depth_total_sum","Length_total_sum","EucDistance_max","PathDistance_max","Branch_Order_max", "Contraction_avg", "Partition_asymmetry_avg","Bif_ampl_local_avg", "Fractal_Dim_avg")
pcaMetrics <- c("neuron_name","Soma_Surface_total_sum","N_bifs_total_sum","N_tips_total_sum","Width_total_sum","Height_total_sum","Depth_total_sum","Length_total_sum","EucDistance_avg","EucDistance_max","PathDistance_avg","PathDistance_max","Branch_Order_max","Terminal_degree_avg","Branch_pathlength_avg","Branch_pathlength_max","Contraction_avg","Partition_asymmetry_avg","Bif_ampl_local_avg","Bif_ampl_local_max","Bif_ampl_remote_avg","Bif_ampl_remote_max","Bif_tilt_local_avg","Bif_tilt_local_max","Bif_tilt_remote_avg","Bif_tilt_remote_max","Bif_torque_local_avg", "Bif_torque_local_max","Bif_torque_remote_avg","Bif_torque_remote_max","Helix_max","Fractal_Dim_avg","Fractal_Dim_max")
length(reducedMetrics)
length(pcaMetrics)
reducedMetrics
#corticalMetaData <- subset(usedMetaData,region1=="Neocortex")
#hippocampalMetaData <- subset(usedMetaData,region1=="Hippocampus")
#primateCorticalMetaData <- subset(corticalMetaData,order=="Primate")
#OtherCorticalMetaData <- subset(corticalMetaData,order!=c("Primate","Rodent"))
#rodentCorticalMetaData <- subset(corticalMetaData,order=="Rodent")

#######
# Just looking at the basic numbers and relationships between types
#######
classCounts <- function(type,typeMetaData=usedMetaData,minSize=0,...){

  counts <- aggregate(typeMetaData[,1],by=list(typeMetaData[,type]),FUN=length) 
  #get percentage of each unique value for the given metadata type
  perc <- counts$x/sum(counts$x)*100
  #add percentage column to the data.frame counts
  counts$p <- perc  
 
  
  #if count<minSize, then name the group as 'Other'
  #for(i in 1:nrow(counts)){
    #if(counts$x[i] < minSize)
      #counts$Group.1[i] <- 'Other'
  #}
  #sum aggregate rows on 'x' and 'p' and there are no duplicate 'other' values in 'Group.1'
  #counts_other <- aggregate(x = counts[,c('x','p')], by = list(counts$Group.1), FUN = sum)
  
  #arrange in descending order
  counts <- counts[order(counts$x,decreasing=T),]
  
  
  
  #return the groups that atleast minSize
  #return(subset(counts,counts$x >= minSize))
  return(counts)
}

# Choose a list of metatypes to for looking into how a given grouping breaks down along other unspecified types
usedMetaCols <- c("species","strain","protocol","gender","archive_name","cellclass1","cellclass2","region1","corticalRegion","corticalLobe","layer","age_class")
usedMetaColsDetailed <- c("species","strain","protocol","gender","archive_name","cellclass1","cellclass2","cellclass3","region1","corticalRegion","region2","region3","layer","age_class")


#print all the groups in HierarchyGrps
printGrps <- function(hierarchy){
  print(length(hierarchy))
  for(i in 1:length(hierarchy)){
    grp <- hierarchy[[i]]
    Grpttl <- paste(grp$GroupSpec,collapse=",")
    Grpcount <- grp$count
    Grpperc <- round(grp$percent,digits = 2)
    print(paste(i, Grpttl, Grpcount, Grpperc, sep = ' | '))
    #print(grp)
  }
}

#print HierarchyGrps in newick format
printNewickGrps <- function(hierarchy){
  print(length(hierarchy))
  for(i in 1:length(hierarchy)){
    grp <- hierarchy[[i]]
    Grpttl <- paste(grp$GroupSpec,collapse=",")
    Grpcount <- grp$count
    Grpttlwnum <- paste(Grpttl,Grpcount)
    Grpperc <- round(grp$percent,digits = 2)
    print(paste("(",Grpttlwnum,":",Grpperc))
  }
}

#to split the last level of main hierarchy by spltBy
Split2Grps <- function(typefrm,typSpec,spltBy,featureMetaData, typeMetaData=usedMetaData,minSize=1){
  #remove rows that have 'Others' as group name
  typefrm <- subset(typefrm,!Group.1%in%c("Others"))
    #need count of parent group to calculate the percentage
    parentSet <- getSubset2(featureMetaData, typSpec, typeMetaData)
    perc <- nrow(parentSet)
    #filter typfrm with rows that has x atleast minSize  
    #typfrmred <- subset(typefrm,typefrm$x>=minSize)
    print("in Split2Grps ")
    print(typefrm)
    #initialize the master list that will store the subgroups
    mstrlst <- list()
    #list to store groups that are merged as 'others'
    othrlst <- list()
    #append the new specification to existing typeSpec list
    nwlstspec <- list()
    specTitles <- list()
    count_o <- 0
    perc_o <- 0
    neuronames_o <- list()
    cnt <- 0
    for (i in 1:nrow(typefrm)){
      typrw <- typefrm[i,]
      #the first column has e.g., colnmvar="Neocortex"
      typval <- typrw[1,1]
      typcnt <- typrw[1,2]
      #paste(typSpec$GroupSpec,collapse="/")
      #if(typval==typSpec[length(typSpec)]){
      #  print(paste("skipping inluded columns",typval,typSpec[length(typSpec)]))
      #}else{
        #add the groups to the return list only if they are not NA and #neurons is >= 40 
        #!(is.na(typval)) && !(typval %in% c("Not reported","NA")) &&
        if(typcnt >= minSize){#} && !(typval %in% c("Not reported","NA"))){
          #print(paste("typval:",typval))
          nwlstspec <- append(typSpec, list(typval))
          #Assign name to newly added value to the specification list
          names(nwlstspec)[length(nwlstspec)] = spltBy #colnmvar
          
          #diff <- setdiff(typSpec, nwlstspec)
          #call subgroups grouped by column variable colnmvar
          x <- getSubset2(featureMetaData, nwlstspec,typeMetaData)
          #print(perc)
          #print(nrow(x))
          #Grpttl <- paste(nwlstspec$GroupSpec,collapse="/")
          #print(paste(nrow(x),Grpttl))
          #split each group by spltBy variable
          #res<- classCounts(spltBy,x,minSize)
          #print(res)
          #store grouptitle, n, % and neuron names for each group 
          #if(!(is.na(typval)) && !(typval %in% c("Not reported","NA","Others"))){  #&& nrow(x)>=minSize){
          subgrp <- list(GroupSpec = nwlstspec, count = nrow(x), percent = nrow(x)/perc*100, neuronames = x$neuron_name)
          print("adding to main list..")
          print(paste(nwlstspec,collapse="/"))
          
          mstrlst <- append(mstrlst,list(subgrp))
          }else{
            nwlstspec <- append(typSpec, list(typval))
            #Assign name to newly added value to the specification list
            names(nwlstspec)[length(nwlstspec)] = spltBy#colnmvar
            x <- getSubset2(featureMetaData, nwlstspec,typeMetaData)
            subgrp <- list(GroupSpec = nwlstspec, count = nrow(x), percent = nrow(x)/perc*100, neuronames = x$neuron_name)
            cnt <- cnt + nrow(x)
            #print(paste(paste(nwlstspec,collapse="/"),nrow(x),cnt,sep="-"))
            othrlst <- append(othrlst,list(subgrp))
          }
      #}
    }
    #group others only if more than one group is present 
    if(length(othrlst) >= 1){
      print("enter only once..")
      neuronames_o <- list()
      #add parent spec ending with "other" as GroupSpec
      nwlstspec <- append(typSpec, list('Others'))
      #Assign name to newly added value to the specification list
      names(nwlstspec)[length(nwlstspec)] = spltBy #colnmvar
      nwlstspecttl_o <- paste(nwlstspec,collapse="/")
      for(i in 1:length(othrlst)){
        x <- othrlst[[i]]
        #summate x for all rows 
        count_o <- count_o + x$count
        #summate all perc 
        perc_o <- perc_o + x$percent
        #print(x$neuronames)
        #union all neuron names
        neuronames_o <- union(neuronames_o,unlist(x$neuronames))
      }  
      subgrp_o <- list(GroupSpec = nwlstspec, count = count_o, percent = perc_o, neuronames = neuronames_o)
      #"Others" group will be only when there is an additional child branch 
      if(length(mstrlst)>=1){
        print(paste("adding other row..",nwlstspecttl_o, count_o, perc_o,length(neuronames_o)))
        mstrlst <- append(mstrlst,list(subgrp_o))
      }
    }
    #finally outside the forloop add the list to the end of the mstrlst
    return (mstrlst)
}

#groups given dataframe (datadf) into separate hierarchies (returnhgrp) based on the ordereing of level1, level2 and level3 
makeHierarchyGroups <- function(initialSpec,datadf=usedMetaData, mdatadf =usedMetaData,level1,level2,level3=null,level4=null,minSize=1){
  #list aggregates for each selection
  #move the title,#s and neurons to resultgrp
  firstsubset <- getSubset2(datadf, initialSpec, mdatadf)
  ttl <- paste(initialSpec,collapse="/")
  #add first record to the resulthgrp
  recordGrprow <- list(GroupSpec = initialSpec, count = nrow(firstsubset), percent = nrow(firstsubset)/nrow(mdatadf)*100,neuronames = firstsubset$neuron_name)
  #initialize variables before the loop
  parentsubset <- firstsubset
  tmpSpec <- initialSpec
  level1Speclst <- list()
  level2Speclst <- list()
  level3Speclst <- list()
  level4Speclst <- list()
  resulthgrp <- list()
  #add first parent to resulthgrp
  resulthgrp <- append(resulthgrp,list(recordGrprow))
  #get the first set of groups at the top level of the hierarchy
  expandParentNodeDf <- classCounts(level1,parentsubset)
  print(expandParentNodeDf)
  #loop through hierarchylevels
  #if(expandParentNodeDf[1,1]== "Not reported" && 
  if(nrow(expandParentNodeDf)==1){
      #skip adding this record to the resulthgrp
  }else{
    recordGrpdf <- Split2Grps(expandParentNodeDf, tmpSpec,level1,datadf,mdatadf,minSize)
    #eliminate Others groups from being expanded
    if(length(recordGrpdf) >= 1 && length(grep("Others",paste(tmpSpec,collapse="/")))==0){
      for(i in 1:length(recordGrpdf)){
        recordGrprow <- recordGrpdf[[i]]
        grpttl <- paste(recordGrprow$GroupSpec,collapse="/")
        #print(paste("grpttl:",grpttl))
        #save level1 derived groups
        resulthgrp <- append(resulthgrp,list(recordGrprow))
        grpttltailval <- sapply(strsplit(as.character(grpttl),"/"),tail,1)
        #add value of level1
        if(length(grep("Others",grpttltailval))==0){
          tmpSpec <- append(initialSpec, list(grpttltailval))
          #Assign name to newly added value to the specification list
          names(tmpSpec)[length(tmpSpec)] = level1
          #print(tmpSpec)
          level1Speclst <- append(level1Speclst, list(tmpSpec))
        }
      }
    }
  }
  print(paste("Added...",length(resulthgrp),length(level1Speclst)))
    
  print("*******level2********")
  #expand level2 on level1Speclst
  for(i in 1:length(level1Speclst)){
    tmpSpec <- unlist(level1Speclst[i])
    #print(tmpSpec)
    level2subset <- getSubset2(datadf, tmpSpec, mdatadf)
    expandParentNodeDf <- classCounts(level2,level2subset)
    #if(expandParentNodeDf[1,1]== "Not reported" && 
    if(nrow(expandParentNodeDf)==1){
      #skip adding this record to the resulthgrp
    }else{
      recordGrpdf <- Split2Grps(expandParentNodeDf, tmpSpec,level2,datadf,mdatadf,minSize)
      #eliminate Others groups from being expanded
      if(length(recordGrpdf) >= 1 && length(grep("Others",paste(tmpSpec,collapse="/")))==0){
        for(i in 1:length(recordGrpdf)){
          recordGrprow <- recordGrpdf[[i]]
          grpttl <- paste(recordGrprow$GroupSpec,collapse="/")
          print(paste("grpttl:",grpttl))
          #save level1 derived groups
          resulthgrp <- append(resulthgrp,list(recordGrprow))
          grpttltailval <- sapply(strsplit(as.character(grpttl),"/"),tail,1)
          #add value of level2 to level1
          if(length(grep("Others",grpttltailval))==0){
            tmpSpec_2 <- append(tmpSpec, list(grpttltailval))
            #Assign name to newly added value to the specification list
            names(tmpSpec_2)[length(tmpSpec_2)] = level2
            level2Speclst <- append(level2Speclst, list(tmpSpec_2))
          }
        }
      }
    }
  }#level1Speclst loop
  print("******level3*******")
  if(!missing(level3)){
    for(i in 1:length(level2Speclst)){
      tmpSpec <- unlist(level2Speclst[i])
      #print(tmpSpec)
      level3subset <- getSubset2(datadf, tmpSpec, mdatadf)
      expandParentNodeDf <- classCounts(level3,level3subset)
      #if(expandParentNodeDf[1,1]== "Not reported" && 
      if(nrow(expandParentNodeDf)==1){
        #skip adding this record to the resulthgrp
      }else{
        recordGrpdf <- Split2Grps(expandParentNodeDf, tmpSpec,level3,datadf,mdatadf,minSize)
        #eliminate Others groups from being expanded
        if(length(recordGrpdf) >= 1 && length(grep("Others",paste(tmpSpec,collapse="/")))==0){
          for(i in 1:length(recordGrpdf)){
            recordGrprow <- recordGrpdf[[i]]
            grpttl <- paste(recordGrprow$GroupSpec,collapse="/")
            print(paste("grpttl:",grpttl))
            #save level1 derived groups
            resulthgrp <- append(resulthgrp,list(recordGrprow))
            grpttltailval <- sapply(strsplit(as.character(grpttl),"/"),tail,1)
            #add value of level2 to level1
            if(length(grep("Others",grpttltailval))==0){
              tmpSpec_3 <- append(tmpSpec, list(grpttltailval))
              #Assign name to newly added value to the specification list
              names(tmpSpec_3)[length(tmpSpec_3)] = level3
              #print(tmpSpec_3)
              level3Speclst <- append(level3Speclst, list(tmpSpec_3))
            }
          }
        }
      }
    }#level2Speclst loop
  }
  print("******level4*******")
  
  if(!missing(level4)){
    for(i in 1:length(level3Speclst)){
      tmpSpec <- unlist(level3Speclst[i])
      print(tmpSpec)
      level4subset <- getSubset2(datadf, tmpSpec, mdatadf)
      expandParentNodeDf <- classCounts(level4,level4subset)
      #if(expandParentNodeDf[1,1]== "Not reported" && 
      if(nrow(expandParentNodeDf)==1){
        #skip adding this record to the resulthgrp
      }else{
        recordGrpdf <- Split2Grps(expandParentNodeDf, tmpSpec,level4,datadf,mdatadf,minSize)
        #eliminate Others groups from being expanded
        if(length(recordGrpdf) >= 1 && length(grep("Others",paste(tmpSpec,collapse="/")))==0){
          for(i in 1:length(recordGrpdf)){
            recordGrprow <- recordGrpdf[[i]]
            grpttl <- paste(recordGrprow$GroupSpec,collapse="/")
            print(paste("grpttl:",grpttl))
            #save level1 derived groups
            resulthgrp <- append(resulthgrp,list(recordGrprow))
            grpttltailval <- sapply(strsplit(as.character(grpttl),"/"),tail,1)
            #add value of level2 to level1
            if(length(grep("Others",grpttltailval))==0){
              tmpSpec_4 <- append(tmpSpec, list(grpttltailval))
              #Assign name to newly added value to the specification list
              names(tmpSpec_4)[length(tmpSpec_4)] = level4
              level4Speclst <- append(level4Speclst, tmpSpec_4)
            }
          }
        }
      }
    }#level2Speclst loop
  }
  print(paste("Added",length(resulthgrp),"groups"))
  return (resulthgrp)
}#end of makeHierarchyGroups

threeLevelHierarchyGrping <- function(primary, featureMetaData, typeMetaData = usedMetaData, level1="null", level2="null", level3="null", minSize=40){
  HierarchyGrps <- list()
  retGrps <- list()
  #starts with order spec
  PrimarySpec <-  primary
  print("first subset")
  x <- getSubset2(featureMetaData, PrimarySpec, typeMetaData)
  #dim(x)
  ttl <- paste(PrimarySpec,collapse="/")
  print("enter three..")
  print(ttl)
  #store specification list, n, % and neuron names for each group 
  g1 <- list(GroupSpec = PrimarySpec, count = nrow(x), percent = nrow(x)/nrow(typeMetaData)*100,neuronames = x$neuron_name)
  print(paste("Primary:",length(g1)))
  #store groups that belong to the main hierarchy
  HierarchyGrps <- append(HierarchyGrps,list(g1))
  #names(HierarchyGrps[[1]])
  #str(HierarchyGrps)
  #data frame (Group.1, x, p) group by level1
  df <- classCounts(level1,x)
  print(df)
  
  #to split the dataframe into subgroups on the given metadata type, each group is atleast minSize
  grpBy <- level1
  #subgroups by level2
  lstBy <- level2
  grp <- g1
  #print(df)
  #percGrp <- Split2Grps(df, PrimarySpec, grpBy, lstBy,featureMetaData, typeMetaData,minSize=1)
  #if(df[1,1]== "Not reported" && 
  if(nrow(df)==1){
      #print("ELIminiated groups..@level1")
      #print(df)
  }else{
     percGrp <- Split2Grps(df, PrimarySpec, grpBy, lstBy,featureMetaData, typeMetaData,minSize)
    print("class1 level******")
    #str(Primarylevel1Bylevel2Grps)
    #length(percGrp)

    if(length(percGrp) >= 1 && length(grep("Others",paste(grp$GroupSpec,collapse="/")))==0){
      for(i in 1:length(percGrp)){
        g <- percGrp[[i]]
        grpttl <- paste(g$GroupSpec,collapse="/")
        #print(grpttl)
        #get the tail part of the title
        #tailname <- sapply(strsplit(as.character(grpttl),"/"),tail,1)
        HierarchyGrps <- append(HierarchyGrps,list(g))
      }
    }
  }
  
 
  print(level2)
  if(!missing(level2)){
    #groups in the second level wiht GroupSpec size >= 3
    HierarchyGrps1 <- HierarchyGrps
    grpBy <- lstBy
    for(i in 1:length(HierarchyGrps1)){
      grp <- HierarchyGrps1[[i]]
      
      #if GroupSpec has >= three variables in the specification, then split further by given splitBy column name
      if(length(grp$GroupSpec)>=2 && length(grep("Others",paste(grp$GroupSpec,collapse="/")))==0){
      
        #print(paste(grp$GroupSpec,collapse="/"))
        x <- getSubset2(featureMetaData, grp$GroupSpec, typeMetaData)
        df <-  classCounts(grpBy,x)
        #percGrp <- Split2Grps(df, grp$GroupSpec, grpBy, lstBy, featureMetaData, typeMetaData,minSize=40)
        #if(df[1,1]== "Not reported" &&
        if(nrow(df)==1){
            #print("ELIminiated groups..@level2")
            #print(df)
        }else{
          percGrp <- Split2Grps(df, grp$GroupSpec, grpBy, lstBy, featureMetaData, typeMetaData,minSize)
          print("class2 level******")
          lstGrpSize <- length(percGrp)
          if(lstGrpSize >= 1){
            for(i in 1:length(percGrp)){
              g <- percGrp[[i]]
              #print(paste(g$GroupSpec,collapse="/"))
              HierarchyGrps1 <- append(HierarchyGrps1,list(g))
            }
          }
          #print(length(HierarchyGrps2))
        }
      }
      retGrps <- HierarchyGrps1
    }
  }
  print(level3)
  if(!missing(level3)){
    #groups in the third level with GroupSpec size >= 4 
    HierarchyGrps2 <- HierarchyGrps1
    lstBy <- level3 
    grpBy <- lstBy
    for(i in 1:length(HierarchyGrps2)){
      grp <- HierarchyGrps2[[i]]
      #if GroupSpec has >= three variables in the specification, then split further by given splitBy column name
      if(length(grp$GroupSpec)>=3 && length(grep("Others",paste(grp$GroupSpec,collapse="/")))==0){
        #print(paste(grp$GroupSpec,collapse="/"))
          print(paste(grp$GroupSpec,collapse="/"))
          x <- getSubset2(featureMetaData, grp$GroupSpec, typeMetaData)
          df <-  classCounts(grpBy,x)
          #percGrp<- Split2Grps(df, grp$GroupSpec, grpBy, lstBy, featureMetaData,typeMetaData,minSize=40)
          #if(df[1,1]== "Not reported" && 
          if(nrow(df)==1){
            
              print("ELIminiated groups..@level3")
              print(df)
          }else{
              percGrp<- Split2Grps(df, grp$GroupSpec, grpBy, lstBy, featureMetaData,typeMetaData,minSize)
              print("class3 level******")
              if(length(percGrp) >= 1){
                for(i in 1:length(percGrp)){
                  g <- percGrp[[i]]
                  print(paste(g$GroupSpec,collapse="/"))
                  HierarchyGrps2 <- append(HierarchyGrps2,list(g))
                }
              }
            }
          }
        }
    retGrps <- HierarchyGrps2
  }
  return (retGrps)
}#end of threeLevelHierarchyGrping

#Testing the columns for skewness and if skewness <0 or >1, then do logtransform and change the column name accordingly
par(ask = FALSE)
plotTesting<- function(dataframe){
  dft <- dataframe
  numcols <- ncol(dft)
  ctr <- 0
  #skip metadata columns
  for(i in 33:numcols){
    if(is.numeric(dft[,i])){
      mycol <- dft[,i]
      mycol <- subset(mycol, !is.na(mycol))
      oldname <- names(dft)[i]
      newname <- paste("ln(",oldname,")",sep="")
      hist(mycol,breaks=30,xlab=oldname,main=paste("Histogram of ",oldname," before transformation"))
      #condition to be satisfied to perfrom log transformation
      if(skewness(mycol)< -0.8 || skewness(mycol)> 0.8){
        ctr <- ctr+1 
        print(paste("performing ln transformation ",newname))
        t <- logTransformData(mycol)
        hist(t, breaks=30, xlab=paste("ln(",newname,")"),main=paste("Histogram of ",newname," after transformation"))
        #replace back the non-NA values in the dataframe
        dft[,i][!is.na(dft[,i])]<- t
        #print(dft[,1][!is.na(dft[,i])])
        #setdiff(dft[,i])
        #dft[,i] <- t
        names(dft)[i] = newname
      }else{
        print(paste("skipping ln transformation ",oldname))
      } 
    }
  }
  print("transformed:"+ctr)
  return(dft)
}
#compute total scatter distance from a given two column data frame and centroid point
scatter <- function(p,xc,yc){
  #compute and return total sum of scatter from the centroid point
  #print(paste(names(p)[1],names(p)[2]))
  ss <- 0   
  for(i in 1:nrow(p)){
    #first column will be considered as feature 1 on X-axis
    x <- p[i,1]
    #second column will be considered as feature 2 on Y-axis
    y <- p[i,2]
    #distance between the point and the centroid
    #print(paste(x,y))
    #print(paste(xc,yc))
    s <- sqrt((x-xc)^2 + (y-yc)^2)
    ss <- ss + s
  }
  #print(ss)
  return (ss)    
}

#compute scatter distance b/w group pairs
grppairscatterDistance <- function(p,m1colname,m2colname){
  #print(paste(m1colname,m2colname))
  #create a dataframe to return with a dummy row for group pairs, and their respective scatter, distance and separation ratios for a given pair of features.
  grppairdf <- data.frame("f1"= 'f1',"f2"='f2',"G1scatter"=0.0,"G2scatter"=0.0,"G1G2Dist"=0.0,"avgsep"=0.0)
  grppairdf <- grppairdf[-1,]
  #makesure there are no NAs in both columns
  xcol <- p[,m1colname]
  ycol <- p[,m2colname]
  xpos <- is.na(xcol)
  ypos <- is.na(ycol)
  print(paste(length(xpos[xpos==TRUE]),length(ypos[ypos==TRUE])))
  #checking and subsetting both columns to same size after removing NA values
  if(length(xpos[xpos==TRUE])>0){
    pdf <- subset(p,!xpos)
    ycol <- pdf[,m2colname]
    ypos <- is.na(ycol)
    if(length(ypos[ypos==TRUE])>0)
      pdf <- subset(pdf,!ypos)
    print(nrow(pdf))
  }else if(length(ypos[ypos==TRUE])>0){
    pdf <- subset(p,!ypos)
    xcol <- pdf[,m1colname]
    xpos <- is.na(xcol)
    if(length(xpos[xpos==TRUE])>0)
      pdf <- subset(pdf,!xpos)
    print(nrow(pdf))
  }else{
    #if neither columns have NA values make sure the groupnum has no NA values
    pdf <- subset(p,!is.na(p$groupnum))
    print(nrow(pdf))
  }
  
  #pdf <- subset(p, !is.na(p$groupnum))
  #get unique groupnum values
  grp <- unique(pdf$groupnum)
  #print(p[,c(m1colname,m2colname)])
  #aggregate for feature pair m1 and m2, each row has centroid point for each group
  pmeandf <- aggregate(pdf[,c(m1colname,m2colname)],by=list("Group" = pdf$groupnum),mean)
 
  #group 1 for a given feature pair
  for(i in 1:(nrow(pmeandf)-1)){
    pRawdf <- subset(pdf[,c(m1colname,m2colname)],pdf$groupnum==pmeandf$Group[i])
    #print(paste("pRawdf",nrow(pRawdf)))
    #xc <- pmeandf$m1colname[i]
    xc <- pmeandf[i,2]
    #yc <- pmeandf$m2colname[i]
    yc <- pmeandf[i,3]
 
    #compute scatter for group 1
    g1scatter <- scatter(pRawdf,xc,yc)
    #print(paste("g1scatter:",g1scatter))
    #group 2 for a given feature pair
    for(i2 in (i+1):nrow(pmeandf)){
      pRawdf2 <- subset(pdf[,c(m1colname,m2colname)],pdf$groupnum==pmeandf$Group[i2])
      #print(paste("pRawdf2",nrow(pRawdf2)))
      #xc2 <- pmeandf$m1colname[i2]
      xc2 <- pmeandf[i2,2]
      #yc2 <- pmeandf$m2colname[i2]
      yc2 <- pmeandf[i2,3]
      #print(paste("centroid2:",xc2,yc2))
      #compute scatter for group 2
      g2scatter <- scatter(pRawdf,xc2,xc2)
      #print(paste("g2scatter:",g2scatter))
      g1g2distance <- sqrt((xc-xc2)^2 + (yc-yc2)^2)
      if(!(g1scatter == 0) && !(g2scatter == 0)){
        #create a new row to be added to the dataframe
        newrow <- data.frame("f1"=m1colname, "f2"=m2colname, "G1scatter"=g1scatter,"G2scatter"=g2scatter,"G1G2Dist"=g1g2distance,"avgsep"=g1g2distance/((g1scatter+g2scatter)/2))
        
      }
      grppairdf <- rbind(grppairdf, newrow)
    }
  }
  #average scatter distances over all groups for a given feature pair
  #both mean and median of average separate distance is added to each pair
  print(paste("global average of",nrow(grppairdf)))
  avgdistanceperpair <- aggregate(grppairdf$avgsep,by=list("f1"=grppairdf$f1,"f2"=grppairdf$f2),function(x) cbind(mean(x),median(x)))
  #print(paste("added",nrow(grppairdf)))
  #print(grppairdf)
  #return (grppairdf)
  return(avgdistanceperpair)
}


#plotting the crosshair in scatter plot using ggplot2
library(ggplot2)
#plots crosshair of mean and SD on xy axes. Also computes if pair of groups are separated without overlapping.
crosshairscatter<- function(pdataframe,m1,m2,hierarchyname){
  #only for brain regions convert the width to normal scale
  pdataframe[,m1] <- exp(pdataframe[,m1])
  #makesure there are no NAs in both columns
  xcol <- pdataframe[,m1]
  ycol <- pdataframe[,m2]
  print(ycol[order(ycol)])
  xpos <- is.na(xcol)
  ypos <- is.na(ycol)
  print(paste(length(xpos[xpos==TRUE]),length(ypos[ypos==TRUE])))
  #checking and subsetting both columns to same size after removing NA values
  if(length(xpos[xpos==TRUE])>0){
    pdataframe <- subset(pdataframe,!xpos)
    ycol <- pdataframe[,m2]
    ypos <- is.na(ycol)
    if(length(ypos[ypos==TRUE])>0)
      pdataframe <- subset(pdataframe,!ypos)
    print(nrow(pdataframe))
  }else if(length(ypos[ypos==TRUE])>0){
    pdataframe <- subset(pdataframe,!ypos)
    xcol <- pdataframe[,m1]
    xpos <- is.na(xcol)
    if(length(xpos[xpos==TRUE])>0)
      pdataframe <- subset(pdataframe,!xpos)
    print(nrow(pdataframe))
  }else{
    #if neither columns have NA values make sure the groupnum has no NA values
    pdataframe <- subset(pdataframe,!is.na(pdataframe$groupnum))
    print(nrow(pdataframe))
  }
  
  #selects a subset where groupttl is != 'Others'
  #p <- subset(pdataframe,!sapply(strsplit(as.character(pdataframe$groupttl),"Control/"),tail,1)%in%c("/Others"))
  p <- subset(pdataframe,!grepl("Other", sapply(strsplit(as.character(pdataframe$groupttl),"Control/"),tail,1)))
  group <- sapply(strsplit(as.character(p$groupttl),"Control/"),tail,1)
  group1 <- group
  print(dim(p))
  print(unique(group))
  print(paste(length(p[,m1]),length(p[,m2])))
  print(paste(m1,"vs",m2))
  ctr <- 0
  
  #xlabel = column1,
  #ylabel = column2,
  #first compute the summary metrics N, mean, and SD on X and y axes
  z1 <- ddply(p,.(groupnum,group),function(dfp,column1,column2) {
    c(N = nrow(dfp),
       x = median(dfp[,column1]), #mean(dfp[,column1]),
       y = median(dfp[,column2]), #mean(dfp[,column2]),
       xQ1 = quantile(dfp[,column1], 0.25,names=FALSE)[1],
       xQ3 = quantile(dfp[,column1], 0.75,names=FALSE)[1],
       yQ1 = quantile(dfp[,column2], 0.25,names=FALSE)[1],
       yQ3 = quantile(dfp[,column2], 0.75,names=FALSE)[1],      
       #xSE = sqrt(var(dfp[,column1]))/nrow(dfp),
       #ySE = sqrt(var(dfp[,column2]))/nrow(dfp),
       #xSD = sd(dfp[,column1]),
       #ySD = sd(dfp[,column2]),
       #xmin = log(mean(dfp[,column1]) - sd(dfp[,column1])),
       #xmax = log(mean(dfp[,column1]) + sd(dfp[,column1])),     
       #ymin = log(mean(dfp[,column2]) - sd(dfp[,column2])),
       #ymax = log(mean(dfp[,column2]) + sd(dfp[,column2])),       
       othr = 1)}, m1, m2)
  print(z1)
  
  #ctr <- crosshairSeparation(z,group)
  
  #geom_text(aes(label = paste("ln(",m1,")"))) +
  #geom_text(aes(label = paste("ln(",m2,")"))) +
  #xlab(paste("ln(",m1,")")) +
  #ylab(paste("ln(",m2,")")) +
  #,xlab=xlabel,ylab=ylabel
  #plot the crosshairs
 
  #selects a subset of p, where groupttl is 'Others
  #p <- subset(pdataframe,sapply(strsplit(as.character(pdataframe$groupttl),"Control/"),tail,1)%in%c("Others"))
  p <- subset(pdataframe,grepl("Other", sapply(strsplit(as.character(pdataframe$groupttl),"Control/"),tail,1)))
  #get the parent description for "Others" by adding "(O) at the end
  #t <- sapply(strsplit(as.character(p$groupttl),"Others"),head,1)
  t <- as.character(p$groupttl)
  #group <- paste(t,"(O)",sep="")
  group <- t
  group2 <- group
  print(unique(group2))
  #to minimize the length, choose the tail part of the parent description
  group <- sapply(strsplit(as.character(group2),"Control/"),tail,1)
  #group[group =="Control"] <- "Others"
  print(unique(group))
  #group <- sapply
  print(dim(p))
  #compute summary of 'Others' groups
  z2 <- ddply(p,.(groupnum,group),function(dfp,column1,column2) {
    c(N = nrow(dfp),
      x = median(dfp[,column1]),
      y = median(dfp[,column2]),
      xQ1 = quantile(dfp[,column1], 0.25, names=FALSE)[1], # first quartile
      xQ3 = quantile(dfp[,column1], 0.75, names=FALSE)[1], # third quartile
      yQ1 = quantile(dfp[,column2], 0.25, names=FALSE)[1],
      yQ3 = quantile(dfp[,column2], 0.75, names=FALSE)[1],      
      #xSE = sqrt(var(dfp[,column1]))/nrow(dfp),
      #ySE = sqrt(var(dfp[,column2]))/nrow(dfp),
      #xSD = sd(dfp[,column1]),
      #ySD = sd(dfp[,column2]),
      #xmin = log(mean(dfp[,column1]) - sd(dfp[,column1])),
      #xmax = log(mean(dfp[,column1]) + sd(dfp[,column1])),     
      #ymin = log(mean(dfp[,column2]) - sd(dfp[,column2])),
      #ymax = log(mean(dfp[,column2]) + sd(dfp[,column2])),     
      othr = 2)}, m1, m2)
  print(z2)
  
  #group <- list()
  #group <- append(group1,group2)
  #print(unique(group))
  #group <- append(group, list(group2))
  #print(unique(group))
  
  #a1 = annotate("text", x = zz$x, y = zz$y, label = unique(group2))
  z <- rbind(z1,z2)
  #groupttl <- z$group #reassigning group argument with Z
  groupnum <- z$groupnum #reassigning the group num with z
  print(unique(groupnum))
  othr <- z$othr
  print(paste(nrow(z),nrow(z2)))
  #zp <- z
  #zp$group <- c("a1","a2","a3","a4","a5","a6","a7","a8","a9","b1")
  #write.table(zp, file = "zp.txt",sep = ",", col.names=TRUE)
  print(length(unique(z$group)))
  print(unique(z$group))
  z[,"group"] <- factor(z[,"group"], levels = c("Salamander","Goldfish", "Other Bony fish","Cat","Blowfly","Human","Monkey","Elephant","Mouse","Rat","Other Rodents","Other Species"))
  #z[,"group"] <- factor(z[,"group"], levels = c("Hippocampus/CA1","Other Hippocampus","Neocortex/Frontal/Motor","Neocortex/Frontal/Prefrontal","Other Neocortex/Frontal","Neocortex/Insula","Neocortex/Occipital","Neocortex/Parietal/S1/L4","Other Neocortex/Parietal/S1","Other Neocortex/Parietal","Other Neocortex","Olfactory bulb","Retina","Other Brain regions"))
  #z[,"group"] <- factor(z[,"group"], levels = c("Basket","Other Interneurons","Ganglion/Monostratified","Other Ganglions","Granule","Medium spiny","Motoneuron","Magnopyramidal","Other Pyramidals","Other Principals"))
  print(z)
  #require(reshape2)
  #mytext <- melt(z, id.vars = "group")
  #print(dim(mytext))
  #print(colnames(mytext))
  #print(mytext)
  #axpos <- ln(seq(100,400,25))
  #axpos <- ln(c(100,125,150,200,250,300,400,500))
  axpos <- ln(c(5,10,20,40,80,160,320,640))
  #axpos <- ln(c(300, 1000, 3000, 10000, 30000, 100000))
  
  #axval <- seq(100,400,25)
  #axval <- c(100,125,150,200,250,300,400,500)
  #axval <- c(300, 1000, 3000, 10000, 30000, 100000)
  axval <- c(5,10,20,40,80,160,320,640)
  print(axval)
  length(axval)
  round2 = 0
  chp <- ggplot(data=z,aes(x=x,y=y,colour = group,cex.lab=1.5,linetype=factor(othr))) + 
    #scale_x_continuous(breaks=axpos,labels=round(axval,round2)) +
    #scale_y_continuous(breaks=axpos,labels=axval) +
    #xlim=c(floor(min(x)),ceiling(max(x))), ylim=c(floor(min(y)),ceiling(max(y))),
    #opts(aspect.ratio = 2/(1+sqrt(5)) ) +
    #geom_point(aes(colour = factor(groupnum))) +
    xlab("Tree Z span (um)") + #"Tree arbor length"
    ylab("Topological Asymmetry") + #"Tree arbor fractal"
    #xlab("Average Bifurcation angle ( )") + #"avg bif angle"
    #ylab("Tree Width (um)") + #"total width"
    #xlab("Total arbor length (um)") + #"total length"
    #ylab("Fractal dimension") + #"average fractal"
    #xlab(m1) +
    #ylab(m2) +
    geom_point(shape=19,size=1.5)+
    geom_errorbar(aes(ymin = yQ1, ymax = yQ3), height=0, width = 0) + # colour = factor(groupnum), linetype = factor(othr))) + 
    geom_errorbarh(aes(xmin = xQ1, xmax = xQ3), height = 0, width = 0) + #, colour = factor(groupnum), linetype = factor(othr))) +
    
    #geom_text(data = mytext, aes(x = group, y = value , label = value), size=4) +
    guides(linetype=FALSE) + #remove the legend added by linetype
    theme_bw(10) + #remove background and increase the size
    #scale_colour_manual(name=hierarchyname, values=factor(unique(groupnum))) + #, breaks=factor(unique(group))) +
    scale_colour_hue(name=hierarchyname, # Legend label, use darker colors
                     l=50,       # Use darker colors, lightness=50
                     c=100) +   #chroma (intensity of color)
    #change line type in the legend manually
    
    #guides(colour=guide_legend(override.aes=list(lwd=4)))
    guides(colour=guide_legend(override.aes=list(linetype=c(rep(1,times=2),rep(2,times=1),rep(1,times=7),rep(2,times=2)),lwd=2)))
    #guides(colour=guide_legend(override.aes=list(linetype=c(rep(1,times=1),rep(2,times=1),rep(1,times=2),rep(2,times=1),rep(1,times=3),rep(2,times=3),rep(1,times=2),rep(2,times=1)),lwd=2)))
    #guides(colour=guide_legend(override.aes=list(linetype=c(rep(1,times=1),rep(2,times=1),rep(1,times=1),rep(2,times=1),rep(1,times=4),rep(2,times=2))),lwd=2))
    #c("Uniglomerular Projecting","Other Uniglomerular Projections","Other Axonal terminals","Basket","Martinotti","Other Interneurons","Ganglion/Monostratified","Other Ganglions","Granule","Medium spiny","Motoneuron","Magnopyramidal","Other Pyramidals","Other Principals")
    #labels=c("Salamander","Goldfish","Cat","Blowfly","Cricket","Drosophila","Human","Monkey","Elephant","Mouse","Rat","Other Rodents")
    #guides(colour=guide_legend(override.aes=list(linetype=c(rep(2,times=1),rep(1,times=1),rep(2,times=2),rep(1,times=2),rep(2,times=2),rep(1,times=4),rep(2,times=1),rep(1,times=1))),lwd=2)) 
                     #breaks=unique(group),
                     #labels=unique(group),#add the legend for all groups
                 
    #scale_linetype_manual(values = c(1,2))
    #geom_errorbar(data=z[(nrow(z)-nrow(z2)):nrow(z2),], aes(x=x, y=y, ymin = y - ySD, ymax = y + ySD, colour = group),linetype="dashed") + 
    #geom_errorbarh(data=z[(nrow(z)-nrow(z2)):nrow(z2),], aes(x=x, y=y, xmin = x - xSD, xmax = x + xSD, colour = group), linetype="dashed")
    
    #scale_colour_manual(values=unique(group))
    #geom_errorbar(aes(ymin = y - ySE, ymax = y + ySE, colour = group)) + 
    #geom_errorbarh(aes(xmin = x - xSE, xmax = x + xSE, colour = group)) +
    
    #ggplotobj <- chp + a1
  
  #plot(grid_chp)
  pname <- paste(m1,"Vs",m2,"_",nrow(z),hierarchyname,sep="")
  print(pname)
  ggsave(filename=paste(pname,".tiff",sep=""),
         plot=chp +theme(aspect.ratio=1))
  chp
}

library(gridExtra)
g_legend<-function(aggplotobj){
  tmp <- ggplot_gtable(ggplot_build(aggplotobj))
  print(tmp$grobs$name)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  print(leg)
  legend <- tmp$grobs[[leg]]
  legend
}

#plotting the scatter plots amongst metric pairs
scatterPlotTesting<- function(plotDf,m1,m2){
  #par(xpd=FALSE)
  hist(plotDf[,m1],breaks=30,xlab=m1,main=paste("Histogram of ",m1," before transformation"))
  x <- logTransformData(plotDf[,m1])
  hist(x, breaks=30, xlab=paste("log(",m1,")"),main=paste("Histogram of ",m1," after transformation"))
  #print(paste("x:",length(x)))
  #hist(x, width=0.33, offset=0.00, col="blue", main=paste("Histogram of", m1))
  
  hist(plotDf[,m2],breaks=30,xlab=m2,main=paste("Histogram of ",m2," before transformation"))
  y <- logTransformData(plotDf[,m2])
  hist(y, breaks=30, xlab=paste("log(",m2,")"),main=paste("Histogram of ",m2," after transformation"))
 
  print(paste("x:",length(x),"y:",length(y)))
 
  #with(plotDf, which(plotDf$groupnum == Grpnum)
  
  # so turn off clipping:
  #par(xpd=TRUE)
  #legend(2.8,-1,c("group A", "group B"), pch = c(1,2), lty = c(1,2))
  #hist(y, width=0.33, offset=0.00, col="red", main="Histogram of HillmanThreshold max")
  #plot(y, x, pch = 3, col = c(unique(plotDf[!is.na(plotDf$groupnum),"groupnum"])),main=paste(m1," vs ",m2),xlab=m1, ylab=m2)
  #reg1 <- lm(x~y)
  #plot(x, y,col = c(unique(plotDf[!is.na(plotDf$groupnum),"groupnum"])),main=paste(m1," vs ",m2),xlab=m1, ylab=m2)
  plot(x, y,col = p$groupnum, main=paste(m1," vs ",m2), cex = .5, xlab=m1, ylab=m2)
  #legend("topright",legend = sapply(strsplit(as.character(unique(p$groupttl)),"/"),tail,1), col = c(13,14,10,19,5,9,15,4,7,8,3,12,2,16,6,1,18,11), pch = 1, cex=.7)
  #legend("topright",legend = c("sprague-dawley","wistar","C57BL6","human","macaque","trangenic mice","cat","rhesus","CD-1","BC57 black","melanogaster","long-evans","african elephant","goldfish","sjl","salamander","blowfly","fischer 344"), col = c(13,14,10,19,5,9,15,4,7,8,3,12,2,16,6,1,18,11), pch = 1, cex=.7)
  #legend("topright",legend = c("cat","salamander","goldfish","fly","human","elephant","monkey","mouse","rat"), col = c(16,1,15,3,17,2,4,10,14), pch = 1, cex=.7)
  #abline(reg1)
  #legend("bottomright",legend = sapply(strsplit(mydataframe$neuron_name, "\\."), "[[", length), cex=.7,col = c(unique(plotDf[!is.na(plotDf$groupnum),"groupnum"])))
  #text(x, y, plotDf$groupnum, cex=0.7, pos=4, col="red")
  
  #m3 <- mean((x-mean(x))^3)
  #skew <- m3/(sd(x)^3)
}

#add the given Grpnum and grpttl to the data frame that has neuron_name column
add2GrpCol <- function(mydata,grpttl,neuronList,Grpnum){
print(paste("# neurons..",length(neuronList)))

 for(i in 1:length(neuronList)){
   #posids <- with(mydata,grepl("neuron_name",neuronList[i]))
   #mydata[posids,"groupnum"] <- Grpnum
   #print(paste(mydata$neuron_name,",",neuronList[i]))
   pos <- mydata$neuron_name == neuronList[i]
   #pos <-  with(mydata,which(mydata$neuron_name==neuronList[i]))
   if(length(pos[pos==TRUE])==1){
    mydata[pos==TRUE,"groupnum"] <- Grpnum
    mydata[pos==TRUE,"groupttl"] <- grpttl
   }else{
     print(neuronList[i])
   }
   
 }
 #print(mydata$groupnum[mydata$groupnum == Grpnum])
 print(paste("compare..",nrow(mydata),length(with(mydata, which(mydata$groupnum == Grpnum)))))
 #print(paste("Adding..",Grpnum))
 #print(length(pos[pos==TRUE]))
 return (mydata)
}

#filtering the groups by keeping the most specific groups and removing the parent group
addgrps2MetricDF <- function(mydata=usedMetaData,specifications=NULL){
  grpcol <- mydata$neuron_name
  tot_len <- 0
  print(length(grpcol))
  if (!is.null(specifications)){
    #print(str(specifications))
    print(length(specifications))
    GrpnameList <- list()
    FilteredSpecList <- list()
    neuronList <- list()
    filteredGrpList <- list()
    #tmpGrpList <- list()
    #start parsing the specifications list backwards, as the most specific group is the last one.
    for (i in length(specifications):1){
      #specType <- names(specifications)[i]
      #print(specType)
      if (is.list(specifications)){
        listelement <- specifications[[i]]
        #listelement[1] is group specification list. 
        GrpSpec <- listelement[1]$GroupSpec
        Grpname <- paste(GrpSpec,collapse="/")
        #print(paste("New filtered list has",length(FilteredSpecList)))
        if(length(FilteredSpecList)==0){
          print("adding first Group to the list ")
          FilteredSpecList <- append(FilteredSpecList,list(Grpname))
          #listelement[4] is list of neuron_names
          neuronList <- unlist(listelement[4])
          #mydata <- add2GrpCol(mydata,neuronList,length(FilteredSpecList))
          #print(unique(mydata$groupnum))
          tot_len <- length(neuronList)
          tmpList <- list(GroupSpec = Grpname, neuronList = list(neuronList))
          filteredGrpList <- append(filteredGrpList, list(tmpList))
        }else if(length(grep(Grpname,FilteredSpecList))==0){#if Grpname pattern is not found, then add that group to the FilteredSpecList
          FilteredSpecList <- append(unique(FilteredSpecList),list(Grpname))
          neuronList <- unlist(listelement[4])
          print(paste("adding..",Grpname,length(neuronList)))
          #print(neuronList)
          #tmpGrpList <- append(unique(tmpGrpList), unique(neuronList))
          #mydata <- add2GrpCol(mydata,neuronList,length(FilteredSpecList))
          #print(unique(mydata$groupnum))
          tot_len <- tot_len + length(neuronList)
          tmpList <- list(GroupSpec = Grpname, neuronList = list(neuronList))
          filteredGrpList <- append(filteredGrpList, list(tmpList))
        }
        print(paste("neuronList:",Grpname, tot_len))
      }
    }
  }
  print(paste("tot_len",tot_len))
  print(FilteredSpecList)
  return (filteredGrpList)
}

#add groupnum and groupttl according to Grplst to mydf
metricDf_plotting <- function(mydf, Grplst){
  for(i in 1:length(Grplst)){
    ttl <- Grplst[[i]]$GroupSpec
    neuronkey <- unlist(Grplst[[i]]$neuronList)
    print(paste(ttl,"-",length(neuronkey)))
    mydf <- add2GrpCol(mydf,ttl,neuronkey,i)
    print(unique(mydf$groupnum))
  }
  return (mydf)
}

#generates a subset dataframe with metrics and metadata to plot the XY crosshairs
getPlotSubset <- function(metricDf, metrics = reducedMetrics, metadata = usedMetaData){
  #total columns minus neuroname
  maxNAcols <- length(metricDf)-1
  #remove the .CNG.swc extension if present from the neuron_name column
  #metricDf <- clnData(metricDf)
  #include columns that are specified in metrics
  metricDf <- metricDf[, metrics]
  #check for NA values
  metricDf <- cntNRmvNAs(metricDf,maxNAcols)
  #check for negative values
  negMetrics <- chkNegVals(metricDf)
  if(length(negMetrics)>0){
    print("negative columns")
    print(negMetrics[negMetrics==TRUE])
  }else{
    print("No negative metrics")
  }
  print("merging...")
  #merge the edited metaData and metricDf
  metricmeta <- merge(x=metadata,y=metricDf, by = intersect(x$neuron_name, y$neuron_name), by.x="neuron_name",by.y="neuron_name")
  print(dim(metadata))
  print(dim(metricDf))
  print(dim(metricmeta))
  return(metricmeta)
}

pcaAnalysis <- function(df){
  #use reduced PCA metrics
  lmprcomp <- prcomp(LMreduced, center = TRUE, scale = TRUE)
  ls(lmprcomp)
  summary(lmprcomp)
  #choose eigenvalues >= 1
  lmprcomp$sdev ^ 2
  # pick upto 25 components as the change in variance reaches is pretty close to 0.
  scree(lmprcomp, npcs = 97, main = "LM metrics",xlab="Components")
  scree(lmprcomp, npcs = 97, type = "line", main = "LM metrics")
}
library(ggplot2)
library(plyr)

####Generate scatter plots from data######
#arbor type == all

aggregate(metricDendrites$expercond, by = list(metricDendrites$expercond), length)
primary <- list(expercond = "Control")
level1 <- "order"
level2 <- "species"
level3 <- "strain"
spstrainMetricDend <- makeHierarchyGroups(primary, metricDendrites,usedMetaData,level1,level2,minSize=55)

####creating the cv table for metrics used in the final analysis####

#CVs for the neurons from the hierarchies - part I
hierarchyneuronames <- spstrainMetricDend[[1]]$neuronames
colnames(metricDendrites)
hierarchyMetricdf <- subset(metricDendrites, metricDendrites$neuron_name%in%hierarchyneuronames)[33:63]
dim(hierarchyMetricdf)
chkMeanCols(hierarchyMetricdf)
hierarchygrpscvarr <- cvfun(hierarchyMetricdf)/100
#CVs for the neurons from the cluster analysis - part II
dim(Dendritedatamatrix)
pcaclusterMetricdf <- subset(metricDendrites, metricDendrites$neuron_name%in%Dendritedatamatrix$neuron_name)[33:63]
dim(pcaclusterMetricdf)
chkMeanCols(pcaclusterMetricdf)
pcaclustergrpscvarr <- cvfun(pcaclusterMetricdf)/100
#table 1 with two columns with CVs from part I (N=7143) & II(N=5099)
cvframe <- data.frame("metrics"=colnames(hierarchyMetricdf), "CV for part I" = hierarchygrpscvarr, "CV for part II"=pcaclustergrpscvarr)
cvframe <- cvframe[order(cvframe$CV.for.part.I,cvframe$CV.for.part.II,decreasing=TRUE),]
cvframe
write.table(cvframe, file = "Table1_cv.txt",sep = ",", row.names=FALSE)


#spMetricall <- threeLevelHierarchyGrping(primary, metricwhole,typeMetaData=usedMetaData, level1,level2,level3,minSize=0)
#printGrps(spMetricall)
printGrps(spstrainMetricDend)
#log transform necessary columns
plotDend_sp <- plotTesting(metricDendrites)
dim(usedMetaData)
dim(metricDendrites)
dim(plotDend_sp)
plotDend_sp$'ln(Depth_total_sum)'
printGrps(spstrainMetricDend)
Grplst <- addgrps2MetricDF(plotDend_sp,spstrainMetricDend)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}
aggregate(metricDendrites$expercond, by = list(metricDendrites$expercond), length)
plotDend_sp <- metricDf_plotting(plotDend_sp,Grplst)
dim(plotDend_sp)
colnames(plotDend_sp)
unique(plotDend_sp$groupttl)
unique(plotDend_sp$groupnum)
#remove Control/Others from the groupttl
#pos <- grepl("Control/Others", plotDend_sp$groupttl)
#length(pos[pos==TRUE])
#plotDend_sp$groupnum[pos==TRUE] <- NA
NApos <- is.na(plotDend_sp$groupnum)
length(NApos[NApos==TRUE])
length(NApos[NApos==FALSE])
psp <- subset(plotDend_sp,!NApos)
colnames(psp)
dim(psp)
#psp <- plotDf_sp[,reducedMetrics]
colnames(psp)
dim(psp)
psp$'ln(Depth_total_sum)'
psp$'ln(Branch_Order_max)'

unique(psp$groupttl)
corrdend_sp <- leastcorrelated(psp[,33:63],0.80)
#adjust pval with bonferronicorrection
corrdend_sp$pval <- adjustpval(corrdend_sp$pval,nrow(corrdend_sp))
colnames(corrdend_sp)
corrdend_sp$pval
print(nrow(corrdend_sp))
str(corrdend_sp)
#corrdf_sorted <- corrdf[with(corrdf, order(f1, f2)), ] 
#corrdf_sorted

#normalize metric dataframe before computing the master scatter separation matrix
for(i in 33:63){
  #get non-NA values for normalization
  colv <- psp[,i][!is.na(psp[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  psp[,i][!is.na(psp[,i])] <- colv/msd
  print(paste(names(psp)[i],msd))
}
psp_scatter <- psp[,colnames(psp)%in%reducedlnMetrics]
dim(psp_scatter)
colnames(psp_scatter)

colnames(psp_scatter)
dim(psp)
#make the ln(Pk_classic_avg) more meaningful by changing to ln(RallsRatio_avg)
#names(p)[49] <- "ln(RallsRatio_avg)"
#p$'ln(RallsRatio_avg)'
#names(p)[49] <- "ln(Pk_classic_avg)"
unique(psp_scatter$groupttl)
unique(psp_scatter[,c("groupttl","groupnum")])

psp_scatter[psp_scatter$groupttl %in% c("Control/Others","Control/Amphibians/Others") ,"groupttl"] <-  "Control/Other Species"
psp_scatter[psp_scatter$groupnum == 1, "groupnum"] <- 10
#combine the two frog spinal cord neurons with other species
#psp_scatter[psp_scatter$groupttl== "Control/Amphibians/Others" ,"groupttl"] <-  "Control/Other Species"
psp_scatter[psp_scatter$groupttl== "Control/Bony fishes/Others" ,"groupttl"] <-  "Control/Other Bony fish"
psp_scatter[psp_scatter$groupttl== "Control/Rodent/Rat" ,"groupttl"] <-  "Control/Rat"
psp_scatter[psp_scatter$groupttl== "Control/Rodent/Mouse"  ,"groupttl"] <-  "Control/Mouse" 
psp_scatter[psp_scatter$groupttl== "Control/Rodent/Others"  ,"groupttl"] <-  "Control/Other Rodents"
psp_scatter[psp_scatter$groupttl== "Control/Primate/Human","groupttl"] <-  "Control/Human"
psp_scatter[psp_scatter$groupttl== "Control/Primate/Monkey","groupttl"] <-  "Control/Monkey"
psp_scatter[psp_scatter$groupttl== "Control/Carnivora","groupttl"] <-  "Control/Cat"
#psp_scatter[psp_scatter$groupttl== "Control/Insects/Drosophila","groupttl"] <-  "Control/Drosophila"
psp_scatter[psp_scatter$groupttl== "Control/Bony fishes/Goldfish","groupttl"] <-  "Control/Goldfish"
psp_scatter[psp_scatter$groupttl== "Control/Proboscidae","groupttl"] <-  "Control/Elephant"
psp_scatter[psp_scatter$groupttl== "Control/Amphibians/Salamander","groupttl"] <-  "Control/Salamander"
psp_scatter[psp_scatter$groupttl== "Control/Insects","groupttl"] <-  "Control/Blowfly"
#psp_scatter[psp_scatter$groupttl== "Control/Insects/Cricket","groupttl"] <-  "Control/Cricket"


nrow(sp_corr_pval_scattertable)
colnames(sp_corr_pval_scattertable)
  for(i in 1:nrow(sp_corr_pval_scattertable)){
    if(sp_corr_pval_scattertable$r2[i]<0.03 && sp_corr_pval_scattertable$pval[i]<0.05 && sp_corr_pval_scattertable$avgsep[i,1]>0.01){
      print(paste(sp_corr_pval_scattertable$f1[i],sp_corr_pval_scattertable$f2[i]))
      if(length(psp_scatter[,as.character(sp_corr_pval_scattertable$f1[i])])>0 && length(psp_scatter[,as.character(sp_corr_pval_scattertable$f2[i])])>0)
        crosshairscatter(psp_scatter,as.character(sp_corr_pval_scattertable$f1[i]),as.character(sp_corr_pval_scattertable$f2[i]),'species')
    }
  }

x11()
par(mfrow=c(3,1))

crosshairscatter(psp_scatter,'ln(Depth_total_sum)','Partition_asymmetry_avg','Species')


#plot for brain regions also
dim(metricDendrites)
aggregate(metricDendrites$expercond, by = list(metricDendrites$lobes), length)
#brainregion hierarchy

#primary <- list(expercond = "Control",region1="Neocortex")
primary <- list(expercond = "Control")

level1 <- "region1"
level2 <- "lobes"
level3 <- "region2"
level4 <- "region3"
#unique(usedMetaData$lobes)
unique(usedMetaData$lobes)
brMetricHierarchyDend <- makeHierarchyGroups(primary, metricDendrites, usedMetaData, level1, level2,level3,level4,minSize=300)
#brMetricHierarchy <- makeHierarchyGroups(primary, metricDendrites, usedMetaData, level1, level2,level3,level4,minSize=300)

#brMetricHierarchy
printGrps(brMetricHierarchyDend)
#str(brMetricHierarchy)

Grplst <- addgrps2MetricDF(metricDendrites,brMetricHierarchyDend)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}

plotDf_br <- plotTesting(metricDendrites)

dim(plotDf_br)
plotDf_br <- metricDf_plotting(plotDf_br,Grplst)
dim(plotDf_br)
colnames(plotDf_br)
plotDf_br$'Bif_ampl_local_avg'
NApos <- is.na(plotDf_br$groupnum)
length(NApos[NApos==TRUE])
length(NApos[NApos==FALSE])
pbr <- subset(plotDf_br,!NApos)
colnames(pbr)
dim(pbr)
pbr$'ln(Soma_Surface_total_sum)'
corrdf_br <- leastcorrelated(pbr[,33:63],0.99)
#adjust pval with bonferronicorrection
nrow(corrdf_br)
corrdf_br$pval <- adjustpval(corrdf_br$pval,nrow(corrdf_br))
colnames(corrdf_br)
corrdf_br$pval
print(nrow(corrdf_br))
str(corrdf_br)

#normalize metric dataframe before computing the master scatter separation matrix
for(i in 33:63){
  #get non-NA values for normalization
  colv <- pbr[,i][!is.na(pbr[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  pbr[,i][!is.na(pbr[,i])] <- colv/msd
  print(paste(names(pbr)[i],msd))
}
pbr_scatter <- pbr[,colnames(pbr)%in%reducedlnMetrics]
dim(pbr)

unique(pbr_scatter[,c('groupnum','groupttl')])
#combine other Neocortex and Not reported Neocortex groups into one
pbr_scatter[pbr_scatter$groupttl == "Control/Neocortex/Not reported", "groupttl"] <- "Control/Neocortex/Others"
pbr_scatter[pbr_scatter$groupnum == 10, "groupnum"] <- 11

pbr_scatter[pbr_scatter$groupttl == "Control/Hippocampus/Others", "groupttl"] <- "Control/Other Hippocampus"
pbr_scatter[pbr_scatter$groupnum == 9, "groupnum"] <- 3

pbr_scatter[pbr_scatter$groupttl== "Control/Hippocampus/Hippocampus/CA1" ,"groupttl"] <-  "Control/Hippocampus/CA1"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Parietal lobe/Somatosensory/Others"  ,"groupttl"] <-  "Control/Other Neocortex/Parietal/S1"
#pbr_scatter[pbr_scatter$groupttl== "Control/Retina/Retina","groupttl"] <-  "Control/Retina"
pbr_scatter[pbr_scatter$groupttl== "Control/Hippocampus/Hippocampus/Others","groupttl"] <-  "Control/Other Hippocampus"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Parietal lobe/Somatosensory/Barrel","groupttl"] <- "Control/Neocortex/Parietal/S1/L4"
#pbr_scatter[pbr_scatter$groupttl== "Control/Olfactory bulb","groupttl"] <-  "Control/Olfactory bulb"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Parietal lobe/Others","groupttl"] <-"Control/Other Neocortex/Parietal"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Occipital lobe","groupttl"] <-  "Control/Neocortex/Occipital"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Frontal lobe/Prefrontal","groupttl"] <-   "Control/Neocortex/Frontal/Prefrontal"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Frontal lobe/Motor","groupttl"] <-   "Control/Neocortex/Frontal/Motor"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Frontal lobe/Others","groupttl"] <-   "Control/Other Neocortex/Frontal"
pbr_scatter[pbr_scatter$groupttl== "Control/Neocortex/Others","groupttl"] <-   "Control/Other Neocortex"
pbr_scatter[pbr_scatter$groupttl== "Control/Others","groupttl"] <-  "Control/Other Brain regions"


crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(N_bifs_total_sum)','Brain regions')
crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(N_tips_total_sum)','Brain regions')
crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(Width_total_sum)','Brain regions')
crosshairscatter(pbr_scatter,'ln(Soma_Surface_total_sum)','ln(Height_total_sum)','Brain regions')
crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(PathDistance_max)','Brain regions')

    
#plot for celltypes also
aggregate(metricDendrites$expercond, by = list(metricDendrites$expercond), length)
#celltype hierarchy
primary <- list(expercond="Control")
x <- getSubset2(usedMetaData, primary, typeMetaData=usedMetaData)
level1 <- "cellclass1"
level2 <- "cellclass2"
level3 <- "cellclass3"
ctMetricHierarchy <- makeHierarchyGroups(primary, metricDendrites, usedMetaData, level1, level2, level3,minSize=100)
printGrps(ctMetricHierarchy)

Grplst <- addgrps2MetricDF(metricDendrites,ctMetricHierarchy)
#Grplst
tot <- 0
for(i in 1:length(Grplst)){
  tot <- tot + length(unlist(Grplst[[i]]$neuronList))
  print(paste(Grplst[[i]]$GroupSpec, length(unlist(Grplst[[i]]$neuronList)), tot))
}

plotDf_ct <- plotTesting(metricDendrites)
plotDf_ct <- metricDf_plotting(plotDf_ct,Grplst)
dim(plotDf_ct)
colnames(plotDf_ct)
unique(plotDf_ct$groupttl)

pct <- subset(plotDf_ct,!is.na(plotDf_ct$groupnum))
colnames(pct)
unique(pct$groupttl)
dim(pct)

corrdf_ct <- leastcorrelated(pct[,33:63],0.80)
#adjust pval with bonferronicorrection
corrdf_ct$pval <- adjustpval(corrdf_ct$pval,nrow(corrdf_ct))
colnames(corrdf_ct)
corrdf_ct$pval
print(nrow(corrdf_ct))
str(corrdf_ct)
#normalize metric dataframe before computing the master scatter separation matrix
for(i in 33:63){
  #get non-NA values for normalization
  colv <- pct[,i][!is.na(pct[,i])]
  msd <- sd(colv)
  print(paste(length(colv),msd))
  pct[,i][!is.na(pct[,i])] <- colv/msd
  print(paste(names(pct)[i],msd))
}
pct_scatter <- pct[,colnames(pct)%in%reducedlnMetrics]
colnames(pct_scatter)

unique(pct_scatter[,c('groupnum','groupttl')])
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Pyramidal cell/Not reported","groupttl"] <- "Control/Principal cell/Pyramidal cell/Others"
pct_scatter[pct_scatter$groupnum==8,"groupnum"] <- 10

pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Not reported/Not reported","groupttl"] <- "Control/Interneuron/Others"
pct_scatter[pct_scatter$groupnum==2,"groupnum"] <- 11

pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Not reported/Others","groupttl"] <- "Control/Interneuron/Others"
pct_scatter[pct_scatter$groupnum==1,"groupnum"] <- 11

pct_scatter[pct_scatter$groupttl=="Control/Others","groupttl"] <- "Control/Interneuron/Others"
pct_scatter[pct_scatter$groupnum==16,"groupnum"] <- 11

pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Ganglion cell/Not reported","groupttl"] <- "Control/Principal cell/Ganglion cell/Others"
pct_scatter[pct_scatter$groupnum==7,"groupnum"] <- 5

pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Granule cell/Not reported","groupttl"] <- "Control/Principal cell/Granule"
pct_scatter[pct_scatter$groupnum==4,"groupnum"] <- 13

pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Granule cell/Others","groupttl"] <- "Control/Principal cell/Granule"
pct_scatter[pct_scatter$groupnum==3,"groupnum"] <- 13

#trim long groupttl names
#pct_scatter[pct_scatter$groupttl=="Control/Axonal terminal/Uniglomerular projection neuron/iPN","groupttl"] <- "Control/Conventional Projections"
#pct_scatter[pct_scatter$groupttl=="Control/Axonal terminal/Uniglomerular projection neuron/Others","groupttl"] <- "Control/Other Uniglomerular Projections"
#pct_scatter[pct_scatter$groupttl=="Control/Axonal terminal/Others","groupttl"] <- "Control/Other Axonal terminals"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Medium spiny cell","groupttl"] <- "Control/Medium spiny"
#pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Martinotti cell","groupttl"] <- "Control/Martinotti"
#pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Not reported/Not reported","groupttl"] <- "Control/Interneuron/Others"

pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Basket cell","groupttl"] <- "Control/Basket"
#pct_scatter[pct_scatter$groupttl=="Control/Granule/Others","groupttl"] <- "Control/Granule"
#pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Granule cell/Not reported","groupttl"] <- "Control/Granule/Others"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Granule","groupttl"] <- "Control/Granule"

pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Pyramidal cell/Magnopyramidal","groupttl"] <- "Control/Magnopyramidal"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Ganglion cell/Monostratified","groupttl"] <- "Control/Ganglion/Monostratified"
#pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Ganglion cell/Not reported","groupttl"] <- "Control/Other Ganglions"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Pyramidal cell/Others","groupttl"] <- "Control/Other Pyramidals"
#pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Pyramidal cell/Not reported","groupttl"] <- "Control/Other Pyramidals"

pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Ganglion cell/Others","groupttl"] <- "Control/Other Ganglions"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Motoneuron","groupttl"] <- "Control/Motoneuron"
pct_scatter[pct_scatter$groupttl=="Control/Principal cell/Others","groupttl"] <- "Control/Other Principals"
pct_scatter[pct_scatter$groupttl=="Control/Interneuron/Others","groupttl"] <- "Control/Other Interneurons"



colnames(pct_scatter)
unique(pct_scatter$groupttl)
unique(pct_scatter[,c('groupttl','groupnum')])
unique(pbr_scatter[,c('groupttl','groupnum')])

dim(pct_scatter)
dim(pbr_scatter)
length(pct_scatter[pct_scatter$groupttl=="Control/Granule",
            "groupttl"])
x11()

pbr_scatter$Bif_ampl_local_avg

#### cross hair plots for species, brain regions and cell types ####
colnames(pbr_scatter)

pbr_scatter[order(pbr_scatter[,"ln(Width_total_sum)"]),][c("groupnum", "groupttl","ln(Width_total_sum)")]

colnames(psp_scatter)
catdepth <- subset(psp_scatter[,c("neuron_name", "ln(Depth_total_sum)")],psp_scatter[,"groupnum"]==13)
catdepth[,2] <- exp(catdepth[,2])
catdepth[order(catdepth[,2]),]
wval <- subset(psp_scatter[,"ln(Depth_total_sum)"],psp_scatter[,"groupnum"]==13)
exp(wval)
range(exp(wval))
median(exp(wval))
quantile(exp(wval),0.25)
quantile(exp(wval),0.75)
mean(exp(wval))
sd(exp(wval))
[c("groupttl","ln(Width_total_sum)")]
p1 <- crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(Width_total_sum)','Brain regions')
p2 <- crosshairscatter(psp_scatter,'ln(Depth_total_sum)','Partition_asymmetry_avg','Species')
p3 <- crosshairscatter(pct_scatter,"ln(Length_total_sum)","ln(Fractal_Dim_avg)",'Arbor types')
p4 <- crosshairscatter(pct_scatter,"Bif_ampl_local_avg","ln(Length_total_sum)",'Arbor types')
#### cross hair scatter for dendrite only data ####
p1 <- crosshairscatter(pbr_scatter,'Bif_ampl_local_avg','ln(Width_total_sum)','Brain regions')

p1
pt1<-ggplot_gtable(ggplot_build(p1))
pt2<-ggplot_gtable(ggplot_build(p2))
pt3<-ggplot_gtable(ggplot_build(p3))
pt1
pt2
l1 <- g_legend(p1)
l2 <- g_legend(p2)
l3 <- g_legend(p3)
grid.arrange(l1,l2)
arrangeGrob(l1,l2)
l1$grobs[[12]]
#print only legends arranged in one column
library(gridExtra)
pdf("legend.pdf")
grid.arrange(arrangeGrob(l1,l2,l3,nrow=1))
dev.off()
#print only plots arranged in one column
pdf("Onlyplots.pdf")
g <- grid.arrange(p1+theme(aspect.ratio=1, legend.position ='none'), 
                  p2+theme(aspect.ratio=1, legend.position = 'none'), 
                  p3+theme(aspect.ratio=1, legend.position = 'none'),
                  ncol=1, nrow=3, widths=c(3/4,3/4,3/4))
dev.off()

grid.arrange(arrangeGrob(p1,p2,p3,ncol=1,
             sub=textGrob("test sub", gp=gpar(font=1))))

#print both plots and legends
pdf("all3plots.pdf")

#gg <- grid.arrange(arrangeGrob(p1+theme(aspect.ratio=1,legend.position = 'none'),l1,p2+theme(aspect.ratio=1,legend.position = 'none'),l2,
#                               p3+theme(aspect.ratio=1,legend.position = 'none'),l3), 
#                   ncol=3, as.table = T, widths=c(7/8,1/8))
gg <- grid.arrange(arrangeGrob(p1+theme(aspect.ratio=1,legend.position = 'none'),l1,nrow=1),
                   arrangeGrob(p2+theme(aspect.ratio=1,legend.position = 'none'),l2,nrow=1),
                   arrangeGrob(p3+theme(aspect.ratio=1,legend.position = 'none'),l3,nrow=1), 
                   ncol=1, as.table = T)
dev.off()


