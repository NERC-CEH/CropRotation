################################################
###### ASSIST - Predicting Crop Rotations ###### 
################################################


#this code relies on 3 files: Fish10kmPoly ; crop_rotation_2015_2019 ; Ident_Crop10km


###################################
## load the libraries required
###################################

library(shapefiles)
library(reshape2)
library(caret)
library(gdata)


###################################
## 
###################################


#load in the 10km fishnet created in ArcGIS
fishnet <- read.shapefile("Fish10kmPoly") 

#apply centroid function to fishnet to calculate the centre point of each 10km cell
cent10km <- data.frame(matrix(unlist(lapply(fishnet$shp$shp,function(x){cent(x$points)})),ncol=2,nrow=length(fishnet$shp$shp),byrow=T))
names(cent10km)=c("x","y")

#calculate the distances between each pair of 10mk centroids
dist.mat=as.matrix(dist(cent10km,diag=TRUE))

#read in each of the crop maps produced thus far
CM_all <- read.dbf("crop_rotation_2015_2019.dbf")$dbf

#read in info on the 10km cells and their IDs - this was an overlay between the gid's and the 10km fishnet done in ArcGIS
fishID <- read.dbf("Ident_Crop10km.dbf") 

#match the ID of the polygons in the crops data to those in the fishnet
mt_id <- match(CM_all[,2],fishID$dbf$gid)
 
#pull out corresponding fishnet (eg 10km cell) id to assign to each polygon 
CM_all$ID10k <- fishID$dbf$FID_Fish10[mt_id]
 
#filter out those field parcels without crop data up to 2019
CM_all <- CM_all[CM_all$recent_yea==2019,]
 
 
#split the rotation column into individual years 
rots <- mapply(function(x){unlist(strsplit(as.character(x),"_"))},x=CM_all$rotation,SIMPLIFY=FALSE)

#check they all have length 5
summary(unlist(lapply(rots,length)))
 
#split the rotation character into constituent crops for each year.  
CM_all$Crop2015 <- (unlist(lapply(rots,function(x){return(x[1])})))
CM_all$Crop2016 <- (unlist(lapply(rots,function(x){return(x[2])})))
CM_all$Crop2017 <- (unlist(lapply(rots,function(x){return(x[3])})))
CM_all$Crop2018 <- (unlist(lapply(rots,function(x){return(x[4])})))
CM_all$Crop2019 <- (unlist(lapply(rots,function(x){return(x[5])})))
 
 
##################
##################

#remove the entries that do not have an associated 10km grid id. These NAs are because the fishnet overlay needs to be redone with the new collated 2015-2019 dataset
#once redone, this step will not be needed
CM_all <- CM_all[-which(is.na(CM_all$ID10k)),]

##################
################## 

 
#tidy up data frame  
CM_all <- CM_all[,c(2,5,10:14,9)]

#find all the crop types
all.crps <- unique(c(as.character(unique(CM_all[,3])),as.character(unique(CM_all[,4])),as.character(unique(CM_all[,5])),as.character(unique(CM_all[,6])),as.character(unique(CM_all[,7]))))
	
#remove any that have been assigned as NA - we will just ignore these for now	
if(any(all.crps=="NA")){
	all.crps <- all.crps[-which(all.crps=="NA")]
}


#extract the IDs that represent the different 10km cells
id10k=unique(CM_all$ID10k)

#match these to the full set of all 10km cells across the UK that are represented in the full fishnet
mt.id.10k = match(id10k,0:3704)



#################
#################

# this sets the bandwidth for the spatial smoothing of the 10km TPMs - it helps convert the distances between 10km cells to weights which will be used
# in a weighted average of TPMs

mult=25000/0.674

#################
#################


##############################

## main set of functions to produce spatially smoothed Transition Proability Matrices (TPMs). 

##############################

## Get the Transition Probability Matrices (TPMs) 
#apply the TPM function to each 10km cell. this returns a list whereby each entry is a transition probability matrix for that 10km cell

all.out <- mapply(fn,ID=id10k,SIMPLIFY=FALSE)


## Get the weights defined according to distances between 10km cells
#apply the weighting function to each 10km cell. this returns a list whereby each entry is the weighting assigned to other 10km cells, given their proximiity to the 10km cell represented by the list entry

all.wgt <- mapply(get.wght,ID=(id10k+1),SIMPLIFY=FALSE)


## Calculate a weighted average of the TPMS  
#apply the matrix weighted average function to each 10km cell
#### NOTE THAT THIS CAN TAKE A WHILE TO RUN!! #####

output=mapply(wgt.mat,idx=1:length(all.wgt),SIMPLIFY=FALSE)



#######

#find the categories of crops that represent the rows in the TPM and concatenate into single character 
cats=fn12(id10k[1])[,1:2]
merge.cats=paste0(cats[,1],cats[,2])



##############################

### FOR PREDICTION ###

##############################

#to predict the crop assigned at time t+1 we need the observed crop at time t and time  t-1
#in this example we pull out the 2017 and 2018 crops to predict 2019 values
obs.crops=paste0(CM_all[,5],CM_all[,6])

#we now match these observed sequences to the corresponding rows in our TPMs - this is so we know what rows to pull out to make predictions from
mat.cats=match(obs.crops,merge.cats)


##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)


##convert the observed and predicted classes to factors with the same set of levels
CM_all$Obs19=as.factor(CM_all[,7])
CM_all$pred19=factor(prd19,levels=levels(CM_all$Obs19))

##produce confusion matrix to calculate the error in the predicted set. 
confusionMatrix(CM_all$pred19,CM_all$Obs19)


######
#still need to add in the backward prediction, though have a working toy example for this. 



