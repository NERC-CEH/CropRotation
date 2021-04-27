
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
library(data.table)


###################################
## 
###################################

ite <- read.csv("ite.csv")

ite.lst <- split(ite,ite$GRIDREF)

ite.dom <- rbindlist(lapply(ite.lst,function(x){return(x[which.max(x$pc_ite),])}))


mindist.id <- function(x0,y0,x1,y1){

	dst <- sqrt((x0-x1)^2 + (y0-y1)^2)

	return(which.min(dst))

}


ite.id=c()

for(i.10k in 1:length(cent10km[,1])){
	
	ite.id[i.10k] <- mindist.id(x0=cent10km$x[i.10k],y0=cent10km$y[i.10k],x1=ite.dom$east.ll,y1=ite.dom$north.ll) 

}


SQ10k_LC <- ite.dom$LC2007[(ite.id[mt.id.10k])]



##############################

## main set of functions to produce spatially smoothed Transition Proability Matrices (TPMs). 

##############################

## Get the Transition Probability Matrices (TPMs) 
#apply the TPM function to each 10km cell. this returns a list whereby each entry is a transition probability matrix for that 10km cell

all.out <- mapply(fn4,ID=id10k,SIMPLIFY=FALSE)


## Get the weights defined according to distances between 10km cells
#apply the weighting function to each 10km cell. this returns a list whereby each entry is the weighting assigned to other 10km cells, given their proximiity to the 10km cell represented by the list entry

all.wgt <- mapply(get.wght.LC,ID=(id10k+1),SIMPLIFY=FALSE)


## Calculate a weighted average of the TPMS  
#apply the matrix weighted average function to each 10km cell
#### NOTE THAT THIS CAN TAKE A WHILE TO RUN!! #####

output=mapply(wgt.mat,idx=1:length(all.wgt),SIMPLIFY=FALSE)



#######

#find the categories of crops that represent the rows in the TPM and concatenate into single character 
cats=fn12.4(id10k[1])[,1:3]
merge.cats=paste0(cats[,1],cats[,2],cats[,3])



##############################

### FOR PREDICTION ###

##############################

#to predict the crop assigned at time t+1 we need the observed crop at time t and time  t-1
#in this example we pull out the 2017 and 2018 crops to predict 2019 values
obs.crops=paste0(CM_all[,4],CM_all[,5],CM_all[,6])

#we now match these observed sequences to the corresponding rows in our TPMs - this is so we know what rows to pull out to make predictions from
mat.cats=match(obs.crops,merge.cats)


##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)


##convert the observed and predicted classes to factors with the same set of levels
CM_all$Obs19=as.factor(CM_all[,7])
CM_all$pred19=factor(prd19,levels=levels(CM_all$Obs19))

##produce confusion matrix to calculate the error in the predicted set. 
confusionMatrix(CM_all$pred19,CM_all$Obs19)


#######################
#######################


##seperate observation and predictions into a list by 10km cell
CM_list <- split(CM_all,CM_all$ID10k)

##apply confusion matrix to each 10km cell to get a spatially explicit accuracy measure
Acc_10k <- unlist(lapply(CM_list,function(x){return(as.numeric(confusionMatrix(x$pred19,x$Obs19)$overall[1]))}))

##write this accuracy measure back to the dbf file so that this can be viewed easily in arcmap/qgis etc. 
fishnet$dbf$dbf$Accuracy=NA
fishnet$dbf$dbf$Accuracy[match(as.integer(names(Acc_10k)),0:3704)]=unlist(Acc_10k)

write.dbf(fishnet$dbf,out.name="Fish10kmPoly_Acc.dbf")

######
#still need to add in the backward prediction, though have a working toy example for this. 



