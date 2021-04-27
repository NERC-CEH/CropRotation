
###################################
## load the libraries required
###################################

library(shapefiles)
library(reshape2)
library(caret)
library(gdata)
library(data.table)


###################################

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
##############################

all.out <- mapply(fn4,ID=id10k,SIMPLIFY=FALSE)

#######

#find the categories of crops that represent the rows in the TPM and concatenate into single character 
cats=fn12.4(id10k[1])[,1:2]
merge.cats=paste0(cats[,1],cats[,2])


obs.crops=paste0(CM_all[,5],CM_all[,6])

#we now match these observed sequences to the corresponding rows in our TPMs - this is so we know what rows to pull out to make predictions from
mat.cats=match(obs.crops,merge.cats)

CM_all$Obs19=as.factor(CM_all[,7])

##############################
##############################

all.wgt <- mapply(get.wght.LC,ID=1:length(mt.id.10k),SIMPLIFY=FALSE)


output=list()
output[[1]] <- wgt.mat(1)

for(ik in 2:length(all.wgt)){
	if(is.element(SQ10k_LC[ik],SQ10k_LC[1:(ik-1)])){
		lcid <- min(which(SQ10k_LC[1:(ik-1)]==SQ10k_LC[ik]))
		output[[ik]] <- output[[lcid]] 
	}else{
		output[[ik]] <- wgt.mat(ik)
	}

	print(ik)
	
}

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred19LC3_nonGr=factor(prd19,levels=levels(CM_all$Obs19))


##############################
##############################


all.wgt <- mapply(get.wght.All,ID=1:length(mt.id.10k),SIMPLIFY=FALSE)

output=list()
output[[1]] <- wgt.mat(1)

for(ik in 2:length(all.wgt)){
output[[ik]] <- output[[1]] 
}

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred19All3_nonGr=factor(prd19,levels=levels(CM_all$Obs19))


##############################
##############################


output <- all.out

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred1910k3_nonGr=factor(prd19,levels=levels(CM_all$Obs19))



##############################
##############################
##############################
##############################
##############################
##############################

all.out <- mapply(fnr34,ID=id10k,SIMPLIFY=FALSE)

#######

#find the categories of crops that represent the rows in the TPM and concatenate into single character 
cats=fnr312.4(id10k[1])[,1:2]
merge.cats=paste0(cats[,1],cats[,2])


obs.crops=paste0(CM_all[,5],CM_all[,6])

#we now match these observed sequences to the corresponding rows in our TPMs - this is so we know what rows to pull out to make predictions from
mat.cats=match(obs.crops,merge.cats)

CM_all$Obs19=as.factor(CM_all[,7])

##############################
##############################

all.wgt <- mapply(get.wght.LC,ID=1:length(mt.id.10k),SIMPLIFY=FALSE)


output=list()
output[[1]] <- wgt.mat(1)

for(ik in 2:length(all.wgt)){
	if(is.element(SQ10k_LC[ik],SQ10k_LC[1:(ik-1)])){
		lcid <- min(which(SQ10k_LC[1:(ik-1)]==SQ10k_LC[ik]))
		output[[ik]] <- output[[lcid]] 
	}else{
		output[[ik]] <- wgt.mat(ik)
	}

	print(ik)
	
}

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred19LC3_nonGr_r3=factor(prd19,levels=levels(CM_all$Obs19))


##############################
##############################


all.wgt <- mapply(get.wght.All,ID=1:length(mt.id.10k),SIMPLIFY=FALSE)

output=list()
output[[1]] <- wgt.mat(1)

for(ik in 2:length(all.wgt)){
output[[ik]] <- output[[1]] 
}

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred19All3_nonGr_r3=factor(prd19,levels=levels(CM_all$Obs19))


##############################
##############################


output <- all.out

##now run the prediction function using the estimated spatially smoothed TPMs
prd19 <- predfor(output=output,CM_all=CM_all,id10k=id10k,mat.cats=mat.cats)

##convert the observed and predicted classes to factors with the same set of levels
CM_all$pred1910k3_nonGr_r3=factor(prd19,levels=levels(CM_all$Obs19))



##############################
##############################


