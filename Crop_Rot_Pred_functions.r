##############################


#a function to calculate the centroid of a 10km square/rectangular polygon, when it is specified as a data frame of points
cent <- function(df){

	mx = min(df[,1])+5000
	my = min(df[,2])+5000

	return(c(mx,my))

}


##############################

#a function to calculate the centroid of a 10km square/rectangular polygon, when it is specified as a data frame of points
avg_coord <- function(df){

	mx = mean(df[,1])
	my = mean(df[,2])

	return(c(mx,my))

}

##############################


#function to calculate 3d transition probability matrix from the frequencies of observations
fn <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
	
	e.res <- embed(unmatrix(res,byrow=T),3)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]

	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k]))}

	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:3, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	outdat=outdat[,-(1:2)]
	
	outdat
}

##############################


fn12 <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
	
	e.res <- embed(unmatrix(res,byrow=T),3)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]

	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k]))}

	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:3, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	#outdat=outdat[,-(1:2)]
	
	outdat
}


###########################
###########################


#function to calculate 3d transition probability matrix from the frequencies of observations - allowing for 4 year rotations
fn4 <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
		
	e.res <- embed(unmatrix(res,byrow=T),4)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]
	
	
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,NA,all.crps[k]))}
	
	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2", "Ym3")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:4, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym3 + Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	outdat=outdat[,-(1:3)]
	
	outdat
}


############################
############################

fn12.4 <- function(ID){

		#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
		
	e.res <- embed(unmatrix(res,byrow=T),4)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]
	
	
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,NA,all.crps[k]))}
	
	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2", "Ym3")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:4, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym3 + Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	#outdat=outdat[,-(1:3)]
	
	outdat

}


##############################
##############################


#use a normal density kernel to determine weights applied to each corresponding distsance pairing
get.wght=function(ID){

	x=as.numeric(dist.mat[,ID])
	wght=dnorm(x/mult)

	wght=wght[mt.id.10k]

	wght
	
}



#use a classification to determine averaging across sets
get.wght.LC=function(ID){

	x=rep(0,length(mt.id.10k))
	li=which(SQ10k_LC==SQ10k_LC[ID])
	x[li]=1
	
	wght=x

	wght
	
}

#assign all weights equal - equivalent to full average
get.wght.All=function(ID){

	wght=rep(1,length(mt.id.10k))

	wght
	
}


##############################


#use the estimates TPM and weights to calculate a weighted average TPM 
wgt.mat = function(idx){

	mat.by.wgts = mapply(function(x,v){x*v},x=all.out,all.wgt[[idx]],SIMPLIFY=FALSE)

	test.out = Reduce('+',mat.by.wgts)

	rsum = apply(test.out,1,sum)

	scaled.out=diag(1/rsum)%*%as.matrix(test.out)

	scaled.out[is.na(scaled.out)]=0

	scaled.out

}


#############################

#function for estimating the next crop in sequence (hence predict forward) given the previous 2 crops
predfor <- function(output,CM_all,id10k,mat.cats){

	CROPS <- colnames(output[[1]])

	#create empty vector to store predictions in
	pred.crop <- rep(NA, length(CM_all[,1]))

	for(k in 1:length(id10k)){

		#which individual polygons are in the current 10km cell in iteration
		vid=which(CM_all$ID10k==id10k[k])

		#use this simple version if just choosing crop with highest probability
		#pred.crop[vid]=CROPS[apply(as.matrix(output[[k]][mat.cats[vid],]),1,which.max)]

		output[[k]]=as.matrix(output[[k]])
		
		#check the probabilities sum to 1 and where there is an empty row (all 0s) recode to equal 1 so code doesnt fall over
		chk=apply(output[[k]],1,sum)
		if(any(chk==0)){output[[k]][chk==0,]=1}
		
		#extract corresponding rows of TPM
		cat.mat = output[[k]][mat.cats[vid],]
		
		
		#if there is only 1 row extracted, then just sample from this, probabilistically sampling the next crop type
		if(is.null(dim(cat.mat))){
			if(any(is.na(cat.mat))){pred.crop[vid]=NA}else{
			pred.crop[vid]=CROPS[sample(1:length(CROPS),1,prob=as.numeric(cat.mat))]
			}
			
		#if multiple row are extract from TPM then sample across each row iteratively, probabilistically sampling the next crop type	
		}else{
		
			na.id.rm <- which(apply(cat.mat,1,function(x){any(is.na(x))}))
			if(length(na.id.rm)>0){
				pred.crop[vid[na.id.rm]]=NA
				cat.mat=cat.mat[-na.id.rm,]
				vid=vid[-na.id.rm]
			}
			if(length(vid)>length(na.id.rm)){
			if(is.null(dim(cat.mat))){
		
			pred.crop[vid]=CROPS[sample(1:length(CROPS),1,prob=as.numeric(cat.mat))]

			#if multiple row are extract from TPM then sample across each row iteratively, probabilistically sampling the next crop type	
			}else{pred.crop[vid]=CROPS[apply(cat.mat,1,function(x){sample(1:length(CROPS),1,prob=as.numeric(x))})]
				#pred.crop[vid]=CROPS[apply(cat.mat,1,function(x){which.max(x)})]
			}}
		}
		
		#print out iteration of loop for users info
		print(k)
		
	}
	
	return(pred.crop)
	
}

#############################

#function for estimating the previous crop in sequence (hence predict backwards) given the following 2 crops
predback <- function(output,CM_all,id10k,mat.cats){

##still to do

}





#function to calculate 3d transition probability matrix from the frequencies of observations
fnr3 <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),5:7])
	res=cbind(res,rep(NA,length(res[,1])))
	
	e.res <- embed(unmatrix(res,byrow=T),3)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]

	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k]))}

	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:3, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	outdat=outdat[,-(1:2)]
	
	outdat
}

##############################


fnr312 <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),5:7])
	res=cbind(res,rep(NA,length(res[,1])))
	
	e.res <- embed(unmatrix(res,byrow=T),3)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]

	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k]))}

	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:3, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	#outdat=outdat[,-(1:2)]
	
	outdat
}


###########################
#####


#function to calculate 3d transition probability matrix from the frequencies of observations - allowing for 4 year rotations
fnr34 <- function(ID){

	#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
		
	e.res <- embed(unmatrix(res,byrow=T),4)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]
	
	
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,NA,all.crps[k]))}
	
	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2", "Ym3")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:4, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym3 + Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	outdat=outdat[,-(1:3)]
	
	outdat
}


############################
############################

fnr312.4 <- function(ID){

		#subset the relevant 10km square
	res=as.matrix(CM_all[CM_all$ID10k==ID & !is.na(CM_all$ID10k),3:7])
	res=cbind(res,rep(NA,length(res[,1])))
		
	e.res <- embed(unmatrix(res,byrow=T),4)
	e.res <- e.res[-which(apply(e.res,1,function(x){any(is.na(x) | x=="NA")})),]
	
	
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(all.crps[k],NA,NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,all.crps[k],NA,NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,all.crps[k],NA))}
	for(k in 1:length(all.crps)){e.res=rbind(e.res,c(NA,NA,NA,all.crps[k]))}
	
	
	colnames(e.res) <- c("Y0", "Ym1", "Ym2", "Ym3")
	 
	# Count the transitions
	counts <- table( as.data.frame(e.res) )
	 
	# Divide by the total number of transitions, to have probabilities
	probabilities <- counts
	probabilities[] <- as.vector(counts) / rep( as.vector(apply( counts, 2:4, sum )), each=dim(counts)[1] )

	# Check that the probabilities sum up to 1
	#apply( probabilities, 2:3, sum )

	# Convert the 3-dimensional array to a data.frame
	outdat=dcast( melt( probabilities ), Ym3 + Ym2 + `Ym1` ~ Y0 )
	
	#recode any NAs as 0s
	outdat[is.na(outdat)]=0
	
	#outdat=outdat[,-(1:3)]
	
	outdat

}

