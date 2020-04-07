# This is the code for classical GS (BayesC) in each flux. 
# Key: flux prediction via BayesC.
# This code is ruuning in R. 
# Model details see below.
# Contact: tong@mpimp-golm.mpg.de
  
## add path and packages
dir <- "/../netGS/"
setwd(dir)

library(BGLR)
library(R.matlab)

options(digits = 15)

################################################################################
# # # # # STEP 4 # # # # # 
# # # # # BayesC of flux # # # # # 
################################################################################

# read fluxes
datamat <- readMat("fluxgenotype.mat")
fluxall <- datamat$fluxall[1:336,]

# read genotypic data
xxall <- read.table("snp.csv",sep=",",header=F)
xx <- as.matrix(xxall)

varigid <- c(1:336)

# read cross-validation fold id
idsall <- read.table("foldid.csv",sep=",",header=F)
rr <- 50 #replicate number
f <- 3 #fold number

ypredallall <- NULL
corallall <- NULL

# BayesC parameters
nIter=1500;
burnIn=500;
thin=3;
saveAt="3e";
S0=NULL;
weights=NULL;
R2=0.5;
verbose=FALSE;
ETA<-list(list(X=xx,model='BayesC'))
nn <- nrow(xx)

for (r in 1:rr){

print(paste("Replicate ",r," Done!",sep=""))

ids <- idsall[,r]

ypredall <- NULL
corall <- NULL
r2id <- varigid

for (n in r2id){

	#print(paste("Replicate ",r," Reaction ",n," Done!",sep=""))

	ypreda <- rep(NA,nn)
	Y <- as.numeric(fluxall[n,])

	corr <- NULL
	for (p in 1:f){
		trset <- which(ids!=p)
		teset <- which(ids==p)
		y <- rep(NA,nn)
		y[trset] <- Y[trset]
		fitbc <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2,verbose=verbose)
		ypredxall <- fitbc$yHat
		ypredx <- ypredxall[teset]
		yp <- Y[teset]
		cor <- cor(yp,ypredx,method="pearson")
		corr <- rbind(corr,cor)
		ypreda[teset] <- ypredx   #### merge to original place
	}
	
	ypredall <- cbind(ypredall,ypreda)
	corall <- cbind(corall,corr)

}

corallall <- rbind(corallall,corall)

# Output: flux prediction value using BayesC
write.table(ypredall,paste("fluxpredict_BayesC_",r,".csv",sep=""),sep=",",row.names=F,col.names=F)

# Output: flux prediction accuracy using BayesC
write.table(corallall,"fluxpredict_BayesC_cor.csv",sep=",",row.names=F,col.names=F)

}

################################################################################
# end
################################################################################


