# This is the code for classical GS (rrBLUP) in each flux. 
# Key: flux prediction via rrBLUP.
# This code is ruuning in R. 
# Model details see below.
# Contact: tong@mpimp-golm.mpg.de
 
## add path and packages
dir <- "/../netGS/"
setwd(dir)

library(rrBLUP)
library(R.matlab)

options(digits = 15)

################################################################################
# # # # # STEP 4 # # # # # 
# # # # # rrBLUP of flux # # # # # 
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
for (r in 1:rr){

print(paste("Replicate ",r," Done!",sep=""))

ids <- idsall[,r]

ypredall <- NULL
corall <- NULL
r2id <- varigid

for (n in r2id){

	#print(paste("Replicate ",r," Reaction ",n," Done!",sep=""))

	ypreda <- rep(NA,nrow(xx))
	Y <- as.numeric(fluxall[n,])

	corr <- NULL
	for (p in 1:f){
		trset <- which(ids!=p)
		teset <- which(ids==p)
		xxtr <- xx[trset,]
		xxte <- xx[teset,]
		y <- Y[trset]
		sol <- mixed.solve(y,Z=xxtr,K=NULL,SE=F)
		ucoef <- as.matrix(sol$u)
		ypredx <- sol$beta + as.numeric(xxte %*% ucoef)
		yp <- Y[teset]
		cor <- cor(yp,ypredx,method="pearson")
		corr <- rbind(corr,cor)
		ypreda[teset] <- ypredx   #### merge to original place
	}
	
	ypredall <- cbind(ypredall,ypreda)
	corall <- cbind(corall,corr)

}

corallall <- rbind(corallall,corall)

# Output: flux prediction value using rrBLUP
write.table(ypredall,paste("fluxpredict_rrBLUP_",r,".csv",sep=""),sep=",",row.names=F,col.names=F)

# Output: flux prediction accuracy using rrBLUP
write.table(corallall,"fluxpredict_rrBLUP_cor.csv",sep=",",row.names=F,col.names=F)

}

################################################################################
# end
################################################################################


