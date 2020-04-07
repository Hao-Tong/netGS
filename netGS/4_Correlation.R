# This is the code for checking correlation coefficient of biomass. 
# This code is ruuning in R. 
# Model details see below.
# Contact: tong@mpimp-golm.mpg.de
 
dir <- "/../netGS/"
setwd(dir)

################################################################################
# # # # # STEP 6 # # # # # 
# # # # # correlation coefficient of biomass # # # # # 
################################################################################

# read measured biomass
biom <- read.table("biomass_optN.csv",sep=",",header=F)[,2]
idsall <- read.table("foldid.csv",sep=",",header=F)

corr0 <- NULL
for (t in 1:50){
for (p in 1:3){
	ff <- read.table(paste("biomasspredict_r",t,"_f",p,".csv",sep=""),sep=",",header=F)[336,]
	ff <- as.numeric(as.matrix(ff))

	fid <- idsall[,t]
	ids <- which(fid==p)
	bioms <- biom[ids]

	corr <- cor(bioms,ff) 

	corr0 <- rbind(corr0,corr)
}
}

corm <- mean(corr0[,1],na.rm=T)
corr0 <- rbind(corr0,corm)

write.table(corr0,"biomcorr.csv",sep=",",row.names=F,col.names=F)

################################################################################
# end
################################################################################


