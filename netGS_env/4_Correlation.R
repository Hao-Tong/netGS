# This is the code for checking correlation coefficient of biomass. 
# This code is ruuning in R. 
# Model details see below.
# Contact: tong@mpimp-golm.mpg.de
 
dir <- "/../netGS_env/" 
setwd(dir)

################################################################################
# # # # # STEP 4 # # # # # 
# # # # # correlation coefficient of biomass # # # # # 
################################################################################

# read measured biomass
biom <- read.table("biomass_optNlowN.csv",sep=",",header=F)
biom1 <- biom[,2]  ## high N
biom1 <- matrix(biom1,length(biom1),1)
biom2 <- biom[,3]  ## low N
biom2 <- matrix(biom2,length(biom2),1)

idsall <- read.table("foldid.csv",sep=",",header=F)

corr0 <- NULL
for (t in 1:50){
for (p in 1:3){
	ff <- read.table(paste("biomasspredict_lowN_r",t,"_f",p,".csv",sep=""),sep=",",header=F)[336,]
	ff <- as.numeric(as.matrix(ff))

	fid <- idsall[,t]
	ids <- which(fid==p)
	bioms <- biom2[ids,]

	corr <- cor(bioms,ff,use="complete.obs") 

	corr0 <- rbind(corr0,corr)
}
}

corm <- mean(corr0[,1],na.rm=T)
corr0 <- rbind(corr0,corm)

write.table(corr0,"biomcorr_lowN.csv",sep=",",row.names=F,col.names=F)

################################################################################
# end
################################################################################


