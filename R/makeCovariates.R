source("R/randomize.R")

# read counts and offset
ytxt=commandArgs()[5]
ktxt=commandArgs()[6]
#ltxt=commandArgs()[7]

# read tables
Y=read.table(ytxt,as.is=T)
fid=Y[[1]]
Y=as.matrix(Y[,-1])
K=as.matrix(read.table(ktxt,as.is=T)[,-1])
n=ncol(Y)

# feature length
#len=scan(ltxt)

# fpm calculation
fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9

# Singular value decomposition
fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))

# Covariate selection
sf=log(apply(Y,2,sum))
covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
if(cor(sf,covs[,1])^2<0.9){covs=cbind(sf, covs)}

# Write covariates
write.table(covs,col=F,row=F,sep="\t",quote=F,file="data/your.X.txt")


