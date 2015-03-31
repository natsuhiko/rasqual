
#source("~/Project/R/Util/gcCor.R")
source("~/Project/R/Util/randomize.R")



gene=read.table("../Peaks/peaks.bed.gz",as.is=T)
gene=gene[gene[[1]]%in%1:22,]
gcc=gene[[8]]
len=gene[[7]]

Y=matrix(scan("fc.gz"),length(len))

K=gcCor(Y,gcc,T)

browser()

sf=apply(Y,2,sum)
sf=sf/mean(sf)

fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))/len*1e9
hoge=svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
d=svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))$d
covs=cbind(log(apply(Y/K,2,sum)),hoge$v[,c(1:3)])
fpkm2=exp(t(resid(lm(t(log(fpkm))~covs)))+apply(log(fpkm),1,mean))

f=file("log.fpkm.bin","wb")
writeBin(c(t(log(fpkm2))),f)
close(f)

f=file("Y.bin","wb")
writeBin(c(t(Y)),f)
close(f)

f=file("K.bin","wb")
writeBin(c(t(K)*sf),f)
close(f)

f=file("Cov.bin","wb")
writeBin(c(covs),f)
close(f)

f=file("sf.bin","wb")
writeBin(as.double(apply(Y,2,sum)),f)
close(f)


#f=file("Y.bin","wb")
#writeBin(c(t(Y)),f)
#close(f)



