source("R/gcCor.R")

# read count and gc content vector (text)
ytxt=commandArgs()[5]
gcc=commandArgs()[6]

# read tables
Y=read.table(ytxt,as.is=T)
fid=Y[[1]]
Y=as.matrix(Y[,-1])
K=array(1, dim(Y))

# GC correction
if(!is.na(gcc)){
	gcc=scan(gcc)
	K=gcCor(Y,gcc)
}

# library size
sf=apply(Y,2,sum)
sf=sf/mean(sf)

# write K as binary
write.table(data.frame(fid, t(t(K)*sf)), col=F, row=F, sep="\t", file="data/your.K.txt", quote=F)

