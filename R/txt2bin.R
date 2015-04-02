
ytxt=commandArgs()[5]
ktxt=commandArgs()[6]

y=read.table(ytxt,as.is=T)
k=read.table(ktxt,as.is=T)

ybin=gsub("txt", "bin", ytxt)
kbin=gsub("txt", "bin", ktxt)

writeBin(as.double(c(t(y[,-1]))), file(ybin,"wb"))
writeBin(as.double(c(t(k[,-1]))), file(kbin,"wb"))


