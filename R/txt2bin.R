
# read count and offset files (text)
ytxt=commandArgs()[5]
ktxt=commandArgs()[6]
# read tables
y=read.table(ytxt,as.is=T)
k=read.table(ktxt,as.is=T)
# output binary file names
ybin=gsub("txt", "bin", ytxt)
kbin=gsub("txt", "bin", ktxt)
# write tables as binary
writeBin(as.double(c(t(y[,-1]))), file(ybin,"wb"))
writeBin(as.double(c(t(k[,-1]))), file(kbin,"wb"))
