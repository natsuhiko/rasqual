randomize <-
function(x,g=NULL){
	if(is.null(g)){
		n=ncol(x);
		t(apply(x,1,function(xx){xx[order(runif(n))]}))
	}else{
		for(i in unique(g)){
			x[,g==i]=randomize(x[,g==i,drop=F])
		}
		x
	}
}
