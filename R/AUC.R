cat("AUC(xb.hat,y)", "\n")

AUC <- function(xb.hat,y){
	n<-length(xb.hat)
	n1<-sum(y)
	mean.rank <- mean(rank(xb.hat)[y == 1])
	AUC<-(mean.rank - (n1 + 1)/2)/(n - n1)
	return(AUC)
}

