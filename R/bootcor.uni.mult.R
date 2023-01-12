# Bootcor.uni.mult
# Determine correlation between logistic regression coefficients 
# in univariable and multivariable analysis
# Uses stratified bootstrapping (by outcome)
#
# Author: Ewout Steyerberg, July-Dec 1997

bootcor.uni.mult <- function(fit, B=40, maxit=25, group=fit$y, trim=0.1, save.indices=F,...)
{if(is.null(X <- fit$x) | is.null(Y <- fit$y))
  stop("you did not specify x=T and y=T in the fit")
  
  n     	<- nrow(X)
  cases 	<- sum(Y)
  p     	<- length(fit$coef) - 1
  namec 	<- names(fit$coef[2:(p+1)])
  bootsamp <- matrix(0, nrow=B, ncol=(2*p+1)) 	## to collect boot results
  dimnames(bootsamp)	<- list(NULL,c(paste(namec,"u",sep=""),paste(namec,"m",sep=""),"Shrinkage"))
  
  Y 	<- as.matrix(Y)
  order1	<- order(group)		# stratified bootstrapping 
  X	<- X[order1,]
  
  cat("Bootsample: ")
  for(i in 1:B)
   { j <- c(sample(1:(n-cases-1), replace=T),sample( (n-cases):n, replace=T) )
  ## take bootstrap samples rownumbers
  for (m in 1:p)			## p predictors	
  { f <- lrm.fit(X[j,m,drop=F], Y[j,,drop=F], maxit=maxit) ## univariable coef
  if (!f$fail) 	{bootsamp[i,m] <- f$coefficients[2]}
  else		{bootsamp[i,m] <- NA} }
  
  f <- lrm.fit(X[j,,drop=F], Y[j,,drop=F], maxit=maxit)  ## multivariable coef
  if (!f$fail)	{bootsamp[i,(p+1):(2*p)] <- f$coefficients[2:(p+1)]
  lp		<-  X %*% f$coefficients[2:(p+1)]
  f 		<- lrm.fit(lp, Y, maxit=maxit)
  bootsamp[i,(2*p+1)] <- f$coefficients[2]} ## shrinkage
  else	      { bootsamp[i,(p+1):(2*p+1)] <- NA }
  if (i %% 10==0) cat(i,"", fill=) 
  } # end loop over bootstraps
  
  # Select rows without missings in any column
  for (m in 1:(2*p+1)) {bootsamp <- bootsamp[!is.na(bootsamp[,m]), ] }
  
  # remove outliers in coefficient estimates; close to non-convergence
  require(outliers) 
  out <- outlier(bootsamp[,1:(2*p)], logical=TRUE) 
  indices <- which(rowSums(out) > 0)
  bootsamp <- bootsamp[-indices, ]  
  
  # Calculate correlation
  corcoef <- cor(bootsamp[ ,1    : p   ],
                bootsamp[ ,(p+1):(2*p)]) 
  
  corcoef <- corcoef[row(corcoef) == col(corcoef) ]
  
  # Calculate mean shrinkage factor
  shrink  	<- mean(bootsamp[,(2*p+1)], trim=trim) 
  nB		<- nrow(bootsamp)
  cat("\nNumber of valid bootstraps for correlation and shrinkage: ",nB )
  
  names(corcoef)	<- namec
  names(shrink)	<- "shrinkage"
  names(nB)	<- "B"
  
  if (!save.indices) retlist <- list(r=corcoef,shrink=shrink,B=nB) else
                      retlist <- list(r=corcoef,shrink=shrink,B=nB,Bmatrix=bootsamp)
  
  cat("\nCorrelation coefficients:",round(retlist$r,2),
      "\nShrinkage:", round(retlist$shrink,2),"\n")

  retlist
} # end function bootcor.uni.mult
