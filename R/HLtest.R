# 26-4-00, Yvonne Vergouwe
# Jan 07, Ewout Steyerberg
# function to perform Hosmer-Lemeshow test for external validation
# instead of chi square and corresponding p-value, this function provides the number of subjects per group,
# and the mean values of p and y per group.
##########################
# p	: predicted probability
# Y	: outcome variable
# g	: number of groups to calculate H-L (10 is default)
#
# NB: the library Hmisc need to be attached in order to be able to run hl.ext2
##########################

hl.ext2 <- function(p, y, g = 10, df = g - 1) {
  matres <- matrix(NA, nrow = g, ncol = 5)
  sor <- order(p)
  p <- p[sor]
  y <- y[sor]
  group <- cut2(p, g = g) # g more or less equal sized groups

  len <- tapply(y, group, length) # n per group
  sump <- tapply(p, group, sum) # expected per group
  sumy <- tapply(y, group, sum) # observed per group
  meanp <- tapply(p, group, mean) # mean probability per group
  meany <- tapply(y, group, mean) # mean observed per group
  matres <- cbind(len, meanp, meany, sump, sumy)
  contr <- ((sumy - sump)^2) / (len * meanp * (1 - meanp)) # contribution per group to chi square
  chisqr <- sum(contr) # chi square total
  pval <- 1 - pchisq(chisqr, df) # p-value corresponding to chi square with df degrees of freedom
  cat("\nChi-square", chisqr, " p-value", pval, "\n")
  dimnames(matres) <- list(c(1:g), Cs(n, avg(p), avg(y), Nexp, Nobs))
  result <- list(table(group), matres, chisqr, pval)
}

# function for winsorizing
winsorize <- function(x, lower = quantile(x, probs = 0.01), higher = quantile(x, probs = 0.99)) {
  ifelse(x < lower, lower,
    ifelse(x > higher, higher, x)
  )
}
