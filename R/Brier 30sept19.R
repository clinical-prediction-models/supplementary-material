# Dec 08, Ewout Steyerberg
# function to calculate scaled Brier score

Brier <-function(lp,y, digits=3, pr=T)
{
B     <- mean((y) * (1-plogis(lp))^2 + (1-y) * plogis(lp)^2)
Bmax  <- mean(y) * (1-mean(y))^2 + (1-mean(y)) * mean(y)^2
Bscaled <- 1 - B/Bmax
result  <- c(B,Bmax,Bscaled)
names(result) <- Cs(B,Bmax,Bscaled)
if (pr) {print(result, digits=digits)}
result
}

