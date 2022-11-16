# Yvonne Vergouwe, 06.05.02
# Ewout Steyerberg, June 2007: 
# for graphs with calibration, discrimination, and clinical usefulness

# lp to p

InvLogit	<-function(lp)
{	p	<-1/(1+exp(-lp)) }

####
#				        ter			  nec
# p(ter)>cut		nac[2]		nbd[2]
# p(ter)<cut		nac[1]		nbd[1]
#

# se,sp, for cutoff value, e.g. 0.7
sesp	<-function(p,y,cutoff=.3)
{
tab	<-matrix(NA,1,13)	
	d	<- cut(p,c(0,cutoff,1))
	n01	<-as.vector(tapply(y,d,length))	# n in 2 groups
	nac	<-as.vector(tapply(y,d,sum))		# events in 2 groups
	nbd	<-n01-nac								# non-events in 2 groups
	acc			<-(nac[2]+nbd[1])/length(y)
	err			<-1-acc
	sens		<-nac[2]/(sum(nac))			#sens
	spec		<-nbd[1]/(sum(nbd))			#spec
	ppv			<-nac[2]/(nac[2]+nbd[2])	#ppv
	npv			<-nbd[1]/(nbd[1]+nac[1])	#npv
	you			<-sens+spec-1
	psep		<-ppv+npv-1
	or			<-(nac[2]*nbd[1])/(nac[1]*nbd[2])	
	w       <- cutoff / (1-cutoff)
  NBmodel <- (nac[2] - w * nbd[2]) / length(y)  # NB = (TP - w FP) / N
  NBtreat <-  (nac[1]+nac[2] - w*(nbd[1]+nbd[2])) / length(y)
  NBnotreat <- 0
  NB      <- NBmodel - (max(NBtreat,0))

	tab[1,]<-c(acc,err,sens,spec,ppv,npv,you,psep,or,NBmodel, NBtreat, NBnotreat, NB)
	dimnames(tab)<-	list(c(cutoff),c("acc","err","sens","spec","ppv","npv","youden","PSEP","OR",
                    "NBmodel", "NBtreat", "NBnotreat", "NB"))
	tab
}

#################
#functie val.prob aangepast:
#- p verdeling opgesplitst naar uitkomst
#- 0 en 1 om aan te geven  welke uitkomst onder/boven lijn staat. je kan in commando andere text opgeven:
#  d1lab="..", d0lab=".."
#- y as: "Observed Frequency"
#- driehoek: "Grouped patients"
#- een horizontale cut-off lijn kan in de plot gezet worden, je hoeft alleen de y-co?rdinaat op te geven:
#  bijv. cutoff=0.7, NB: stijl van lijnen enigszins gewijzigd.

# Jan 2019: add colors

val.prob.yvon<-
function(p, y, logit, group, weights = rep(1, length(y)), normwt = F, pl = T, 
	smooth = T, logistic.cal = F, xlab = "Predicted Probability", ylab = 
	"Observed Frequency", xlim = c(-0.02, 1),ylim = c(-0.15,1), m, g, cuts, emax.lim = c(0, 1), 
	legendloc =  c(0.55 , 0.27), statloc = c(0,1),
	riskdist = "predicted", cex = 0.75, mkh = 0.02, connect.group = 
	F, connect.smooth = T, g.group = 4, evaluate = 100, nmin = 0, d0lab="0", d1lab="1",cutoff,
	col.point='#00468B', col.logcal='#00468B', col1='#ED0000', col0='#42B540')
{
	if(missing(p))
		p <- 1/(1 + exp( - logit))
	else logit <- log(p/(1 - p))
	if(length(p) != length(y))
		stop("lengths of p or logit and y do not agree")
	names(p) <- names(y) <- names(logit) <- NULL
	if(!missing(group)) {
		if(length(group) == 1 && is.logical(group) && group)
			group <- rep("", length(y))
		if(!is.factor(group))
			group <- if(is.logical(group) || is.character(group)) 
				  as.factor(group) else cut2(group, g = 
				  g.group)
		names(group) <- NULL
		nma <- !(is.na(p + y + weights) | is.na(group))
		ng <- length(levels(group))
	}
	else {
		nma <- !is.na(p + y + weights)
		ng <- 0
	}
	logit <- logit[nma]
	y <- y[nma]
	p <- p[nma]
	if(ng > 0) {
		group <- group[nma]
		weights <- weights[nma]
		return(val.probg(p, y, group, evaluate, weights, normwt, nmin)
			)
	}
	if(length(unique(p)) == 1) {
#22Sep94
		P <- mean(y)
		Intc <- log(P/(1 - P))
		n <- length(y)
		D <- -1/n
		L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
		L.cal <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = T)
		U.chisq <- L01 - L.cal
		U.p <- 1 - pchisq(U.chisq, 1)
		U <- (U.chisq - 1)/n
		Q <- D - U

		stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y - p[
			1])^2), Intc, 0, rep(abs(p[1] - P), 2))
		names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq", 
			"D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", 
			"Intercept", "Slope", "Emax", "Eavg")
		return(stats)
	}
	i <- !is.infinite(logit)
	nm <- sum(!i)
	if(nm > 0)
		warning(paste(nm, 
			"observations deleted from logistic calibration due to probs. of 0 or 1"
			))
	f <- lrm.fit(logit[i], y[i])
	f2<-	lrm.fit(offset=logit[i], y=y[i])
	stats <- f$stats
	n <- stats["Obs"]
	predprob <- seq(emax.lim[1], emax.lim[2], by = 0.0005)
	lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
	calp <- 1/(1 + exp( - lt))
	emax <- max(abs(predprob - calp))
	if(pl) {
		plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = "",ylab = "")
		mtext(xlab, side=1, line=2.5)
    mtext(ylab, side=2, line=2.5)
		abline(0, 1, lty = 3)
		lt <- 3
		leg <- "Ideal"
		marks <- -1
		if(logistic.cal) {
			lt <- c(lt, 1)
			leg <- c(leg, "Logistic calibration")
			marks <- c(marks, -1)
		}
		if(smooth) {
			Sm <- lowess(p, y, iter = 0)
			if(connect.smooth) {
				lines(Sm, lty = 1)
				lt <- c(lt, 1)
				marks <- c(marks, -1)
			}
			else {
				points(Sm)
				lt <- c(lt, 0)
				marks <- c(marks, 1)
			}
			leg <- c(leg, "Nonparametric")
			cal.smooth <- approx(Sm, xout = p)$y
			eavg <- mean(abs(p - cal.smooth))
		}
		if(!missing(m) | !missing(g) | !missing(cuts)) {
			if(!missing(m))
				q <- cut2(p, m = m, levels.mean = T, digits = 
				  7)
			else if(!missing(g))
				q <- cut2(p, g = g, levels.mean = T, digits = 
				  7)
			else if(!missing(cuts))
				q <- cut2(p, cuts = cuts, levels.mean = T, 
				  digits = 7)
			means <- as.single(levels(q))
			prop <- tapply(y, q, function(x)
			mean(x, na.rm = T))
			points(means, prop, pch = 2, col=col.point, lwd=2)
			if(connect.group) {
				lines(means, prop)
				lt <- c(lt, 1)
			}
			else lt <- c(lt, 0)
			leg <- c(leg, "Grouped patients")
			marks <- c(marks, 2)
		}
	}
	lr <- stats["Model L.R."]
	p.lr <- stats["P"]
	D <- (lr - 1)/n
	L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
	U.chisq <- L01 - f$deviance[2]
	p.U <- 1 - pchisq(U.chisq, 2)
	U <- (U.chisq - 2)/n
	Q <- D - U
	Dxy <- stats["Dxy"]
	C <- stats["C"]

#several R2 measures set to zero June 2007 (nopt necessary now)
	n		<-length(y)
	r2.pear		<-0
	r2.nagel	<-0
	r2.harr	 <- 0

	B <- sum((p - y)^2)/n
		stats <- c(Dxy, C, r2.pear,r2.nagel,r2.harr, D, lr, p.lr, U, U.chisq, p.U, Q, B, f2$coef[1],f$coef[2], 
		emax)
	names(stats) <- c("Dxy", "C (ROC)", "R2 Pearson", "R2 Nagelkerke", "R2 Harrell", "D", "D:Chi-sq", "D:p", "U", 
		"U:Chi-sq", "U:p", "Q", "Brier", "Int|slope=1", "Slope", "Emax")
	if(smooth)
		stats <- c(stats, c(Eavg = eavg))
	if(!missing(cutoff)) {
			lines(c(cutoff),c(1),lty=2,type="h")
			lt <- c(lt, 2)
			leg <- c(leg, Cs(Cut-off))
		}	
	if(pl) {
		logit <- seq(-7, 7, length = 200)
		prob <- 1/(1 + exp( - logit))
		pred.prob <- f$coef[1] + f$coef[2] * logit
		pred.prob <- 1/(1 + exp( - pred.prob))
		if(logistic.cal) lines(prob, pred.prob, lty = 1, col=col.logcal, lwd=2)	
	#	pc <- rep(" ", length(lt))
#	pc[lt==0] <- "."
		lp <- legendloc
		if(!is.logical(lp)) {
			if(!is.list(lp))
				lp <- list(x = lp[1], y = lp[2])
			legend(lp, leg, lty = lt, marks = marks, mkh = mkh, 
				cex = cex, bty = "n")	#, pch=pc)
		}
		if(!is.logical(statloc)) {
			dostats <- c(1, 2, 3, 4, 7, 10, 11, 12, 13, 14)
			leg <- format(names(stats)[dostats])	#constant length
			leg <- paste(leg, ":", format(stats[dostats]), sep = 
				"")
			if(!is.list(statloc))
				statloc <- list(x = statloc[1], y = statloc[2]
				  )
			text(statloc, paste(format(names(stats[dostats])), 
				collapse = "\n"), adj = 0, cex = cex)
			text(statloc$x + 0.225 , statloc$y, paste(
				format(round(stats[dostats], 3)), collapse = 
				"\n"), adj = 1, cex = cex)	
	#	legend(statloc, leg, lty=rep(0, length(dostats)))
		}
		if(is.character(riskdist)) {
			if(riskdist == "calibrated") {
				x <- f$coef[1] + f$coef[2] * log(p/(1 - p))
				x <- 1/(1 + exp( - x))
				x[p == 0] <- 0
				x[p == 1] <- 1
			}
			else x <- p
			bins <- seq(0, 1, length = 99)
			x <- x[x >= 0 & x <= 1]
#08.04.01,yvon: verdeling van predicted opgesplitst naar uitkomst
			f0	<-table(cut(x[y==0],bins))
			f1	<-table(cut(x[y==1],bins))
			j0	<-f0 > 0
			j1	<-f1 > 0
			bins0 <-(bins[-99])[j0]
			bins1 <-(bins[-99])[j1]
			f0	<-f0[j0]
			f1	<-f1[j1]
			maxf <-max(f0,f1)
			f0	<-(0.1*f0)/maxf
			f1	<-(0.1*f1)/maxf
			segments(bins1,-0.05,bins1,f1-0.05, col=col1)
			segments(bins0,-0.05,bins0,-f0-0.05, col=col0)
			lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(-0.05,-0.05))
			text(max(bins0,bins1)+0.02,-0.025,d1lab,cex=0.7)
			text(max(bins0,bins1)+0.02,-0.08,d0lab,cex=0.7)
					}
			}
	stats
}

#################
stat.ES <- function(p,y,cutoff=.7, rounddec1=2, rounddec2=2, rounddec3=2, rounddec4=2,cex.text=1,
      pl=T, g=10, d0lab="", d1lab="", logistic.cal=F, smooth=T, cutline=T, main="", cex.main=1) {
	r <- c(val.prob.yvon(p=p,y=y, pl=pl, g=g, smooth=smooth, logistic.cal=logistic.cal, d0lab=d0lab, d1lab=d1lab,
      xlab="Predicted Probability", ylab="Actual Proportion", 
      statloc=F, legendloc=F)[c(14,15,2,3)],
			sesp(p=p,y=y,cutoff=cutoff)[c(3,4,1,10:13)])
	text(x=rep(0,4),y=c(0.95,0.85,0.75,0.65), c("a|b=1", "slope b", "c stat", "NB"), adj=0, cex=cex.text)
	text(x=rep(0.25,4),y=c(0.95,0.85,0.75,0.65), 
    c(round(r[1],rounddec1),round(r[2],rounddec2),round(r[3],rounddec3),round(r[11],rounddec4)), 
    adj=0, cex=cex.text) 
	if (cutline) arrows(x0=cutoff, y0=.05, x1=cutoff, y1=-0.05, length = 0.08)
	text(x=cutoff,y=.1, "Cutoff", adj=0.5, cex=cex.text)
  title(main, cex.main=cex.main)
names(r)	<-list('Int|slope=1',"Slope", "C (ROC)", "R2 Pearson", 
            "sens", "spec", "acc","NBmodel", "NBtreat", "NBnotreat", "NB")
}
# stat.ES(p=pred1,y=y)
#################

