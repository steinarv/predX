# -------------------------------- Basic variance models -------------------------------------------
stdNormVar <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n), HtOut=rep(h, n.ahead), convergence=0)
}

#optTdf(x, dfrange=c(1,00))sum(log(dt(x, df=df)))
#stdTVar <- function(x, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
#	n <- length(x)
#	h <- var(x)
#	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
#}

stdHS <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}

stdEWMA <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}

filterHS <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}
# --------------------------   Estimates simple GARCHish models -------------------------------------

#OPTstdGARCH <- function(x, h0, param){
#	.Call("LLstdGARCH", x, h0, param, PACKAGE = "predX" )
#}

#OPTeGARCH <- function(x, h0, param){
	#.Call("LLeGARCH", x, h0, param, PACKAGE = "predX" )
#}

#testEgarch <- function(x){
#	.Call("HTeGARCH", x, var(x), c( log(var(x))*0.25-0.016, 0.02, 0, 0.75 ), 1, PACKAGE = "predX" )
#}

INVunityf <- function(x){
	x <- sapply(x, FUN=function(z)min(z,0.9999)) #Ensure numeric value in return
	log(x/(1-x))
	}

fMSE <- function(y, x, trim=0)mean((y-x)^2, trim=trim)
fABS <- function(y, x, trim=0)mean(abs(y-x), trim=trim)

OPTexsm1 <- function(y, param, startVal, scorefunc=fMSE){
	n <- length(y)
	fitval <- .Call("EXPSMOOTH1", Y=y, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(y[10:n], fitval[10:n])
}

exsm1 <- function(y, param=NULL, doOptim=TRUE, nout=1, scorefunc=fMSE, solver.method="Brent", solver.control=list()){
	
	n <- length(y)
	if(is.null(param)){param <- 0.5}else{param <- param}
	if(doOptim){
		opt <- optim(param, OPTexsm1, y=y, scorefunc=scorefunc, lower=0, upper=1, 
						startVal=mean(y), method=solver.method, control=solver.control)
		param <- opt$par
		names(opt$par) <- "lambda"
	}
	
	fit <- .Call("EXPSMOOTH1", Y=y, PARAM=param, STARTVAL=mean(y), NOUT=nout, PACKAGE = "predX")

	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}

OPTexsm2 <- function(y, param, thold, startVal, scorefunc=fMSE){
	n <- length(y)
	# Order of arguments given to .Call function does matter!!
	fitval <- .Call("EXPSMOOTH2", Y=y, LAMBA=param, THOLD=thold, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(y[2:n], fitval[1:(n-1)])
}

exsm2 <- function(y, param=NULL, thold=NULL, doOptim=TRUE, nout=1, scorefunc=fMSE, solver.method="Brent", solver.control=list()){
	
	n <- length(y)
	if(is.null(param)){param <- 0.5}else{param <- param}
	if(is.null(thold)){thold <- rep(sd(y)*3,n)}
	if(doOptim){
		opt <- optim(param, OPTexsm2, y=y, startVal=mean(y), scorefunc=scorefunc, 
						thold=thold, lower=0, upper=1, 
						method=solver.method, control=solver.control)
		param <- opt$par
	}
	
	fit <- .Call("EXPSMOOTH2", Y=y, PARAM=param, THOLD=thold, STARTVAL=mean(y), NOUT=nout, PACKAGE = "predX")

	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}


OPTexsm3 <- function(y, param, startVal, scorefunc){
	n <- length(y)
	fitval <- .Call("EXPSMOOTH3", Y=y, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX" )
	scorefunc(y[10:n], fitval[10:n])
}

exsm3 <- function(y, param=NULL, doOptim=TRUE, nout=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(y)
	if(is.null(param)){param <-  INVunityf(c(0.5,0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTexsm3, y=y, startVal=c(mean(y), 0, y[1]), scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("EXPSMOOTH3", Y=y, PARAM=param, STARTVAL=c(mean(y), 0, y[1]),
					NOUT=nout, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}

# Seasonal smoothing without trend, with possible esplanatory variables
OPTseasexsm <- function(y, s, param, startVal, scorefunc){
	n <- length(y)
	fitval <- .Call("SEASEXPSMOOTH", Y=y, S=s, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(y[(s*2):n], fitval[(s*2):n])
}

seasexsm <- function(y, s, param=NULL, doOptim=TRUE, nout=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(y)
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(y[1:10]); startVal[2:(s+1)] <- aggregate(y, by=list(rep_len(1:s,n)), mean)$x/mean(y)
	startVal[s+2] <- startVal[2]
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTseasexsm, y=y, s=s, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("SEASEXPSMOOTH", Y=y, S=s, PARAM=param, STARTVAL=startVal,
					NOUT=nout, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}

OPTseasregexsm <- function(y, s, X, param, startVal, scorefunc){
	n <- length(y); ind <- (s*2):n
	fitval <- .Call("SEASEXPSMOOTH", Y=y, S=s, PARAM=param, STARTVAL=startVal, NOUT=0, PACKAGE = "predX")
	lmfit <- lm(y[ind]~fitval[ind]+X[ind, ])
	scorefunc(y[ind], lmfit$fitted.values)
}

seasregexsm <- function(y, s, X, param=NULL, doOptim=TRUE, nout=0, Xout=NULL, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(y)
	if(nrow(X)!=n)stop("Provide a matrix with explanatory variables and number of rows equal to length y")
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(y[1:10]); startVal[2:(s+1)] <- aggregate(y, by=list(rep_len(1:s,n)), mean)$x/mean(y)
	startVal[s+2] <- startVal[2]
	
	if(!is.null(Xout))
	if(nrow(Xout)!=nout || ncol(Xout)!=ncol(X))
		stop("Explanatory variables provided for predictions does not match")
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTseasregexsm, y=y, s=s, X=X, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	
	
	ind <- (s*2):n
	filterfit <- .Call("SEASEXPSMOOTH", Y=y, S=s, PARAM=param, STARTVAL=startVal,
						NOUT=nout, PACKAGE = "predX" )
	lmfit <- lm(y[ind]~filterfit[ind]+X[ind,])
	
	fit <- cbind(1, filterfit, rbind(X,Xout))%*%matrix(coef(lmfit), ncol=1)
	
	
	
	lOut <- list(fitIn=fit[1:n, 1], fitOut=fit[ifelse(nout==0,n,n+1):(n+nout), 1], lmfit=lmfit)
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}


# Similar day smoothing without trend, with possible esplanatory variables
OPTsimdexsm <- function(y, days, s, param, startVal, scorefunc){
	n <- length(y)
	fitval <- .Call("SIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, STARTVAL=startVal, PACKAGE = "predX")
	scorefunc(y[(s*2):n], fitval[(s*2):n])
}

simdexsm <- function(y, days, param=NULL, doOptim=TRUE, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(y); nout <- length(days)-n; s <- length(unique(days))
	startVal = rep(NA, s+1) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(y[1:10])
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	startVal[2:(s+1)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
		
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTsimdexsm, y=y, days=days[1:n], s=s, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("SIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, STARTVAL=startVal,
					PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}


# RiskMetric adjusted similar day smoothing without trend, with possible esplanatory variables
OPTrmsimdexsm <- function(y, days, s, thold, param, startVal, scorefunc){
	n <- length(y)

	fitval <- .Call("RMSIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, 
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX")
	scorefunc(y[(s*2):n], fitval[(s*2):n])

}

rmsimdexsm <- function(y, days, param=NULL, doOptim=TRUE, thold=2, 
			solver.method="Nelder-Mead", solver.control=list()){
			
	n <- length(y); nout <- length(days)-n; s <- length(unique(days))
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- var(y); startVal[2] <- mean(y[1:10]);
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	startVal[3:(s+2)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))
	}else{param <- INVunityf(param)}
	
	if(doOptim){
	if(is.null(thold)){
		
		opt <- list(value=Inf)
		for(i in seq(from=2, to=4, by=0.2)){
		newopt <- optim(param, OPTrmsimdexsm, y=y, days=days[1:n], s=s, startVal=startVal, scorefunc=fMSE,
				thold=i, method=solver.method, control=solver.control)
		if(newopt$value < opt$value){opt <- newopt; thold <- i}
		}
	
	}else{
		opt <- optim(param, OPTrmsimdexsm, y=y, days=days[1:n], s=s, startVal=startVal, scorefunc=fMSE,
				thold=thold, method=solver.method, control=solver.control)
	}
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("RMSIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, 
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt, list(thold=thold))
	
	lOut

}

###############################################################################################
# RiskMetric adjusted similar day smoothing without trend, with possible esplanatory variables#
###############################################################################################

#For regression with intercept include a columns of ones in X

OPTrmsimdregexsm <- function(y, X, days, s, thold, param, startVal, scorefunc){
	n <- length(y)
	ind <- (s*2):n
	
	fitval <- .Call("RMSIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, 
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX")
	
	
	lmfit <- lm(y[ind]~fitval[ind]+X[ind, ]-1)
	scorefunc(y[ind], lmfit$fitted.values)
}

rmsimdregexsm <- function(y, X, days, param=NULL, doOptim=TRUE, thold=2, 
			solver.method="Nelder-Mead", solver.control=list()){
			
	n <- length(y); nout <- length(days)-n; s <- length(unique(days))
	if(nrow(X)!=(n+nout))stop("Provide a matrix with explanatory variables and number of rows 
				equal to length days")
	
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- var(y); startVal[2] <- mean(y[1:10]);
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	startVal[3:(s+2)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))
	}else{param <- INVunityf(param)}
	
	if(doOptim){
	if(is.null(thold)){
		
		opt <- list(value=Inf)
		for(i in seq(from=2, to=4, by=0.2)){
		newopt <- optim(param, OPTrmsimdregexsm, y=y, X=X[1:n, ,drop=FALSE], days=days[1:n], s=s, 
				startVal=startVal, scorefunc=fMSE, thold=i, method=solver.method, 
				control=solver.control)
		if(newopt$value < opt$value){opt <- newopt; thold <- i}
		}
	
	}else{
		opt <- optim(param, OPTrmsimdregexsm, y=y, X=X[1:n, ,drop=FALSE], days=days[1:n], s=s, 
				startVal=startVal, scorefunc=fMSE, thold=thold, method=solver.method,
				control=solver.control)
	}
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	
	
	ind <- (s*2):n
	# Filter data and predict nout days ahead
	filterfit <- .Call("RMSIMDAYEXPSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, 
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX" )
	# Linear regression with filtered values and X matrix as explanatory variables
	lmfit <- lm(y[ind]~filterfit[ind]+X[ind, , drop=FALSE]-1)
	# Fitted values
	fit <- cbind(filterfit, X)%*%matrix(coef(lmfit), ncol=1)
	
	lOut <- list(fitIn=fit[1:n, 1], lmfit=lmfit)
	if(nout>0)lOut <- c(lOut, list(fitOut=fit[(n+1):(n+nout), 1]))
	if(doOptim)lOut <- c(lOut, opt, list(thold=thold))
	
	lOut

}


#########################################################################################################
#		Similar day-seasonal exponential smoothing with model  					#
#		optimization over a defined number of forecasts. 					#
#		A riskmetric type of time varying volatilyt of the					#
#		prediction error of the model is estimated. Observations that deviates 			#
# 		from the predicted value with a defined number of standard deviations 			#
#		may be ignored.										#
#########################################################################################################
OPTsimdaysmooth <- function(y, ymat, days, s, opt.nout, thold, param, startVal, scorefunc){
	n <- length(y)

	fitval <- .Call("SIMDAYSMOOTH", Y=y, DAYS=days, S=s, OPTNOUT=opt.nout, PARAM=param,
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX")
			
	scorefunc(ymat[-(1:(s*2)), ], fitval[-(1:(s*2)), ])

}

simdaysmooth <- function(y, days, param=NULL, doOptim=TRUE, thold=2, opt.nout=7,
			solver.method="Nelder-Mead", solver.control=list()){
			
	n <- length(y); nout <- length(days)-n; s <- length(unique(days))
	
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- var(y); startVal[2] <- mean(y[1:10]);
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	startVal[3:(s+2)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))
	}else{param <- INVunityf(param)}
	
	if(doOptim){
		#Matrix used for efficient estimation of model predictions errors at each step in filtration
		ymat <- t(sapply(1:(n-opt.nout+1), FUN=function(x)y[x:(x+opt.nout-1)]))
		
		opt <- optim(param, OPTsimdaysmooth, y=y, ymat, days=days[1:n], s=s, opt.nout=opt.nout, 
				startVal=startVal, scorefunc=fMSE, thold=thold, 
				method=solver.method, control=solver.control)
	
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("SIMDAYSMOOTH", Y=y, DAYS=days, S=s, PARAM=param, 
			THOLD=thold, STARTVAL=startVal, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt, list(thold=thold))
	
	lOut

}






































stdGARCH <- function(x, h0=NULL, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	
	if(is.null(h0))h0 <- var(x);
	if(is.null(param)){param <-  log(c(var(x)*0.05, 0.2, 0.75))}else{param <- log(param)}
	if(doOptim){
		opt <- optim(param, OPTstdGARCH, x=x, h0=h0, method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- exp(opt$par)
		opt$value <- -opt$value
	}
	
	Ht <- 1#.Call("HTstdGARCH", x, h0, param, n.ahead, PACKAGE = "predX" )
	
	
	lOut <- list(HtIn=Ht[1:n], HtOut=Ht[(n+1):(n+n.ahead)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
	
}

eGARCH <- function(x, h0=NULL, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	
	if(is.null(h0))h0 <- var(x);
	if(is.null(param)) param <- c( log(var(x))*0.25-0.016, 0.02, 0, 0.75 ) #a0=log(var(x))*(1-b0)-a1*sqrt(2/pi)
	if(doOptim){
		opt <- optim(param, OPTeGARCH, x=x, h0=h0, method=solver.method, control=solver.control)
		param <- opt$par
		opt$value <- -opt$value
	}

	Ht <- 1#exp(.Call("HTeGARCH", x, h0, param, n.ahead, PACKAGE = "predX" ))
	
	
	lOut <- list(HtIn=Ht[1:n], HtOut=Ht[(n+1):(n+n.ahead)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
	
}
