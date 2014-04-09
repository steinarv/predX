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

INVunityf <- function(x)log(x/(1-x))

fMSE <- function(y,x)mean((y-x)^2)

OPTexsm1 <- function(x, param, startVal, scorefunc=fMSE){
	n <- length(x)
	fitval <- .Call("EXPSMOOTH1", X=x, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(x[10:n], fitval[10:n])
}

exsm1 <- function(x, param=NULL, doOptim=TRUE, nout=1, scorefunc=fMSE, solver.method="Brent", solver.control=list()){
	
	n <- length(x)
	if(is.null(param)){param <- 0.5}else{param <- param}
	if(doOptim){
		opt <- optim(param, OPTexsm1, x=x, scorefunc=scorefunc, lower=0, upper=1, 
						startVal=mean(x), method=solver.method, control=solver.control)
		param <- opt$par
		names(opt$par) <- "lambda"
	}
	
	fit <- .Call("EXPSMOOTH1", X=x, PARAM=param, STARTVAL=mean(x), NOUT=nout, PACKAGE = "predX")

	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}

OPTexsm2 <- function(x, param, thold, startVal, scorefunc=fMSE){
	n <- length(x)
	# Order of arguments given to .Call function does matter!!
	fitval <- .Call("EXPSMOOTH2", X=x, LAMBA=param, THOLD=thold, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(x[2:n], fitval[1:(n-1)])
}

exsm2 <- function(x, param=NULL, thold=NULL, doOptim=TRUE, nout=1, scorefunc=fMSE, solver.method="Brent", solver.control=list()){
	
	n <- length(x)
	if(is.null(param)){param <- 0.5}else{param <- param}
	if(is.null(thold)){thold <- rep(sd(x)*3,n)}
	if(doOptim){
		opt <- optim(param, OPTexsm2, x=x, startVal=mean(x), scorefunc=scorefunc, 
						thold=thold, lower=0, upper=1, 
						method=solver.method, control=solver.control)
		param <- opt$par
	}
	
	fit <- .Call("EXPSMOOTH2", X=x, PARAM=param, THOLD=thold, STARTVAL=mean(x), NOUT=nout, PACKAGE = "predX")

	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}


OPTexsm3 <- function(x, param, startVal, scorefunc){
	n <- length(x)
	fitval <- .Call("EXPSMOOTH3", X=x, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX" )
	scorefunc(x[10:n], fitval[10:n])
}

exsm3 <- function(x, param=NULL, doOptim=TRUE, nout=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	if(is.null(param)){param <-  INVunityf(c(0.5,0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTexsm3, x=x, startVal=c(mean(x), 0, x[1]), scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("EXPSMOOTH3", X=x, PARAM=param, STARTVAL=c(mean(x), 0, x[1]),
					NOUT=nout, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
}

# Seasonal smoothing without trend, with possible esplanatory variables
OPTseasexsm <- function(x, s, param, startVal, scorefunc){
	n <- length(x)
	fitval <- .Call("SEASEXPSMOOTH", X=x, S=s, PARAM=param, STARTVAL=startVal, NOUT=1, PACKAGE = "predX")
	scorefunc(x[(s*2):n], fitval[(s*2):n])
}

seasexsm <- function(x, s, param=NULL, doOptim=TRUE, nout=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(x[1:10]); startVal[2:(s+1)] <- aggregate(x, by=list(rep_len(1:s,n)), mean)$x/mean(x)
	startVal[s+2] <- startVal[2]
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTseasexsm, x=x, s=s, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("SEASEXPSMOOTH", X=x, S=s, PARAM=param, STARTVAL=startVal,
					NOUT=nout, PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}

OPTseasregexsm <- function(x, s, X, param, startVal, scorefunc){
	n <- length(x); ind <- (s*2):n
	fitval <- .Call("SEASEXPSMOOTH", X=x, S=s, PARAM=param, STARTVAL=startVal, NOUT=0, PACKAGE = "predX")
	lmfit <- lm(x[ind]~fitval[ind]+X[ind, ])
	scorefunc(x[ind], lmfit$fitted.values)
}

seasregexsm <- function(x, s, X, param=NULL, doOptim=TRUE, nout=0, Xout=NULL, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	if(nrow(X)!=n)stop("Provide a matrix with explanatory variables and number of rows equal to length x")
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(x[1:10]); startVal[2:(s+1)] <- aggregate(x, by=list(rep_len(1:s,n)), mean)$x/mean(x)
	startVal[s+2] <- startVal[2]
	
	if(!is.null(Xout))
	if(nrow(Xout)!=nout || ncol(Xout)!=ncol(X))
		stop("Explanatory variables provided for predictions does not match")
	
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTseasregexsm, x=x, s=s, X=X, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	
	
	ind <- (s*2):n
	filterfit <- .Call("SEASEXPSMOOTH", X=x, S=s, PARAM=param, STARTVAL=startVal,
						NOUT=nout, PACKAGE = "predX" )
	lmfit <- lm(x[ind]~filterfit[ind]+X[ind,])
	
	fit <- cbind(1, filterfit, rbind(X,Xout))%*%matrix(coef(lmfit), ncol=1)
	
	
	
	lOut <- list(fitIn=fit[1:n, 1], fitOut=fit[ifelse(nout==0,n,n+1):(n+nout), 1], lmfit=lmfit)
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}


# Similar day smoothing without trend, with possible esplanatory variables
OPTsimdexsm <- function(x, days, s, param, startVal, scorefunc){
	n <- length(x)
	fitval <- .Call("SIMDAYEXPSMOOTH", X=x, DAYS=days, S=s, PARAM=param, STARTVAL=startVal, PACKAGE = "predX")
	scorefunc(x[(s*2):n], fitval[(s*2):n])
}

simdexsm <- function(x, days, param=NULL, doOptim=TRUE, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x); nout <- length(days)-n; s <- length(unique(days))
	startVal = rep(NA, s+1) #Level0 and Seas1:(s+1)
	startVal[1] <- mean(x[1:10])
	nn <- min(5*s, n) #Number of days used of initializing seasonal component
	startVal[2:(s+1)] <- aggregate(x[1:nn], by=list(days[1:nn]), mean)$x/mean(x)
		
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5))}else{param <- INVunityf(param)}
	if(doOptim){
		opt <- optim(param, OPTsimdexsm, x=x, days=days[1:n], s=s, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("SIMDAYEXPSMOOTH", X=x, DAYS=days, S=s, PARAM=param, STARTVAL=startVal,
					PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}


# RiskMetric adjusted similar day smoothing without trend, with possible esplanatory variables
OPTrmsimdexsm <- function(x, days, s, param, startVal, scorefunc){
	n <- length(x)
	fitval <- .Call("RMSIMDAYEXPSMOOTH", X=x, DAYS=days, S=s, PARAM=param, STARTVAL=startVal, PACKAGE = "predX")
	scorefunc(x[(s*2):n], fitval[(s*2):n])
}

rmsimdexsm <- function(x, days, param=NULL, doOptim=TRUE, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x); nout <- length(days)-n; s <- length(unique(days))
	startVal = rep(NA, s+2) #Level0 and Seas1:(s+1)
	startVal[1] <- var(x); startVal[2] <- mean(x[1:10]);
	nn <- min(5*s, n) #Number of days used of initializing seasonal component
	startVal[3:(s+2)] <- aggregate(x[1:nn], by=list(days[1:nn]), mean)$x/mean(x)
		
	if(is.null(param)){param <-  INVunityf(c(0.5, 0.5, 0.5))
		}else{param <- INVunityf(param/c(1,1,4))}
	
	if(doOptim){
		opt <- optim(param, OPTrmsimdexsm, x=x, days=days[1:n], s=s, startVal=startVal, scorefunc=fMSE,
						method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- 1/(1+exp(-opt$par))*c(1,1,4)
	}
	
	fit <- .Call("RMSIMDAYEXPSMOOTH", X=x, DAYS=days, S=s, PARAM=param, STARTVAL=startVal,
					PACKAGE = "predX" )
	
	lOut <- list(fitIn=fit[1:n], fitOut=fit[(n+1):(n+nout)])
	if(doOptim)lOut <- c(lOut, opt)
	
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
