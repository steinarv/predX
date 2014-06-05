#INVunityf <- function(x){
#x <- sapply(x, FUN=function(z)min(z,0.9999)) #Ensure numeric value in return
#log(x/(1-x))
#}

#fMSE <- function(y, x, trim=0)mean((y-x)^2, trim=trim)
#fABS <- function(y, x, trim=0)mean(abs(y-x), trim=trim)

OPThw_triple <- function(y, ymat, s, opt.nout, param, trend, seas, startVal, scorefunc, trim, mult){
	n <- length(y)
	
	if(trend & seas){
	  param_ <- param
	}else{
	  param_ <- c(param[1], -1000, -1000)
	  if(trend)param_[2] <- param[2]
	  if(seas)param_[3] <- param[2]
  }
	

	fitval <- .Call("HW_TRIPLE", Y=y, S=s, OPTNOUT=opt.nout, PARAM=param_,
			             STARTVAL=startVal, NOUT=0, MULT=mult, PACKAGE = "predX")
			
	scorefunc(ymat[(s*2+1):(n-opt.nout+1), ], fitval[(s*2+1):(n-opt.nout+1), ], trim=trim)

}

hw_triple <- function(y, s, nout=0, param=NULL, doOptim=TRUE, opt.nout=7, trend=TRUE, seas=TRUE, 
			mult=FALSE, scorefunc=fMSE, trim=0, solver.method="Nelder-Mead", solver.control=list()){
			
	n <- length(y)
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	
	startVal = rep(NA, s+2) #Level0 Trend0, and Seas1:s
	startVal[1] <- mean(y[1:nn]); startVal[2] <- 0
	
	if(seas){
		if(mult){
		startVal[3:(s+2)] <- aggregate(y[1:nn], by=list(rep(1:s, len=nn)), mean)$x/mean(y[1:nn])
		}else{
		startVal[3:(s+2)] <- aggregate(y[1:nn], by=list(rep(1:s, len=nn)), mean)$x-mean(y[1:nn])
		}
	}else{
		startVal[3:(s+2)] <- 1
	}
	
	nparam <- (1+seas+trend)
	if(is.null(param) || length(param)!=nparam){ 
	   param <-  INVunityf(rep(0.25, nparam))
	}else{param <- INVunityf(param)}
	

	if(doOptim){
		#Matrix used for efficient estimation of model predictions errors at each step in filtration
		if(opt.nout>1){
			ymat <- t(sapply(1:(n-opt.nout+1), FUN=function(x)y[x:(x+opt.nout-1)]))
		}else{
			ymat <- matrix(y, ncol=1)
		}
		
		if(seas || trend){
			opt <- optim(param, OPThw_triple, y=y, ymat=ymat, s=s, opt.nout=opt.nout, 
			    	trend=trend, seas=seas, startVal=startVal, scorefunc=scorefunc, trim=trim, 
				mult=mult, method=solver.method, control=solver.control)
		}else{
			opt <- optim(param, OPThw_triple, y=y, ymat=ymat, s=s, opt.nout=opt.nout, 
			    	trend=trend, seas=seas, startVal=startVal, scorefunc=scorefunc, trim=trim, 
				mult=mult, lower=-10, upper=10, method="Brent", control=solver.control)
		}
		
		
		param <- opt$par
		if(trend & seas){
	  		param_ <- param
		}else{
	  		param_ <- c(param[1], -1000, -1000)
	  		if(trend)param_[2] <- param[2]
	  		if(seas)param_[3] <- param[2]
		}
		
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("HW_TRIPLE", Y=y, S=s, OPTNOUT=1, PARAM=param_, 
			        STARTVAL=startVal, NOUT=nout, MULT=mult, PACKAGE = "predX" )
	
	lOut <- list(startVal=startVal, fitIn=fit[1:n, 1])
	if(nout>0)lOut <- c(lOut, list(fitOut=fit[(n+1):(n+nout), 1]))
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}






#########################################################################################################
#          					                                                        #
#		A more advanced exponential smoothing for times series with				#
#		with irregularity in the seasonal component. This model also				#
#		tracks the prediction error which makes it possible to attempt				#
#		identifying outliers based on the size of the prediction error.				#
#		Finaly this model also allows for passing external guesses 				#
#		for the level and weights to be assigned to it.						#
#													#
#########################################################################################################

OPThw_simday <- function(y, ymat, days, l=l, s, opt.nout, param, trend, w1, w2, setw, thold, startVal, scorefunc, trim, mult){
	n <- length(y)
	
	if(trend & setw){
  		param_ <- param
	}else if(trend){
  		param_ <- c(param[1], param[2], param[3], w1, w2)
  	}else if(setw){
  		param_ <- c(param[1], INVunityf(0), param[2], param[3], param[4])
  	}else{
  		param_ <- c(param[1], INVunityf(0), param[2], w1, w2)
  	}	


	fitval <- .Call("HW_SIMDAY", Y=y, DAYS=days, L=l, S=s, OPTNOUT=opt.nout, PARAM=param_, THOLD=thold,
			             STARTVAL=startVal, MULT=mult, PACKAGE = "predX")
			
	scorefunc(ymat[(s*2+1):(n-opt.nout+1), ], fitval[(s*2+1):(n-opt.nout+1), ], trim=trim)

}

hw_simday <- function(y, days, l=NULL, param=NULL, doOptim=TRUE, opt.nout=7, trend=TRUE, thold=3, setw=FALSE,
			mult=FALSE, scorefunc=fMSE, trim=0, solver.method="Nelder-Mead", solver.control=list()){
	
	#setw needs to be TRUE if w1 and w2 is either provided through param or to be optimized
	#if Trend is True a parameter for the trend coefficient must be provided
	
	n <- length(y); nout <- length(days)-n; s <- length(unique(days));
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	
	nparam <- 2+trend+setw*2
	# Parameter vector, is transformed trough 1/(1+exp(-x)) in c++ to ensure 0<>1
	if(length(param)<nparam & doOptim){	#Bad param vector for optimization 
		param <-  rep(0.25, nparam) #Create start values
	}else if(length(param)<nparam){ #Not enough parameters and no optimization planed.
		stop("Not enough parameters to run model")	
	}
	
	param <- INVunityf(param)
	
	# Level vector and its weights
	if(is.null(l)){
		l <- numeric(n+nout) #Pass null vector
		w1 <- INVunityf(1); w2 <- INVunityf(0); #Lock weights
		setw=FALSE
	}else if(!setw & length(param)<(nparam+2)){
		w1 <- w2 <- INVunityf(0.5)
	}else{
		w1 <- param[length(param)-1]; w2 <- param[length(param)]
	}
	
	#print(paste0("w1 = ", w1, ", w2 = ", w2))
	
	# Start values
	startVal = rep(NA, s+3) #Level0 Trend0, and Seas1:s
	startVal[1] <- sd(y); startVal[2] <- mean(y[1:nn]); startVal[3] <- 0
	
	if(mult){
		startVal[4:(s+3)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
	}else{
		startVal[4:(s+3)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x-mean(y[1:nn])
	}

	
	# Optimization
	if(doOptim){
		#Matrix used for efficient estimation of model predictions errors at each step in filtration
		if(opt.nout>1){
			ymat <- t(sapply(1:(n-opt.nout+1), FUN=function(x)y[x:(x+opt.nout-1)]))
		}else{
			ymat <- matrix(y, ncol=1)
		}
		

		opt <- optim(param[1:nparam], OPThw_simday, y=y, ymat=ymat, days=days[1:n], l=l[1:n], s=s, opt.nout=opt.nout, 
			setw=setw, trend=trend, w1=w1, w2=w2, thold=thold, startVal=startVal, scorefunc=scorefunc, 
			trim=trim, mult=mult, method=solver.method, control=solver.control)

		
		param <- opt$par
		
	}else{
		opt <- list(value=NA, par=numeric(5))
	}
	
	
	if(trend & setw){
		param_ <- param
	}else if(trend){
		param_ <- c(param[1], param[2], param[3], w1, w2)
	}else if(setw){
		param_ <- c(param[1], INVunityf(0), param[2], param[3], param[4])
	}else{
		param_ <- c(param[1], INVunityf(0), param[2], w1, w2)
	}
			
	
	

	#print(paste0("w1 = ", w1, ", w2 = ", w2))
	#print(param_)
	fit <- .Call("HW_SIMDAY", Y=y, DAYS=days, L=l, S=s, OPTNOUT=1, PARAM=param_, THOLD=thold, 	
			        STARTVAL=startVal, MULT=mult, PACKAGE = "predX" )
	#Parameters is passed by address and param_ is altered (1/(1+exp(-x)))

  	opt$par <- param_
  	names(opt$par) <- c("alpha", "beta", "gamma", "w1", "w2")
	
	lOut <- c(opt, list(startVal=startVal, fitIn=fit[1:n, 1]))
	if(nout>0)lOut <- c(lOut, list(fitOut=fit[(n+1):(n+nout), 1]))	
	if(!doOptim)lOut$value <-
	scorefunc(y[(s*2+1):(n-opt.nout+1)], lOut$fitIn[(s*2+1):(n-opt.nout+1)], trim=trim)
	
	lOut

}
