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
#													#
#########################################################################################################

OPThw_simday <- function(y, ymat, days, s, opt.nout, param, trend, thold, startVal, scorefunc, trim, mult){
	n <- length(y)
	
	if(trend){
	  param_ <- param
	}else{
	  param_ <- c(param[1], -1000, param[2])
	}
	

	fitval <- .Call("HW_SIMDAY", Y=y, DAYS=days, S=s, OPTNOUT=opt.nout, PARAM=param_, THOLD=thold,
			             STARTVAL=startVal, MULT=mult, PACKAGE = "predX")
			
	scorefunc(ymat[(s*2+1):(n-opt.nout+1), ], fitval[(s*2+1):(n-opt.nout+1), ], trim=trim)

}

hw_simday <- function(y, days, param=NULL, doOptim=TRUE, opt.nout=7, trend=TRUE, thold=3, 
			mult=FALSE, scorefunc=fMSE, trim=0, solver.method="Nelder-Mead", solver.control=list()){
			
	n <- length(y); nout <- length(days)-n; s <- length(unique(days));
	nn <- min(10*s, n) #Number of days used of initializing seasonal component
	
	
	startVal = rep(NA, s+3) #Level0 Trend0, and Seas1:s
	startVal[1] <- sd(y); startVal[2] <- mean(y[1:nn]); startVal[3] <- 0
	
	if(mult){
		startVal[4:(s+3)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x/mean(y[1:nn])
	}else{
		startVal[4:(s+3)] <- aggregate(y[1:nn], by=list(days[1:nn]), mean)$x-mean(y[1:nn])
	}

	
	nparam <- (2+trend)
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
		

		opt <- optim(param, OPThw_simday, y=y, ymat=ymat, days=days, s=s, opt.nout=opt.nout, 
		    	trend=trend, thold=thold, startVal=startVal, scorefunc=scorefunc, trim=trim, 
			mult=mult, method=solver.method, control=solver.control)

		
		param <- opt$par
		if(trend){
	  		param_ <- param
		}else{
	  		param_ <- c(param[1], -1000, param[2])
	  	}
		
		opt$par <- 1/(1+exp(-opt$par))
	}
	
	fit <- .Call("HW_SIMDAY", Y=y, DAYS=days, S=s, OPTNOUT=1, PARAM=param_, THOLD=thold, 	
			        STARTVAL=startVal, MULT=mult, PACKAGE = "predX" )
	
	lOut <- list(startVal=startVal, fitIn=fit[1:n, 1])
	if(nout>0)lOut <- c(lOut, list(fitOut=fit[(n+1):(n+nout), 1]))
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut

}
