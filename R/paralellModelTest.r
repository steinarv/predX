
simpleDoRoll <- function(x, f=hw2, n.ahead=10, win.size=500, dates=NULL, ... ){
	n <- length(x)
	if(!is.null(dates))if(length(dates)!=n){
							warning("Length of dates array do not match length of data array")
							dates=NULL
						}
	
	num <- floor((n-win.size)/n.ahead)
  
	dfOut <- data.frame(Convergence = logical(num), InSamp_score = numeric(num), Emp_Level = numeric(num), 
							Exp_Level = numeric(num), MSEin = numeric(num), MSEout = numeric(num), Date = as.Date(NA))
	
  
	a <- 1
	b <- win.size
	i <- 1
	while(b <= n-n.ahead){
    
		# f() standard optmizing function for package predX
		fit <- f(x[a:b], n.ahead=n.ahead, ... )
		if(fit$convergence!=0){
			warning(paste0("Convergence is: ", fit$convergence, 
							" for window: ", a, " to ", b, ". Using paramters from previous optimization."))
			if(a==1)stop("Convergance not possible for first window")
			else fit <- f(x[a:b], n.ahead=n.ahead, param=prevPar, doOptim=FALSE)
			bConv = FALSE
		}else{
			prevPar <- fit$par
			bConv = TRUE
		}
    
		Exp_Level <- fit$vFilOut[n.ahead]
		
		if(!is.null(dates))dfOut[i, "Date"] <- dates[b]		#The date for the last obs in sample
		dfOut[i, c("Convergence", "InSamp_score", "Emp_Level", "Exp_Level", "MSEin", "MSEout")] <- 
								c(	bConv, fit$score, x[b+n.ahead], fit$vFilOut[n.ahead], mean((x[a:b]-fit$vFilIn)^2),
									mean((x[(b+1):(b+n.ahead)]-fit$vFilOut)^2) )
                      
		i <- i+1
		a <- a+n.ahead
		b <- b+n.ahead
	}
  
	dfOut
}


# Takes a series and do rolling window optmiziation and out of sample predictions and evaluation with (possible) several models
parallelTest.byFunc <- function(ser, serName="series1", vf=c(hw2, hw2), vfNames=paste0("model", 1:length(vf)), dates=NULL,
								n.ahead=10, win.size=500, ret.type="level", freq="freq", 
								scoreFunc=paste0("scoreFunc", 1:length(vf)), ... ){
	n <- length(vf)
	if(length(vfNames)!=n){
		warning("Number of function names do not match the number of functions to be evaluated")
		vfNames=rep("", n)
	}
	
	if(!is.null(dates))if(length(dates)!=length(ser)){
							warning("Length of dates array do not match length of data array")
							dates=NULL
						}else{
							daterange <- format(range(dates), "%Y%m%d")
						}
	else{
		daterange <- rep("00000000" ,2)
	}
	if(any(is.na(daterange)))daterange <- rep("00000000" ,2)
	
	
	aggf = switch(ret.type, "level"=function(x){x}, "arithmetic"=aggprod, "logarithmic"=aggsum) #Need to work with this to allow for predicting returns
	
	ptm <- proc.time()
	cl <- makeCluster(rep("localhost", n), type = "SOCK") #n should not exceed number of cores
	registerDoSNOW(cl)

	lOut <- foreach(i=1:n, .packages="predX", .inorder=TRUE) %dopar% {  
		
		f <- vf[[i]]
		
		dfX <- simpleDoRoll(ser, f=f, n.ahead=n.ahead, win.size=win.size, dates, ...)

		medAPE <- median(abs((dfX$Emp_Level[-1]-dfX$Exp_Level[-1])/dfX$Emp_Level[-1]))
		AMAPE <- mean(abs((dfX$Emp_Level[-1]-dfX$Exp_Level[-1])/(dfX$Emp_Level[-1]+dfX$Exp_Level[-1])))
		
		# Identification string for this analysis. 
		rObjName <- paste0( serName, "_", vfNames[i], "_",
							paste(daterange, collapse=""), "_",
							n.ahead, "_", win.size, "_", freq, "_",
							ret.type)
  
		list(series=serName, model=vfNames[i], rObj=rObjName,
			num.pred= nrow(dfX), n.non.conv=sum(!dfX[ ,"Convergence"]), 
			MAE=mean(abs(dfX[, "Emp_Level"]-dfX[, "Exp_Level"])), medAPE=medAPE, AMAPE=AMAPE, 
			daterange=daterange, win.size=win.size, ret.type=ret.type, n.ahead=n.ahead, 
			tot.obs=length(ser), freq=freq, scoreFunc=scoreFunc[i], allData=dfX) 
	}#foreach

	stopCluster(cl)
	print(proc.time() - ptm)
	
	names(lOut) <- vfNames
	lOut
}
