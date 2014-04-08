
doRoll <- function(x, f=hw2, n.ahead=10, win.size=500, dates=NULL, ... ){
	n <- length(x)
	
	if(is.null(dates) || length(dates)!=n){
							warning("Length of dates array do not match length of data array")
							dates <- as.Date(1:n, origin="1900-01-01")
	}
	
	num <- floor((n-win.size)/n.ahead)
  
	a <- 1
	b <- win.size
	i <- 1
	while(b <= (n-n.ahead)){
    
		# f() standard optmizing function for package predX
		fit <- f(x[a:b], nout=n.ahead, ... )
		
		#If first time, define output data frame
		if(b==win.size){
			df1 <- data.frame(Date = as.Date(NA), Convergence = logical(num), InSamp_score = numeric(num),
								MSEin = numeric(num))
			mat1 <- matrix(NA, ncol=length(fit$par), nrow=num) #Parameters of fitted model
			colnames(mat1) <- names(fit$par)
			mat2 <- matrix(NA, ncol=n.ahead*2, nrow=num) #Empirical and predicted values 1 to n.ahead
			colnames(mat2) <- c(paste0("Emp_t+", 1:n.ahead), paste0("Pred_t+", 1:n.ahead))
			#dfOut <- cbind(df1, mat1, mat2)
			#Combine after loop finished
										
		}
		
		
		if(fit$convergence!=0){
			warning(paste0("Convergence is: ", fit$convergence, 
							" for window: ", dates[a], " to ", dates[b], ". Using paramters from previous optimization."))
			if(a==1)stop("Convergance not possible for first window")
			else fit <- f(x[a:b], nout=n.ahead, param=prevPar, doOptim=FALSE)
			bConv = FALSE
			mat1[i, ] <- prevPar
		}else{
			prevPar <- fit$par
			bConv = TRUE
			mat1[i, ] <- fit$par
		}
		
    
		df1[i, "Date"] <- dates[b]
		df1[i, c("Convergence", "InSamp_score", "MSEin")] <- c(bConv, fit$value, mean((x[a:b]-fit$fitIn)^2))
		
		mat2[i, ] <- c(x[b+(1:n.ahead)], fit$fitOut)
		                      
		i <- i+1
		a <- a+n.ahead
		b <- b+n.ahead
	}
  
	cbind(df1, mat1, mat2)
}


# Takes a series and do rolling window optmiziation and out of sample predictions and evaluation with (possible) several models
parDoRoll <- function(ser, serName="series1", vf=c(hw2, hw2), vfNames=paste0("model", 1:length(vf)), dates=NULL,
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
