# ------------------------- Support functions --------------------------------
rollavg <- function(x, n){
	.Call("rollavg", x, n, PACKAGE = "predX" )
}

pow2 <- function(x)x^2

# -----------------------------------------------------------------------------

# ------- Functions for aggregating data across time --------
aggsum <- function(x)sum(x)
aggprod <- function(x)prod(1+x)-1
# -----------------------------------------------------------------------------


#------------ Inverse unity function ------------------------------------------
INVunityf <- function(x){
	x <- sapply(x, FUN=function(z)min(z,0.9999)) #Ensure numeric value in return
	max(-100000000, log(x/(1-x)))
	}
# .............................................................................

# ----------- Score functiona -------------------------------------------------
fMSE <- function(y, x, trim=0)mean((y-x)^2, trim=trim)
fABS <- function(y, x, trim=0)mean(abs(y-x), trim=trim)
# ------------------------------------------------------------------------------
