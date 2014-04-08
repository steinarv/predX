# ------------------------- Support functions --------------------------------
rollavg <- function(x, n){
	.Call("rollavg", x, n, PACKAGE = "predX" )
}

pow2 <- function(x)x^2

# -----------------------------------------------------------------------------

# ------- Functions for aggregating data across time --------
aggsum <- function(x)sum(x)
aggprod <- function(x)prod(1+x)-1
# -----------------------------------------------------------