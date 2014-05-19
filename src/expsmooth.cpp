#include "expsmooth.h"
#include "miscf.h"

using namespace Rcpp ;


// ---------------------------- Holt-Winters triple exponential smoothing (multiplicative) ----------------------------------
SEXP HW_TRIPLE(SEXP Y, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP STARTVAL, SEXP NOUT, SEXP MULT) {
	NumericVector nvX(Y); int n = nvX.size(); int f = as<int>(NOUT);
	int s = as<int>(S); int o = as<int>(OPTNOUT); int m = as<int>(MULT);
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double beta = nvPARAM(1); double gamma = nvPARAM(2);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	
	double dL = nvSTARTVAL(0);
	double dL1 = dL;
	double dT = nvSTARTVAL(1);
	for(int i=0;i<s;i++)nvS(i)=nvSTARTVAL(i+2);
	
	int j;
	
	for(int i=1;i<(n+f);i++){
	if(m==1){	// If multiplicative	
		
		if(i<n){
			if(i < (n-o)){
				for(j=0; j<o; j++)
				nvFIL(i, j)=(dL+dT*(j+1))*nvS((i+j)%s); //Predicted/Filtered value for "today+j"	
			}else{
				nvFIL(i, 0)=(dL+dT)*nvS(i%s); //Predicted/Filtered value for "today"	
			}
			
			
			dL1=dL;
			dL=alfa*nvX(i)/nvS(i%s)+(1-alfa)*(dL+dT); 	//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(i%s)=gamma*nvX(i)/dL+(1-gamma)*nvS(i%s); 	//Season updated
		}else{
			nvFIL(i, 0) = (dL+dT)*nvS(i%s);
		}
	
	}else{ 		// If additive
	
		if(i<n){
			if(i < (n-o)){
				for(j=0; j<o; j++)
				nvFIL(i, j)=(dL+dT*(j+1))+nvS((i+j)%s); //Predicted/Filtered value for "today+j"	
			}else{
				nvFIL(i, 0)=(dL+dT)+nvS(i%s); //Predicted/Filtered value for "today"	
			}
			
			
			dL1=dL;
			dL=alfa*(nvX(i)-nvS(i%s))+(1-alfa)*(dL+dT); 	//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(i%s)=gamma*(nvX(i)-dL)+(1-gamma)*nvS(i%s); 	//Season updated
		}else{
			nvFIL(i, 0) = (dL+dT)+nvS(i%s);
		}
	
	
	}
	}
	
	return(wrap(nvFIL));
}

