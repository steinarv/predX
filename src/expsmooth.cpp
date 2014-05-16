#include "expsmooth.h"

using namespace Rcpp ;

// Function that ensures smoothing parameters between 0 and 1
void unityFunc(NumericVector &y){
  for(int i=0;i<y.size();i++)
    y(i) = 1/(1+::exp(-y(i)));

}


// ---------------------------- Holt-Winters double exponential smoothing -------------------------------
SEXP HW_1(SEXP Y, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) {
	NumericVector nvX(Y); int n = nvX.size(); int f = as<int>(NOUT);
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double beta = nvPARAM(1);
		  
	NumericVector nvL(n); NumericVector nvT(n); NumericVector nvFIL(n+f);
	NumericVector nvSTARTVAL(STARTVAL);
	nvL(0) = nvSTARTVAL(0); nvT(0) = nvSTARTVAL(1); nvFIL(0) = nvX(0);
	  
	for(int i=1;i<(n+f);i++){
		if(i<n){
			nvFIL(i)=nvL(i-1)+nvT(i-1); //Predicted/Filtered value for "today"
			nvL(i)=alfa*nvX(i)+(1-alfa)*nvL(i); //Level updated with value of today
			nvT(i)=beta*(nvL(i)-nvL(i-1))+(1-beta)*nvT(i-1); //Trend updated with change in level

		}else{
			nvFIL(i) = nvL(n-1)+nvT(n-1)*(i-n+1);
		}
	}
	
	return(wrap(nvFIL));
}

// ---------------------------- Holt-Winters triple exponential smoothing (multiplicative) ----------------------------------
SEXP HW_2(SEXP Y, SEXP S, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) {
	NumericVector nvX(Y); int n = nvX.size(); int f = as<int>(NOUT);
	int s = as<int>(S);
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvL(n); NumericVector nvS(n+s); NumericVector nvFIL(n+f);
	NumericVector nvSTARTVAL(STARTVAL);
	nvL(0) = nvSTARTVAL(0); 
	
	for(int i=0;i<(s+1);i++)nvS(i)=nvSTARTVAL(i+1);
	nvFIL(0) = nvX(0);

	for(int i=1;i<(n+f);i++){
		if(i<n){
			nvFIL(i)=nvL(i-1)*nvS(i); //Predicted/Filtered value for "today"
			nvL(i)=alfa*nvX(i)/nvS(i)+(1-alfa)*nvL(i-1); //Level updated with value of today
			nvS(i+s)=gamma*nvX(i)/nvL(i)+(1-gamma)*nvS(i); //Trend updated with change in level
		}else{
			nvFIL(i) = nvL(n-1)*nvS(n+(i-n)%s);
		}
	}
	
	return(wrap(nvFIL));
}

// ---------------------------- Holt-Winters triple exponential smoothing (multiplicative) ----------------------------------
SEXP HW_TRIPLE_M(SEXP Y, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) {
	NumericVector nvX(Y); int n = nvX.size(); int f = as<int>(NOUT);
	int s = as<int>(S); int o = as<int>(OPTNOUT); 
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1); double beta = nvPARAM(2);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	
	double dL = nvSTARTVAL(0);
	double dL1 = dL;
	double dT = nvSTARTVAL(1);
	for(int i=0;i<s;i++)nvS(i)=nvSTARTVAL(i+2);
	
	int j;
	
	for(int i=1;i<(n+f);i++){
		
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
			nvFIL(i) = (dL+dT)*nvS(i%s);
		}
	}
	
	return(wrap(nvFIL));
}

