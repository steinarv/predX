#include "expsmooth.h"
#include "miscf.h"

using namespace Rcpp ;


// --------------------------------- Holt-Winters triple exponential smoothing ----------------------------------------
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


//  ---------------- Holt Winters similar day exponential smoothing with error tracking and outlier detection -----------------
SEXP HW_SIMDAY(SEXP Y, SEXP DAYS, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL, SEXP MULT) {
	
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0; int m = as<int>(MULT);
	int o = as<int>(OPTNOUT);

	
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double beta = nvPARAM(1); double gamma = nvPARAM(2);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); 	//variance
	double dT = nvSTARTVAL(1);	//trend
	double dL = nvSTARTVAL(2); 	//level
	double dL1 = dL;		//holds previous level
	

	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+3);
	

	for(int i=1;i<(n+f);i++){
	std::cout << i << std::endl;

		d = nvDAYS(i);
		
		// If multiplicative........................................................................
		if(m==1){
			
			
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=(dL+dT*(j+1))*nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=(dL+dT)*nvS(d);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvX(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvX(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvX(i) < nvFIL(i, 0) ? 
					xhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					xhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				xhat=nvX(i);
			}
			
			dL1=dL;
			dL=alfa*(xhat/nvS(d))+(1-alfa)*(dL+dT); 		//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(d)=gamma*(xhat/dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			nvFIL(i, 0) = (dL+dT*(i-n))*nvS(d);
		}
		
		
		
		
		// If additive .............................................................................	
		}else{
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=dL+dT*(j+1)+nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=dL+dT+nvS(d);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvX(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvX(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvX(i) < nvFIL(i, 0) ? 
					xhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					xhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				xhat=nvX(i);
			}
			
			dL1=dL;
			dL=alfa*(xhat-nvS(d))+(1-alfa)*(dL+dT); 		//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(d)=gamma*(xhat-dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			nvFIL(i, 0) = dL+dT*(i-n)+nvS(d);
		}
		
	} // En if multiplicative/additive
	
	} // End for
	
	return(wrap(nvFIL));
}
