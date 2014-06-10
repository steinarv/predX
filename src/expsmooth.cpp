#include "expsmooth.h"
#include "miscf.h"

using namespace Rcpp ;


// --------------------------------- Holt-Winters triple exponential smoothing ----------------------------------------
SEXP HW_TRIPLE(SEXP Y, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP STARTVAL, SEXP NOUT, SEXP MULT) {
	NumericVector nvY(Y); int n = nvY.size(); int f = as<int>(NOUT);
	int s = as<int>(S); int o = as<int>(OPTNOUT); int m = as<int>(MULT);
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alpha = nvPARAM(0); double beta = nvPARAM(1); double gamma = nvPARAM(2);
	
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
			dL=alpha*nvY(i)/nvS(i%s)+(1-alpha)*(dL+dT); 	//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(i%s)=gamma*nvY(i)/dL+(1-gamma)*nvS(i%s); 	//Season updated
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
			dL=alpha*(nvY(i)-nvS(i%s))+(1-alpha)*(dL+dT); 	//Level updated with value of today
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			nvS(i%s)=gamma*(nvY(i)-dL)+(1-gamma)*nvS(i%s); 	//Season updated
		}else{
			nvFIL(i, 0) = (dL+dT)+nvS(i%s);
		}
	
	
	}
	}
	
	return(wrap(nvFIL));
}


//  -- Holt Winters similar day exponential smoothing with external level information, error tracking and outlier detection ----------
SEXP HW_SIMDAY(SEXP Y, SEXP DAYS, SEXP L, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL, SEXP MULT) {
	
	NumericVector nvY(Y); NumericVector nvDAYS(DAYS); NumericVector nvL(L);
	
	int n = nvY.size(); int f = nvDAYS.size()-n; int s = as<int>(S); 
	int d = 0; int m = as<int>(MULT); int o = as<int>(OPTNOUT);

	
	double yhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alpha = nvPARAM(0); double beta = nvPARAM(1); double gamma = nvPARAM(2);
	double w1 = nvPARAM(3); double w2 = nvPARAM(4); //w1 = 1 and w2 = 0 when no external level is supplied
	
	/*std::cout << "alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma <<
	", w1: " << w1 << ", w2: " << w2 << std::endl;*/
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	
	double dVAR = nvSTARTVAL(0); 	//variance
	double dL = nvSTARTVAL(1); 	//level
	double dLfil = dL;		//filtered level
	double dL1 = dL;		//holds previous level
	
	double dT = nvSTARTVAL(2);	//trend
	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+3); //Seasonal effects
	

	for(int i=1;i<(n+f);i++){

		d = nvDAYS(i);
		
		// If multiplicative........................................................................
		if(m==1){
			
			
		if(i<n){
			dVAR = 0.06*pow(nvY(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			nvFIL(i, 0)=(w1*dLfil+w2*nvL(i)+dT)*nvS(d);
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=1; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=(w1*dLfil+w2*nvL(i+j)+dT*(j+1))*nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvY(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvY(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvY(i) < nvFIL(i, 0) ? 
					yhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					yhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				yhat=nvY(i);
			}
			
			dL1=dL;
			dLfil=alpha*(yhat/nvS(d))+(1-alpha)*(dLfil+dT);	//Filtered level updated with value of today
			dL=w1*dLfil+w2*nvL(i);
			
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			
			nvS(d)=gamma*(yhat/dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			nvFIL(i, 0) = (w1*dLfil+w2*nvL(i)+dT*(i-n))*nvS(d);
		}
		
		
		
		
		// If additive .............................................................................	
		}else{
		
		if(i<n){
			dVAR = 0.06*pow(nvY(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			nvFIL(i, 0)=w1*dLfil+w2*nvL(i)+dT+nvS(d);
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=1; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=w1*dLfil+w2*nvL(i+j)+dT*(j+1)+nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvY(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvY(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvY(i) < nvFIL(i, 0) ? 
					yhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					yhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				yhat=nvY(i);
			}
			
			dL1=dL;
			dLfil=alpha*(yhat-nvS(d))+(1-alpha)*(dLfil+dT); 	//Filtered level updated with value of today
			dL=w1*dLfil+w2*nvL(i);
			
			dT=beta*(dL-dL1)+(1-beta)*dT;			//Trend updated with change in level
			
			nvS(d)=gamma*(yhat-dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			std::cout << "w1: " << w1 << ", dLfil: " << dLfil << ", w2: " << w2 << ", nvL(i): " << nvL(i) <<
			", dT: " << dT << ", (i-n): " << (i-n) << ", nvS(d): " << nvS(d) <<
			", alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << std::endl;
			nvFIL(i, 0) = w1*dLfil+w2*nvL(i)+dT*(i-n)+nvS(d);
		}
		
	} // En if multiplicative/additive
	
	} // End for
	
	return(wrap(nvFIL));
}



// --- Holt Winters similar day exponential smoothing with external level information, explanatory variable for level, ---
// --- error tracking and outlier detection ------------------------------------------------------------------------------
SEXP HW_SIMDAY_REG(SEXP Y, SEXP DAYS, SEXP L, SEXP S, SEXP X, SEXP OPTNOUT, SEXP PARAM, 
			SEXP THOLD, SEXP STARTVAL, SEXP MULT) {
	
	NumericVector nvY(Y); NumericVector nvDAYS(DAYS); NumericVector nvL(L);
	NumericVector nvX(X);
	
	int n = nvY.size(); int f = nvDAYS.size()-n; int s = as<int>(S); 
	int d = 0; int m = as<int>(MULT); int o = as<int>(OPTNOUT);

	
	double yhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alpha = nvPARAM(0); double beta = nvPARAM(1); double gamma = nvPARAM(2);
	double w1 = nvPARAM(3); double w2 = nvPARAM(4); //w1 = 1 and w2 = 0 when no external level is supplied
	
	/*std::cout << "alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma <<
	", w1: " << w1 << ", w2: " << w2 << std::endl;*/
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	
	double dVAR = nvSTARTVAL(0); 	//variance
	double dL = nvSTARTVAL(1); 	//level
	double dLfil = dL;		//filtered level
	double dLfil_ = dLfil;		//filtered level used in "o" steps prediction for optimization

	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+3); //Seasonal effects
	

	for(int i=1;i<(n+f);i++){

		d = nvDAYS(i);
		
		// If multiplicative........................................................................
		if(m==1){
			
			
		if(i<n){
			dVAR = 0.06*pow(nvY(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			dLfil_ = dLfil+(1-alpha)*(beta*nvX(i)); //filter updated with exvar but not todays obs
			nvFIL(i, 0)=(w1*dLfil_+w2*nvL(i))*nvS(d);
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=1; j<o; j++){
					dLfil_ = dLfil_+(1-alpha)*(beta*nvX(i+j));
					d=nvDAYS(i+j);
					nvFIL(i, j)=(w1*dLfil_+w2*nvL(i+j))*nvS(d); 
					//Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvY(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvY(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvY(i) < nvFIL(i, 0) ? 
					yhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					yhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				yhat=nvY(i);
			}
			
			
			dLfil=alpha*(yhat/nvS(d))+(1-alpha)*(dLfil+beta*nvX(i));	//Filtered level updated with value of today
			dL=w1*dLfil+w2*nvL(i);					//and explanatory variable. The a weighted
										//average between filtered level and
										//level provided by user is computed
										
			nvS(d)=gamma*(yhat/dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			dLfil = dLfil+(1-alpha)*(beta*nvX(i));		//Makes it possible to exclude explanatory var for some i
			nvFIL(i, 0) = (w1*dLfil+w2*nvL(i))*nvS(d);
		}
		
		
		
		
		// If additive .............................................................................	
		}else{
		
		if(i<n){
			dVAR = 0.06*pow(nvY(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			

			dLfil_ = dLfil+(1-alpha)*(beta*nvX(i)); //filter updated with exvar but not todays obs
			nvFIL(i, 0)=(w1*dLfil_+w2*nvL(i))+nvS(d);
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=1; j<o; j++){
					dLfil_ = dLfil_+(1-alpha)*(beta*nvX(i+j));
					d=nvDAYS(i+j);
					nvFIL(i, j)=(w1*dLfil_+w2*nvL(i))+nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}
			
			// If x is more than two standard deviation of one step ahead prediction value we set it 
			// equal to predicted value + thold*sd when updating equations
			if( nvY(i) < (nvFIL(i, 0)-thold*sqrt(dVAR)) || 
							nvY(i) > (nvFIL(i, 0)+thold*sqrt(dVAR)) ){
				
				nvY(i) < nvFIL(i, 0) ? 
					yhat = (nvFIL(i, 0)-thold*sqrt(dVAR)) : 
					yhat = (nvFIL(i, 0)+thold*sqrt(dVAR)) ;
								
			
			}else{
				yhat=nvY(i);
			}
			
			dLfil=alpha*(yhat-nvS(d))+(1-alpha)*(dLfil+beta*nvX(i));	//Filtered level updated with value of today
			dL=w1*dLfil+w2*nvL(i);					//and explanatory variable. The a weighted
										//average between filtered level and
										//level provided by user is computed
			
			nvS(d)=gamma*(yhat-dL)+(1-gamma)*nvS(d); 	//Seasonal component updated 
			
		}else{
			std::cout << "w1: " << w1 << ", dLfil: " << dLfil << ", w2: " << w2 << ", nvL(i): " << nvL(i) <<
			", nvX(i): " << nvX(i) << ", (i-n): " << (i-n) << ", nvS(d): " << nvS(d) <<
			", alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << std::endl;
			dLfil = dLfil+(1-alpha)*(beta*nvX(i));		//Makes it possible to exclude explanatory var for some i
			nvFIL(i, 0) = w1*dLfil+w2*nvL(i)+nvS(d);
		}
		
	} // En if multiplicative/additive
	
	} // End for
	
	return(wrap(nvFIL));
}
