#include "forecastingModels.h"
#include "miscf.h"

using namespace Rcpp ;

// ---------------------------- Exponential smoothing 1 --------------------------------------------------------------
SEXP EXPSMOOTH1(SEXP Y, SEXP PARAM, SEXP STARTVAL, SEXP NOUT){
	NumericVector nvX(Y); int n=nvX.size();
	int f=as<int>(NOUT); NumericVector nvFIL(n+f); //f = number of forecasts to be returned
	
	double dLAMBDA=as<double>(PARAM);
	
	nvFIL(0)=as<double>(STARTVAL);
	
	for(int i=1;i<(n+f);i++){
		if(i<=n){
		nvFIL(i)=(1-dLAMBDA)*nvX(i-1)+dLAMBDA*nvFIL(i-1);
		}else{
		nvFIL(i)=nvFIL(n-1);
		}
	
	}

	return(wrap(nvFIL));
}

// ---------------------------- Exponential smoothing 2 --------------------------------------------------------------
SEXP EXPSMOOTH2(SEXP Y, SEXP PARAM, SEXP THOLD, SEXP STARTVAL, SEXP NOUT){
	NumericVector nvX(Y); int n=nvX.size();
	int f=as<int>(NOUT); NumericVector nvFIL(n+f); //f = number of forecasts to be returned
	
	NumericVector nvTHOLD(THOLD);
			
	double dLAMBDA=as<double>(PARAM);
	
	nvFIL(0)=as<double>(STARTVAL);
	
	for(int i=1;i<(n+f);i++){
		if(i<n){

			if(nvX(i)>0){
				if(nvX(i)>(nvFIL(i-1)+nvTHOLD(i))){
					//std::cout << "Value before: " << nvX(i) << " E: " << nvFIL(i-1) << "  THOLD: " << nvTHOLD(i) << std::endl;
					nvFIL(i)=(1-dLAMBDA)*(nvFIL(i-1)+nvTHOLD(i))+dLAMBDA*nvFIL(i-1);
					//nvX(i)=nvFIL(i-1)+nvTHOLD(i);
					//std::cout << "Value before: " << nvX(i) << " E: " << nvFIL(i-1) << "  THOLD: " << nvTHOLD(i) << std::endl;
				}else{
					nvFIL(i)=(1-dLAMBDA)*nvX(i-1)+dLAMBDA*nvFIL(i-1);
				}
			}else{
				if(nvX(i)<(nvFIL(i-1)-nvTHOLD(i))){
					//std::cout << "Value before: " << nvX(i) << " E: " << nvFIL(i-1) << "  THOLD: " << nvTHOLD(i) << std::endl;
					nvFIL(i)=(1-dLAMBDA)*(nvFIL(i-1)-nvTHOLD(i))+dLAMBDA*nvFIL(i-1);
					//nvX(i)=nvFIL(i-1)-nvTHOLD(i);
					//std::cout << "Value before: " << nvX(i) << " E: " << nvFIL(i-1) << "  THOLD: " << nvTHOLD(i) << std::endl;
				}else{
					nvFIL(i)=(1-dLAMBDA)*nvX(i-1)+dLAMBDA*nvFIL(i-1);
				}
			}
				
		//nvFIL(i)=(1-dLAMBDA)*nvX(i-1)+dLAMBDA*nvFIL(i-1);
		}else{
		nvFIL(i)=nvFIL(n-1);
		}
	
	}

	return(wrap(nvFIL));
}

// ---------------------------- Holt Winters model 1 -----------------------------------------------------------------
SEXP EXPSMOOTH3(SEXP Y, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) {
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
// -------------------------------------------------------------------------------------------------------------------


// ---------------------------- simple seasonal Holt Winters model ---------------------------------------------------
SEXP SEASEXPSMOOTH(SEXP Y, SEXP S, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) {
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
// -------------------------------------------------------------------------------------------------------------------

// ------------------------------ simelar day Holt Winters model -----------------------------------------------------
SEXP SIMDAYEXPSMOOTH(SEXP Y, SEXP DAYS, SEXP S, SEXP PARAM, SEXP STARTVAL) {
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvL(n); NumericVector nvS(s); NumericVector nvFIL(n+f);
	NumericVector nvSTARTVAL(STARTVAL);
	nvL(0) = nvSTARTVAL(0); 
	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+1);
	nvFIL(0) = nvX(0);

	for(int i=1;i<(n+f);i++){
		d = nvDAYS(i);
		
		if(i<n){
			nvFIL(i)=nvL(i-1)*nvS(d); //Predicted/Filtered value for "today"
			nvL(i)=alfa*nvX(i)/nvS(d)+(1-alfa)*nvL(i-1); //Level updated with value of today
			nvS(d)=gamma*nvX(i)/nvL(i)+(1-gamma)*nvS(d); //Trend updated with change in level
		}else{
			nvFIL(i) = nvL(n-1)*nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------


// ------------------------------ robust simelar day Holt Winters model -----------------------------------------------
SEXP RMSIMDAYEXPSMOOTH(SEXP Y, SEXP DAYS, SEXP S, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) {
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvS(s); NumericVector nvFIL(n+f);
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); //variance
	double dL = nvSTARTVAL(1); //level
	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+2);
	nvFIL(0) = nvX(0);

	for(int i=1;i<(n+f);i++){
	
		d = nvDAYS(i);
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1), 2)+0.94*dVAR;
			nvFIL(i)=dL*nvS(d); //Predicted/Filtered value for "today"
			
			// If x is more than two standard deviation of filtered value we set it 
			// equal to filtered value when updating equations
			if( nvX(i) < (nvFIL(i)-thold*sqrt(dVAR)) || 
							nvX(i) > (nvFIL(i)+thold*sqrt(dVAR)) ){
								
				nvX(i) < nvFIL(i) ? 
					xhat = (nvFIL(i)-thold*sqrt(dVAR)) : 
					xhat = (nvFIL(i)+thold*sqrt(dVAR));
								
				dL=alfa*xhat/nvS(d)+(1-alfa)*dL; 
				nvS(d)=gamma*xhat/dL+(1-gamma)*nvS(d); 
			
			}else{
				dL=alfa*nvX(i)/nvS(d)+(1-alfa)*dL; //Level updated with value of today
				nvS(d)=gamma*nvX(i)/dL+(1-gamma)*nvS(d); //Trend updated with change in level
			
			}
		}else{
			nvFIL(i) = dL*nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------


// ------------------------------ robust simelar day Holt Winters model -----------------------------------------------
// ----------------- This model calculates OPTNOUT predictions at each step in filtration -----------------------------
// ------------------------ allowing for more flexible optimization routines ------------------------------------------
SEXP SIMDAYSMOOTH(SEXP Y, SEXP DAYS, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) {
	
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	int o = as<int>(OPTNOUT);
	
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); //variance
	double dL = nvSTARTVAL(1); //level
	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+2);
	

	for(int i=1;i<(n+f);i++){
	
		d = nvDAYS(i);
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			//nvFIL(i, 0)=dL*nvS(d); //Predicted/Filtered value for "today"
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=dL*nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=dL*nvS(d);
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
			
			
			dL=alfa*xhat/nvS(d)+(1-alfa)*dL; 	//Level updated with value of today
			nvS(d)=gamma*xhat/dL+(1-gamma)*nvS(d); 	//Trend updated with change in level
			
			
		}else{
			nvFIL(i, 0) = dL*nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------



// ------------------------------ robust simelar day Holt Winters model -----------------------------------------------
// ----------------- This model calculates OPTNOUT predictions at each step in filtration -----------------------------
// ------------------------ allowing for more flexible optimization routines ------------------------------------------
SEXP SIMDAYSMOOTHLEVELINPUT(SEXP Y, SEXP L, SEXP DAYS, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) {
	
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	int o = as<int>(OPTNOUT);
	
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); //variance
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+2);
	
	
	// One smoothed level and one external level input is combined based on historical accuracy
	double dL = nvSTARTVAL(1); //level
	double dL1 = dL; // smoothed level
	NumericVector nvL(L);
	double dL2 = 0; // input level
	double r1=0.5; // accuracy of smoothed level
	double r2=0.5; // accuracy of input level
	double w1=0.5; // Weight assigned to smoothed level
	double w2=0.5; // Weight assigned to input level
	
	
	for(int i=1;i<(n+f);i++){
	
		d = nvDAYS(i);
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			//nvFIL(i, 0)=dL*nvS(d); //Predicted/Filtered value for "today"
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=dL*nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=dL*nvS(d);
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
			
			//Calculate the optimal combination of the two estimates of the level
			r1=0.1*pow(xhat/nvS(d)-dL1, 2)+0.9*r1; 	//Accuracy of smoothed level
			r2=0.1*pow(xhat/nvS(d)-dL2, 2)+0.9*r2;	//Accuracy of input level
			
			w1=(1-r1/(r1+r2)); w2=1-w1; //Weights
			
			dL1=alfa*xhat/nvS(d)+(1-alfa)*dL1; 	//Smoothed level updated with value of today
			dL2=nvL(i);
			dL=w1*dL1+w2*dL2;
			
			nvS(d)=gamma*xhat/dL+(1-gamma)*nvS(d); 	//Trend updated with change in level
			
		}else{
			dL=w1*dL1+w2*nvL(i);
			nvFIL(i, 0) = dL*nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------


// --------------------------------------------------------------------------------------------------------------------
// ------------------------------ slight modification on filtering technique ------------------------------------------
// --------------------------------------------------------------------------------------------------------------------

SEXP SIMDAYSMOOTH2(SEXP Y, SEXP DAYS, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) {
	
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	int o = as<int>(OPTNOUT);
	
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); //variance
	double dL = nvSTARTVAL(1); //level
	
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+2);
	

	for(int i=1;i<(n+f);i++){
	
		d = nvDAYS(i);
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			//nvFIL(i, 0)=dL*nvS(d); //Predicted/Filtered value for "today"
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=dL+nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=dL+nvS(d);
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
			
			
			dL=alfa*(xhat-nvS(d))+(1-alfa)*dL; 	//Level updated with value of today
			nvS(d)=gamma*(xhat-dL)+(1-gamma)*nvS(d); 	//Trend updated with change in level
			
			
		}else{
			nvFIL(i, 0) = dL+nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------



// ------------------------------ robust simelar day Holt Winters model -----------------------------------------------
// ----------------- This model calculates OPTNOUT predictions at each step in filtration -----------------------------
// ------------------------ allowing for more flexible optimization routines ------------------------------------------
SEXP SIMDAYSMOOTHLEVELINPUT2(SEXP Y, SEXP L, SEXP DAYS, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) {
	
	NumericVector nvX(Y); NumericVector nvDAYS(DAYS); int n = nvX.size(); 
	int f = nvDAYS.size()-n; int s = as<int>(S); int d = 0;
	int o = as<int>(OPTNOUT);
	
	double xhat = 0; // Normalized x when outliers detected
	double thold = as<double>(THOLD); //Number of standard deviations for treshold (0 < > 4)
	
	NumericVector nvPARAM(PARAM); unityFunc(nvPARAM);
	double alfa = nvPARAM(0); double gamma = nvPARAM(1);
	
	//std::cout << "alfa: " << alfa << ", gamma: " << gamma << std::endl;
	
	NumericVector nvS(s); NumericMatrix nvFIL(n+f, o);
	
	NumericVector nvSTARTVAL(STARTVAL);
	double dVAR = nvSTARTVAL(0); //variance
	for(int i=0;i<(s);i++)nvS(i)=nvSTARTVAL(i+2);
	
	
	// One smoothed level and one external level input is combined based on historical accuracy
	double dL = nvSTARTVAL(1); //level
	double dL1 = dL; // smoothed level
	NumericVector nvL(L);
	double dL2 = 0; // input level
	double r1=0.5; // accuracy of smoothed level
	double r2=0.5; // accuracy of input level
	double w1=0.5; // Weight assigned to smoothed level
	double w2=0.5; // Weight assigned to input level
	
	
	for(int i=1;i<(n+f);i++){
	
		d = nvDAYS(i);
		
		if(i<n){
			dVAR = 0.06*pow(nvX(i-1)-nvFIL(i-1, 0), 2)+0.94*dVAR;
			
			//nvFIL(i, 0)=dL*nvS(d); //Predicted/Filtered value for "today"
			
			if(i<=(n-o)){ //Make predictions "o" steps ahead
				for(int j=0; j<o; j++){
					d=nvDAYS(i+j);
					nvFIL(i, j)=dL+nvS(d); //Predicted/Filtered value for "today"
				}
				d = nvDAYS(i);
			}else{
				nvFIL(i, 0)=dL+nvS(d);
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
			
			//Calculate the optimal combination of the two estimates of the level
			r1=0.1*pow(xhat-nvS(d)-dL1, 2)+0.9*r1; 	//Accuracy of smoothed level
			r2=0.1*pow(xhat-nvS(d)-dL2, 2)+0.9*r2;	//Accuracy of input level
			
		//	w1=(1-r1/(r1+r2)); w2=1-w1; //Weights
			
			// ------------ delete ------------
			//std::cout << "w1 and w2: " << w1 << ", " << w2 << std::endl;
			//Probably best to predefine at average or optimize w1 and w2
			
			dL1=alfa*(xhat-nvS(d))+(1-alfa)*dL1; 	//Smoothed level updated with value of today
			dL2=nvL(i);
			dL=w1*dL1+w2*dL2;
			
			nvS(d)=gamma*(xhat-dL)+(1-gamma)*nvS(d); 	//Trend updated with change in level
			
		}else{
			dL=w1*dL1+w2*nvL(i);
			nvFIL(i, 0) = dL+nvS(d);
		}
	}
	
	return(wrap(nvFIL));
}
// -------------------------------------------------------------------------------------------------------------------































/*
// ---------------------------- Returns conditional variance for standard GARCH model --------------------------------------
SEXP HTstdGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) {
	double dH0 = as<double>(H0); double dNout = as<double>(Nout);
	NumericVector nvX(X); NumericVector nvParam(Param);
	
	int n = nvX.size();
	double a0 = ::exp(nvParam(0)); double a1 = ::exp(nvParam(1));
	double b0 = ::exp(nvParam(2));
  
	NumericVector H(n+dNout);
	H(0) = dH0;
  
	for(int i=1;i<(n+dNout);i++){
		//std::cout << "a0: " << a0 << " a1: " << a1 << " b0: " << b0 << std::endl;
		//std::cout << "pow(nvX(i-1),2)   " << pow(nvX(i-1),2) << std::endl;
		//std::cout << "a0+a1*pow(nvX(i-1),2)+b0*H(i-1)   " << a0+a1*pow(nvX(i-1),2)+b0*H(i-1) << std::endl;
		if(i<=n)H(i)=a0+a1*pow(nvX(i-1),2)+b0*H(i-1);
		else H(i)=a0+(a1+b0)*H(i-1); //Out of sample predictions
	}
  
	return wrap(H);
}
// -------------------------------------------------------------------------------------------------------------------

// --------------------------------- Returns log likelihood for EGARCH model -----------------------------------------
SEXP LLeGARCH(SEXP X, SEXP H0, SEXP Param) {
	double dH0 = as<double>(H0);
	NumericVector nvX(X); NumericVector nvParam(Param);
	
	int n = nvX.size();
	double a0 = nvParam(0); double a1 = nvParam(1);
	double gamma = nvParam(2); double b0 = nvParam(3);
  
	NumericVector logH(n);
	logH(0) = ::log(dH0);
	
	double logL=0;
  
	for(int i=1;i<n;i++){
		logH(i)=a0+a1*(std::abs(nvX(i-1))+gamma*nvX(i-1))/sqrt(::exp(logH(i-1)))+b0*logH(i-1);
		logL+=logH(i)+pow(nvX(i),2)/::exp(logH(i));		//The correct way is to subtract rather than add 
	}
  
	return wrap(logL);
}

// -------------------------------- Returns times series of conditional EGARCH variance ---------------------------------
SEXP HTeGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) {
	double dH0 = as<double>(H0); double dNout = as<double>(Nout);
	NumericVector nvX(X); NumericVector nvParam(Param);	
	
	int n = nvX.size();
	double a0 = nvParam(0); double a1 = nvParam(1);
	double gamma = nvParam(2); double b0 = nvParam(3);
  
	NumericVector logH(n+dNout);
	logH(0) = ::log(dH0);
	
	//std::cout << "dH0    " << dH0 << std::endl;
	//std::cout << "::log(dH0)   " << ::log(dH0)  << std::endl;

	for(int i=1;i<(n+dNout);i++){
		if(i<=n) logH(i)=a0+a1*(std::abs(nvX(i-1))+gamma*nvX(i-1))/sqrt(::exp(logH(i-1)))+b0*logH(i-1);
		else logH(i)=a0+a1*sqrt(2/M_PI)+b0*logH(i-1);
	}
  
	return wrap(logH);
}*/
