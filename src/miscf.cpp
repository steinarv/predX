#include "miscf.h"

using namespace Rcpp ;

SEXP rollavg(SEXP x, SEXP dim){
  arma::vec orig = as<arma::vec>(x);
  int n = as<int>(dim);  
  int obs = orig.n_elem;
  NumericVector out(obs-n+1);
                  
  for(int i=n-1; i<obs; i++){
    out(i-n+1) = arma::sum(orig.rows(i-n+1,i))/n;
  }
  
  return(wrap(out));
}