#ifndef _predX_var_H
#define _predX_var_H

#include <RcppArmadillo.h>

RcppExport SEXP EXPSMOOTH1(SEXP X, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) ;
RcppExport SEXP EXPSMOOTH2(SEXP X, SEXP PARAM, SEXP THOLD, SEXP STARTVAL, SEXP NOUT) ;
RcppExport SEXP EXPSMOOTH3(SEXP X, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) ;

RcppExport SEXP SEASEXPSMOOTH(SEXP X, SEXP S, SEXP PARAM, SEXP STARTVAL, SEXP NOUT) ;
RcppExport SEXP SIMDAYEXPSMOOTH(SEXP X, SEXP DAYS, SEXP S, SEXP PARAM, SEXP STARTVAL) ;
RcppExport SEXP RMSIMDAYEXPSMOOTH(SEXP X, SEXP DAYS, SEXP S, SEXP PARAM, SEXP THOLD, SEXP STARTVAL) ;




#endif
