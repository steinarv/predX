#ifndef _predX_expmod_H
#define _predX_expmod_H

#include <RcppArmadillo.h>


RcppExport SEXP HW_TRIPLE(SEXP Y, SEXP S, SEXP PARAM, SEXP OPTNOUT, SEXP STARTVAL, SEXP NOUT, SEXP MULT);
RcppExport SEXP HW_SIMDAY(SEXP Y, SEXP DAYS, SEXP L, SEXP S, SEXP OPTNOUT, SEXP PARAM, SEXP THOLD, 
                            SEXP STARTVAL, SEXP MULT);
RcppExport SEXP HW_SIMDAY_REG(SEXP Y, SEXP DAYS, SEXP L, SEXP S, SEXP X, SEXP OPTNOUT, SEXP PARAM,
                    SEXP THOLD, SEXP STARTVAL, SEXP MULT);








#endif
