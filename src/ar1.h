#ifndef AR1_H
#define AR1_H

#include <RcppArmadillo.h>
using namespace Rcpp;


NumericMatrix ar1_cormatrix_formulated(int param_no, double rho_hat);

#endif
