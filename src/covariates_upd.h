#ifndef COVARIATES_UPD_H
#define COVARIATES_UPD_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List covariates_upd(SEXP covariates, int sample_size, int rows_no, int cols_no);
#endif
