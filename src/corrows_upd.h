#ifndef CORROWS_UPD_H
#define CORROWS_UPD_H

#include <RcppArmadillo.h>
using namespace Rcpp;


NumericMatrix correlation_rows_updated_cpp(NumericVector residuals,
                                           int sample_size,
                                           int rows_no,
                                           int cols_no,
                                           int parameters_no,
                                           std::string corstr_rows,
                                           double dispersion);

#endif
