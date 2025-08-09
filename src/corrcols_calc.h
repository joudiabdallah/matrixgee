#ifndef CORRCOLS_CALC_H
#define CORRCOLS_CALC_H

#include <RcppArmadillo.h>
using namespace Rcpp;


NumericMatrix correlation_cols_calc_cpp(NumericVector residuals,
                                        int sample_size,
                                        int rows_no,
                                        int cols_no,
                                        int parameters_no,
                                        std::string corstr_cols,
                                        double dispersion_parameter);

#endif
