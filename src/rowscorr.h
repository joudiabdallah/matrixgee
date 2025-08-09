#ifndef ROWSCORR_H
#define ROWSCORR_H

#include <RcppArmadillo.h>
using namespace Rcpp;


NumericMatrix correlation_rows_calc_cpp(NumericVector residuals,
                                        int sample_size,
                                        int rows_no,
                                        int cols_no,
                                        int parameters_no,
                                        std::string corstr_rows,
                                        double dispersion_parameter);

#endif

