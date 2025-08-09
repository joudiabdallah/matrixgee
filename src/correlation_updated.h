#ifndef CORRELATION_GENERAL_UPDATED_H
#define CORRELATION_GENERAL_UPDATED_H

#include <Rcpp.h>
using namespace Rcpp;
NumericMatrix correlation_general_updated(NumericVector residuals,
                                          int sample_size,
                                          int rows_no,
                                          int cols_no,
                                          int parameters_no,
                                          std::string corstr,
                                          double dispersion_parameter);

#endif
