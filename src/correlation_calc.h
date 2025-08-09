#ifndef CORRELATION_GENERAL_H
#define CORRELATION_GENERAL_H

#include <Rcpp.h>
using namespace Rcpp;
NumericMatrix correlation_general(NumericVector residuals,
                                  int sample_size,
                                  int rows_no,
                                  int cols_no,
                                  int parameters_no,
                                  std::string corstr,
                                  double dispersion_parameter);

#endif
