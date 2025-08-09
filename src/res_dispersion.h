#ifndef RES_DISPERSION_H
#define RES_DISPERSION_H

#include <Rcpp.h>
using namespace Rcpp;


List residuals_dispersion_helper(NumericVector model_residuals,
                                 int sample_size,
                                 int rows_no,
                                 int cols_no,
                                 int parameters_no);

#endif
