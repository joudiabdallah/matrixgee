#ifndef HELPER_DATA_H
#define HELPER_DATA_H

#include <RcppArmadillo.h>
using namespace Rcpp;


DataFrame matrix_to_long_dataframe_cpp(NumericMatrix data, int sample_size, int rows_no, int cols_no);

#endif
