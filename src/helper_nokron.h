#ifndef VECTOR_TO_LONG_DATAFRAME_H
#define VECTOR_TO_LONG_DATAFRAME_H

#include <RcppArmadillo.h>


Rcpp::DataFrame vector_to_long_dataframe_cpp(Rcpp::NumericMatrix data, int sample_size, int rows_no, int cols_no);

#endif
