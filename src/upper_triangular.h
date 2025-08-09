#ifndef EXTRACT_UPPER_TRIANGULAR_H
#define EXTRACT_UPPER_TRIANGULAR_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


arma::vec extract_upper_from_matrix(const arma::mat& kronecker_matrix);
#endif
