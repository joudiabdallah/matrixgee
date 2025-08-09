// helper.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//#include "helper_header.h"



// [[Rcpp::export]]
NumericMatrix ar1_cormatrix_formulated(int param_no, double rho_hat) {

        NumericMatrix correlation_matrix(param_no, param_no);

        for (int i = 0; i < param_no; i++) {
                for (int j = 0; j < param_no; j++) {
                        correlation_matrix(i, j) = (i == j) ? 1.0 : pow(rho_hat, std::abs(i - j));
                }
        }

        return correlation_matrix;
}
