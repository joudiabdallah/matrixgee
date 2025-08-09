#ifndef DESIGN_MATRIX_CPP_H
#define DESIGN_MATRIX_CPP_H

#include <RcppArmadillo.h>

// Declare the main interface function
arma::mat create_design_matrix_cpp(int N, int r, int c,
                                   Rcpp::CharacterVector intercepts,
                                   arma::mat design_covblock);

// Declare internal helpers (inline or available for separate compilation)
bool contains_intercept(Rcpp::CharacterVector intercepts, const std::string& target);

arma::mat scalar_intercept_block(int N, int r, int c);
arma::mat row_intercept_block(int N, int r, int c);
arma::mat col_intercept_block(int N, int r, int c);
arma::mat matrix_intercept_block(int N, int r, int c);

#endif // DESIGN_MATRIX_CPP_H
