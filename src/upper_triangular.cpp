#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec extract_upper_from_matrix(const arma::mat& kronecker_matrix) {
        int n_rows = kronecker_matrix.n_rows;
        int n_cols = kronecker_matrix.n_cols;

        std::vector<double> upper_values;

        for (int i = 0; i < n_rows; ++i) {
                for (int j = i + 1; j < n_cols; ++j) {
                        upper_values.push_back(kronecker_matrix(i, j));
                }
        }

        return arma::vec(upper_values);
}

