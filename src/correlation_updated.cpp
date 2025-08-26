// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "ar1.h"

// [[Rcpp::export]]
NumericMatrix correlation_general_updated(NumericVector residuals,
                                          int sample_size,
                                          int rows_no,
                                          int cols_no,
                                          int parameters_no,
                                          std::string corstr,
                                          double dispersion_parameter) {
        int dim = rows_no * cols_no;

        if (corstr == "independence") {
                NumericMatrix corr(dim, dim);
                for (int i = 0; i < dim; ++i) {
                        corr(i, i) = 1.0;
                }
                return corr;
        }

        //  reshape residual vector into residuals_matrix_all (r x c*N)
        int total_elements = sample_size * rows_no * cols_no;
        NumericMatrix residuals_matrix_all(rows_no, cols_no * sample_size);

        for (int i = 0; i < total_elements; ++i) {
                int sample_idx = i / (rows_no * cols_no);
                int within_sample = i % (rows_no * cols_no);
                int row_idx = within_sample % rows_no;
                int col_idx = (within_sample / rows_no) + sample_idx * cols_no;
                residuals_matrix_all(row_idx, col_idx) = residuals[i];
        }

        double denominator = 1.0;
        double total_numerator = 0.0;

        if (corstr == "Exchangeable") {
                denominator = ((0.5 * sample_size * dim * (dim - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        arma::mat resid_matrix_i(rows_no, cols_no);
                        for (int r = 0; r < rows_no; ++r) {
                                for (int c = 0; c < cols_no; ++c) {
                                        resid_matrix_i(r, c) = residuals_matrix_all(r, start_col + c);
                                }
                        }

                        arma::vec resid_vec = vectorise(resid_matrix_i);

                        for (int j = 0; j < dim - 1; ++j) {
                                for (int k = j + 1; k < dim; ++k) {
                                        total_numerator += resid_vec(j) * resid_vec(k);
                                }
                        }
                }

        } else if (corstr == "ar1") {
                denominator = ((sample_size * (dim - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        arma::mat resid_matrix_i(rows_no, cols_no);
                        for (int r = 0; r < rows_no; ++r) {
                                for (int c = 0; c < cols_no; ++c) {
                                        resid_matrix_i(r, c) = residuals_matrix_all(r, start_col + c);
                                }
                        }

                        arma::vec resid_vec = vectorise(resid_matrix_i);
                        for (int j = 0; j < dim - 1; ++j) {
                                total_numerator += resid_vec(j) * resid_vec(j + 1);
                        }
                }
        }

        //  Estimate rho
        double rho_hat_trial = total_numerator / denominator;
        double lower_bound = -1.0 / (dim - 1);
        double rho_hat = (rho_hat_trial >= lower_bound && rho_hat_trial <= 1.0) ? rho_hat_trial : 0.87;

        //Return correlation matrix
        NumericMatrix corr(dim, dim);
        if (corstr == "Exchangeable") {
                for (int i = 0; i < dim; ++i) {
                        for (int j = 0; j < dim; ++j) {
                                corr(i, j) = (i == j) ? 1.0 : rho_hat;
                        }
                }
        } else if (corstr == "ar1") {
                corr = ar1_cormatrix_formulated(dim, rho_hat);
        }

        return corr;

}
