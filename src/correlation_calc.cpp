// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "ar1.h"


// [[Rcpp::export]]

NumericMatrix correlation_general(NumericVector residuals,
                                  int sample_size,
                                  int rows_no,
                                  int cols_no,
                                  int parameters_no,
                                  std::string corstr,
                                  double dispersion_parameter) {
        //return identity matrix if the corstr is independence
        if (corstr == "independence") {
                int dim = cols_no * rows_no;
                NumericMatrix corr(dim, dim);
                for (int i = 0; i < dim; ++i) {
                        corr(i, i) = 1.0;
                }
                return corr;
        }
        // Step 1: Reshape residuals into r x (c * N) matrix
        int total_elements = sample_size * rows_no * cols_no;
        NumericMatrix residuals_matrix_all(rows_no, cols_no * sample_size);

        for (int i = 0; i < total_elements; ++i) {
                int sample_idx = i / (rows_no * cols_no);
                int within_sample = i % (rows_no * cols_no);
                int row_idx = within_sample % rows_no;
                int col_idx = (within_sample / rows_no) + sample_idx * cols_no;
                residuals_matrix_all(row_idx, col_idx) = residuals[i];
        }

        // Initialize shared variables
        double total_numerator = 0.0;
        double denominator = 1.0;
        int dim = rows_no * cols_no;

        // Step 2: Compute correlation
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

                        arma::vec resid_vec = arma::vectorise(resid_matrix_i); // rc x 1

                        for (int j = 0; j < dim - 1; ++j) {
                                for (int k = j + 1; k < dim; ++k) {
                                        total_numerator += resid_vec(j) * resid_vec(k);
                                }
                        }
                }

                double rho_hat_trial = total_numerator / denominator;
                double lower_bound = -1.0 / (dim - 1);
                double rho_hat = (rho_hat_trial <= 1.0 && rho_hat_trial >= lower_bound) ? rho_hat_trial : 0.8;

                // Construct exchangeable matrix
                NumericMatrix corr(dim, dim);
                for (int i = 0; i < dim; ++i) {
                        for (int j = 0; j < dim; ++j) {
                                corr(i, j) = (i == j) ? 1.0 : rho_hat;
                        }
                }
                return corr;

        }
        else if (corstr == "ar1") {
                denominator = ((sample_size * (dim - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        arma::mat resid_matrix_i(rows_no, cols_no);
                        for (int r = 0; r < rows_no; ++r) {
                                for (int c = 0; c < cols_no; ++c) {
                                        resid_matrix_i(r, c) = residuals_matrix_all(r, start_col + c);
                                }
                        }

                        arma::vec resid_vec = arma::vectorise(resid_matrix_i); // rc x 1

                        for (int j = 0; j < dim - 1; ++j) {
                                total_numerator += resid_vec(j) * resid_vec(j + 1); // lag-1 product
                        }
                }

                double rho_hat_trial = total_numerator / denominator;
                double lower_bound = -1.0;
                double rho_hat = (rho_hat_trial <= 1.0 && rho_hat_trial >= lower_bound) ? rho_hat_trial : 0.5;

                NumericMatrix corr = ar1_cormatrix_formulated(dim, rho_hat);
                return corr;
        }


        else {
                Rcpp::stop("Unknown correlation structure: " + corstr);
        }

}
