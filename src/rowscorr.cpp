// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "ar1.h"

// [[Rcpp::export]]
NumericMatrix correlation_rows_calc_cpp(NumericVector residuals,
                                        int sample_size,
                                        int rows_no,
                                        int cols_no,
                                        int parameters_no,
                                        std::string corstr_rows,
                                        double dispersion_parameter) {
        if (corstr_rows == "independence") {
                NumericMatrix corr(rows_no, rows_no);
                for (int i = 0; i < rows_no; ++i) {
                        corr(i, i) = 1.0;
                }
                return corr;
        }

        int total_elements = sample_size * rows_no * cols_no;
        NumericMatrix residuals_matrix_all(rows_no, cols_no * sample_size);

        for (int i = 0; i < total_elements; ++i) {
                int sample_idx = i / (rows_no * cols_no);
                int within_sample_idx = i % (rows_no * cols_no);
                int row_idx = within_sample_idx % rows_no;
                int col_idx = (within_sample_idx / rows_no) + sample_idx * cols_no;
                residuals_matrix_all(row_idx, col_idx) = residuals[i];
        }

        double total_numerator = 0.0;
        double denominator = 1.0;

        if (corstr_rows == "Exchangeable") {
                denominator = ((0.5 * sample_size * rows_no * cols_no * (rows_no - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        for (int j = 0; j < cols_no; ++j) {
                                for (int t = 0; t < rows_no - 1; ++t) {
                                        for (int s = t + 1; s < rows_no; ++s) {
                                                total_numerator += residuals_matrix_all(t, start_col + j) * residuals_matrix_all(s, start_col + j);
                                        }
                                }
                        }
                }
        } else if (corstr_rows == "ar1") {
                denominator = ((sample_size * cols_no * (rows_no - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        for (int t = 0; t < rows_no - 1; ++t) {
                                for (int j = 0; j < cols_no; ++j) {
                                        total_numerator += residuals_matrix_all(t, start_col + j) * residuals_matrix_all(t + 1, start_col + j);
                                }
                        }
                }
        }

        double rho_hat_trial = total_numerator / denominator;
        double lower_bound = -1.0 / (rows_no - 1);
        double rho_hat = (rho_hat_trial <= 1.0 && rho_hat_trial >= lower_bound) ? rho_hat_trial : 0.7;

        if (corstr_rows == "Exchangeable") {
                NumericMatrix corr(rows_no, rows_no);
                for (int i = 0; i < rows_no; ++i) {
                        for (int j = 0; j < rows_no; ++j) {
                                corr(i, j) = (i == j) ? 1.0 : rho_hat;
                        }
                }
                return corr;
        } else if (corstr_rows == "ar1") {
                return ar1_cormatrix_formulated(rows_no, rho_hat);
        }

        return NumericMatrix(rows_no, rows_no);
}
