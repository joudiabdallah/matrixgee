// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "ar1.h"

// [[Rcpp::export]]
NumericMatrix correlation_cols_updated_cpp(NumericVector residuals,
                                           int sample_size,
                                           int rows_no,
                                           int cols_no,
                                           int parameters_no,
                                           std::string corstr_cols,
                                           double dispersion_parameter) {

        if (corstr_cols == "independence") {
                NumericMatrix corr(cols_no, cols_no);
                for (int i = 0; i < cols_no; ++i) {
                        corr(i, i) = 1.0;
                }
                return corr;
        }


        int total_elements = sample_size * rows_no * cols_no;
        NumericMatrix residuals_matrix_all(rows_no, cols_no * sample_size);

        for (int i = 0; i < total_elements; ++i) {
                int sample_idx = i / (rows_no * cols_no);
                int within_sample = i % (rows_no * cols_no);
                int row_idx = within_sample % rows_no;
                int col_idx = (within_sample / rows_no) + sample_idx * cols_no;
                residuals_matrix_all(row_idx, col_idx) = residuals[i];
        }

        // Calculate dispersion
        //double sum_sq = sum(residuals * residuals);
        //double a = sample_size * rows_no * cols_no - parameters_no;
        //double dispersion_parameter = sum_sq / a;

        double denominator_cols = 1.0;
        double total_numerator_cols = 0.0;

        // Exchangeable structure
        if (corstr_cols == "Exchangeable") {
                denominator_cols = ((0.5 * sample_size * rows_no * cols_no * (cols_no - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        for (int t = 0; t < rows_no; ++t) {
                                for (int j = 0; j < cols_no - 1; ++j) {
                                        for (int k = j + 1; k < cols_no; ++k) {
                                                total_numerator_cols += residuals_matrix_all(t, start_col + j) *
                                                        residuals_matrix_all(t, start_col + k);
                                        }
                                }
                        }
                }

        } else if (corstr_cols == "ar1") {
                denominator_cols = ((sample_size * rows_no * (cols_no - 1)) - parameters_no) * dispersion_parameter;

                for (int i = 0; i < sample_size; ++i) {
                        int start_col = i * cols_no;
                        for (int t = 0; t < rows_no; ++t) {
                                for (int c = 0; c < cols_no - 1; ++c) {
                                        total_numerator_cols += residuals_matrix_all(t, start_col + c) *
                                                residuals_matrix_all(t, start_col + c + 1);
                                }
                        }
                }
        }

        // Estimate rho
        double rho_hat_trial = total_numerator_cols / denominator_cols;
        double lower_bound = -1.0 / (cols_no - 1);
        double rho_hat = (rho_hat_trial >= lower_bound && rho_hat_trial <= 1.0) ? rho_hat_trial : 0.87;


        NumericMatrix corr(cols_no, cols_no);
        if (corstr_cols == "Exchangeable") {
                for (int i = 0; i < cols_no; ++i) {
                        for (int j = 0; j < cols_no; ++j) {
                                corr(i, j) = (i == j) ? 1.0 : rho_hat;
                        }
                }
        } else if (corstr_cols == "ar1") {
                corr = ar1_cormatrix_formulated(cols_no, rho_hat);
        }

        return corr;
}
