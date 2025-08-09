// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "helper_data.h"
#include "res_dispersion.h"
#include "corrcols_calc.h"
#include "rowscorr.h"
#include "corrcols_upd.h"
#include "fitgee.h"
#include "ar1.h"
#include "corrows_upd.h"
#include "cov.h"
#include "link.h"
#include "nuissance.h"
#include "utils.h"
#include "variance.h"
#include "upper_triangular.h"
#include "covariates_upd.h"
#include "design_matrix_cpp.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List matrixgee_cpp_upd(NumericMatrix data,
                       SEXP covariates,
                       CharacterVector intercepts,
                       int sample_size,
                       int rows_no,
                       int cols_no,
                       int max_iter,
                       double tol,
                       std::string corstr_rows,
                       std::string corstr_cols) {

        DataFrame long_format_dataset = matrix_to_long_dataframe_cpp(data, sample_size, rows_no, cols_no);
        List covariates_fct = covariates_upd(covariates, sample_size, rows_no, cols_no);
        arma::mat covariates_block = covariates_fct["design_covblock"];
        int nrows = covariates_block.n_rows;
        int ncols = covariates_block.n_cols;

        List covariates_df(ncols);
        for (int j = 0; j < ncols; ++j) {
                NumericVector col(nrows);
                for (int i = 0; i < nrows; ++i) {
                        col[i] = covariates_block(i, j);
                }
                covariates_df[j] = col;
        }
        CharacterVector covariate_names(ncols);
        for (int j = 0; j < ncols; ++j) {
                covariate_names[j] = "cov_" + std::to_string(j + 1);
        }
        covariates_df.attr("names") = covariate_names;
        covariates_df.attr("class") = "data.frame";
        covariates_df.attr("row.names") = seq(1, nrows);

        DataFrame covariates_block_df = as<DataFrame>(covariates_df);
        long_format_dataset = Rcpp::as<DataFrame>(Rcpp::Function("cbind")(long_format_dataset, covariates_block_df));

        arma::mat model_matrix = create_design_matrix_cpp(sample_size, rows_no, cols_no,
                                                          intercepts, covariates_block);


        int parameters_no = model_matrix.n_cols;

        NumericVector id_num = long_format_dataset["id"];
        NumericVector rowcol_index_num = long_format_dataset["rowcol_index"];
        arma::vec id = as<arma::vec>(id_num);
        arma::vec rowcol_index = as<arma::vec>(rowcol_index_num);
        NumericVector id_check = long_format_dataset["id"];
        Rcout << "ID min: " << min(id_check) << ", max: " << max(id_check) << std::endl;
        bool converged = false;
        int iteration = 0;
        arma::vec y = as<arma::vec>(long_format_dataset["y"]);

        arma::vec offset = arma::zeros(nrows);
        arma::vec weights = arma::ones(nrows);
        arma::vec beta_start = arma::zeros(parameters_no);
        arma::vec alpha_vector = arma::ones(1);
        double phi = 1.0;
        int mdependence = 1;

        Rcpp::Rcout << "First 5 rows of model_matrix:\n";
        model_matrix.rows(0,4).print();

        Rcpp::Rcout << "First 5 entries of y:\n";
        y.subvec(0, 4).t().print();

        Rcpp::List independence_model = fit_geesolver_cc(
                y, model_matrix, id, rowcol_index, weights, "identity", "gaussian", beta_start,
                offset, 25, 1e-6, 10, 1, 0.0, "gee", 0, alpha_vector, 1, "independence",
                mdependence, phi, 1);

        arma::vec estimate_a = independence_model["beta_hat"];
        arma::vec residuals = independence_model["residuals"];
        double dispersion_param = independence_model["phi"];

        Rcpp::NumericVector residuals_vec(residuals.begin(), residuals.end());
        Rcout << "residuals size: " << residuals.size() << std::endl;

        NumericMatrix correlation_rows = correlation_rows_calc_cpp(
                residuals_vec, sample_size, rows_no, cols_no, parameters_no,
                corstr_rows, dispersion_param);

        NumericMatrix correlation_cols = correlation_cols_calc_cpp(
                residuals_vec, sample_size, rows_no, cols_no, parameters_no,
                corstr_cols, dispersion_param);

        arma::mat correlation_cols_mat(correlation_cols.begin(), correlation_cols.nrow(), correlation_cols.ncol(), false);
        arma::mat correlation_rows_mat(correlation_rows.begin(), correlation_rows.nrow(), correlation_rows.ncol(), false);

        arma::mat kronecker_matrix = arma::kron(correlation_cols_mat, correlation_rows_mat);
        Rcout << "Kronecker matrix size: " << kronecker_matrix.n_rows << " x " << kronecker_matrix.n_cols << std::endl;
        arma::vec upper_triangular_vec = extract_upper_from_matrix(kronecker_matrix);

        Rcpp::List fit_result;
        while (!converged && iteration < max_iter) {
                iteration += 1;
                fit_result = fit_geesolver_cc(
                        y, model_matrix, id, rowcol_index, weights, "identity", "gaussian", beta_start,
                        offset, 25, 1e-6, 10, 1, 0.0, "gee", 0, upper_triangular_vec, 1,
                        "fixed", mdependence, phi, 1);

                residuals = Rcpp::as<arma::vec>(fit_result["residuals"]);
                estimate_a = Rcpp::as<arma::vec>(fit_result["beta_hat"]);
                dispersion_param = fit_result["phi"];
        }

        arma::mat correlation_rows_new = as<arma::mat>(
                correlation_rows_updated_cpp(Rcpp::NumericVector(residuals.begin(), residuals.end()),
                                             sample_size, rows_no, cols_no, parameters_no,
                                             corstr_rows, dispersion_param));

        arma::mat correlation_cols_new = as<arma::mat>(
                correlation_cols_updated_cpp(Rcpp::NumericVector(residuals.begin(), residuals.end()),
                                             sample_size, rows_no, cols_no, parameters_no,
                                             corstr_cols, dispersion_param));

        kronecker_matrix = arma::kron(correlation_cols_new, correlation_rows_new);
        upper_triangular_vec = extract_upper_from_matrix(kronecker_matrix);

        arma::vec new_a = Rcpp::as<arma::vec>(fit_result["beta_hat"]);
        if (new_a.n_elem > 0) {
                double diff = arma::max(arma::abs(new_a - estimate_a));
                if (diff < tol) {
                        converged = true;
                }
        } else {
                Rcpp::stop("Error: null coefficient");
        }

        estimate_a = new_a;
        return fit_result;
}
