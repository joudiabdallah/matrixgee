// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "ar1.h"
// [[Rcpp::export]]
List covariates_upd(SEXP covariates, int sample_size, int rows_no, int cols_no) {


        if (Rf_isNull(covariates)) {
                int total_rows = sample_size * rows_no * cols_no;
                arma::mat D_covmat(total_rows, 0);

                return List::create(
                        Named("design_covblock") = D_covmat,
                        Named("paramscov_no") = 0
                );
        }



        else if (Rf_isVector(covariates)) {
                NumericVector cov_vec(covariates);
                int q = cov_vec.size() / sample_size;

                if (cov_vec.size() != sample_size * q) {
                        stop("Length of covariate vector must equal sample_size × q");
                }

                int total_rows = sample_size * rows_no * cols_no;
                arma::mat D_covmat(total_rows, q, fill::zeros);

                List X_subject(sample_size);

                for (int i = 0; i < sample_size; i++) {
                        arma::rowvec x_i(q);
                        for (int j = 0; j < q; j++) {
                                x_i(j) = cov_vec[i * q + j];
                        }

                        X_subject[i] = NumericVector(x_i.begin(), x_i.end());

                        int subject_start = i * (rows_no * cols_no);
                        for (int j = 0; j < rows_no * cols_no; j++) {
                                D_covmat.row(subject_start + j) = x_i;
                        }
                }

                return List::create(
                        Named("design_covblock") = D_covmat,
                        Named("X_subject") = X_subject,
                        Named("paramscov_no") = q
                );
        }







        else if (Rf_isMatrix(covariates)) {
                NumericMatrix cov_mat(covariates);
                int p = cov_mat.nrow() / sample_size;
                int q = cov_mat.ncol();
                mat ones_q = ones<mat>(1, q);
                mat ones_c = ones<mat>(cols_no, 1);
                mat i_r = eye<mat>(rows_no, rows_no);
                int total_rows = sample_size * rows_no * cols_no;
                mat D_covmat = zeros<mat>(total_rows, rows_no * p);

                mat cov_arm(cov_mat.begin(), cov_mat.nrow(), cov_mat.ncol(), false);

                for (int i = 0; i < sample_size; i++) {
                        int start_index = i * p;
                        mat x_i = cov_arm.rows(start_index, start_index + p - 1);
                        mat term_1 = ones_c * ones_q * x_i.t();
                        mat kron_term = kron(term_1, i_r);
                        int row_start = i * (rows_no * cols_no);
                        D_covmat.rows(row_start, row_start + kron_term.n_rows - 1) = kron_term;
                }

                return List::create(Named("design_covblock") = D_covmat,
                                    Named("paramscov_no") = p * q);

        }else {
                stop("Input must be a vector or matrix.");
        }

















}
