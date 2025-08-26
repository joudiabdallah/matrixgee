// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List covariates(SEXP covariates,
                    int sample_size,
                    int rows_no,
                    int cols_no,
                    Nullable<IntegerVector> dim_covariates = R_NilValue) {

        const int total_rows = sample_size * rows_no * cols_no;

        // Case 0: NULL -> no covariates
        if (Rf_isNull(covariates)) {
                arma::mat D_covmat(total_rows, 0);
                return List::create(
                        Named("design_covblock") = D_covmat,
                        Named("paramscov_no") = 0
                );
        }

        // Helper to replicate a row vector across r*c rows
        auto replicate_subject_block = [&](const arma::rowvec& x_row) {
                return arma::repmat(x_row, rows_no * cols_no, 1);
        };

        // Case 1: VECTOR input (subject-level covariates)
        if (Rf_isVector(covariates) && !Rf_isMatrix(covariates)) {
                NumericVector cov_vecR(covariates);

                // 1A) Multiple blocks per subject using dim_covariates
                if (dim_covariates.isNotNull()) {
                        IntegerVector dims(dim_covariates.get());
                        if (dims.size() == 0) stop("dim_covariates must have at least one positive entry.");
                        int sum_q = 0;
                        for (int v : dims) {
                                if (v <= 0) stop("All entries of dim_covariates must be positive.");
                                sum_q += v;
                        }
                        if (cov_vecR.size() != sample_size * sum_q) {
                                stop("Length(covariates) must equal sample_size * sum(dim_covariates).");
                        }

                        arma::mat D_covmat(total_rows, sum_q, arma::fill::zeros);
                        arma::vec cov_vec(cov_vecR.begin(), cov_vecR.size(), /*copy_aux_mem=*/false);

                        // block column starts
                        std::vector<int> col_start(dims.size());
                        {
                                int acc = 0;
                                for (int k = 0; k < dims.size(); ++k) {
                                        col_start[k] = acc;
                                        acc += dims[k];
                                }
                        }

                        const int subj_stride = sum_q;
                        const int m = rows_no * cols_no;

                        for (int i = 0; i < sample_size; ++i) {
                                const int subj_off  = i * subj_stride;
                                const int row_start = i * m;

                                int offset = 0;
                                for (int k = 0; k < dims.size(); ++k) {
                                        const int qk = dims[k];

                                        arma::rowvec x_i_k(qk);
                                        for (int t = 0; t < qk; ++t) {
                                                x_i_k(t) = cov_vec[subj_off + offset + t];
                                        }
                                        arma::mat block = replicate_subject_block(x_i_k);

                                        D_covmat(arma::span(row_start, row_start + m - 1),
                                                 arma::span(col_start[k], col_start[k] + qk - 1)) = block;

                                        offset += qk;
                                }
                        }

                        // also return column start indices (0-based) for mapping
                        IntegerVector col_starts(dims.size());
                        for (int k = 0; k < dims.size(); ++k) col_starts[k] = col_start[k];

                        return List::create(
                                Named("design_covblock") = D_covmat,
                                Named("paramscov_no")    = sum_q,
                                Named("block_sizes")     = dims,
                                Named("block_col_starts")= col_starts
                        );
                }

                // 1B) Backward-compatible single-block behavior (no dim_covariates)
                int q = cov_vecR.size() / sample_size;
                if (cov_vecR.size() != sample_size * q) {
                        stop("Length of covariate vector must equal sample_size × q.");
                }

                arma::mat D_covmat(total_rows, q, arma::fill::zeros);
                arma::vec cov_vec(cov_vecR.begin(), cov_vecR.size(), false);

                const int m = rows_no * cols_no;
                for (int i = 0; i < sample_size; ++i) {
                        arma::rowvec x_i(q);
                        for (int j = 0; j < q; ++j) x_i(j) = cov_vec[i * q + j];
                        arma::mat block = replicate_subject_block(x_i);
                        const int row_start = i * m;
                        D_covmat(arma::span(row_start, row_start + m - 1),
                                 arma::span(0, q - 1)) = block;
                }

                return List::create(
                        Named("design_covblock") = D_covmat,
                        Named("paramscov_no")    = q
                );
        }

        // Case 2: MATRIX input (row-varying matrix, unchanged)
        if (Rf_isMatrix(covariates)) {
                NumericMatrix cov_matR(covariates);
                int p = cov_matR.nrow() / sample_size;  // rows per subject
                int q = cov_matR.ncol();                // features per row-unit
                if (p * sample_size != cov_matR.nrow()) {
                        stop("For matrix covariates, nrow must be sample_size * p.");
                }

                arma::mat cov_arm(cov_matR.begin(), cov_matR.nrow(), cov_matR.ncol(), /*copy_aux_mem=*/false);
                arma::mat D_covmat(total_rows, rows_no * p, arma::fill::zeros);

                arma::mat ones_q = arma::ones<arma::mat>(1, q);
                arma::mat ones_c = arma::ones<arma::mat>(cols_no, 1);
                arma::mat I_r    = arma::eye<arma::mat>(rows_no, rows_no);

                const int m = rows_no * cols_no;
                for (int i = 0; i < sample_size; ++i) {
                        int start_index = i * p;
                        arma::mat x_i = cov_arm.rows(start_index, start_index + p - 1); // p x q
                        arma::mat term_1 = ones_c * ones_q * x_i.t();                   // c x p
                        arma::mat kron_term = arma::kron(term_1, I_r);                  // (r*c) x (r*p)

                        int row_start = i * m;
                        D_covmat.rows(row_start, row_start + kron_term.n_rows - 1) = kron_term;
                }

                return List::create(
                        Named("design_covblock") = D_covmat,
                        Named("paramscov_no")    = p * q
                );
        }

        stop("Input must be NULL, a vector, or a matrix.");
}
