// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Helper to check if intercept type is included
bool contains(CharacterVector intercepts, std::string key) {
        for (int i = 0; i < intercepts.size(); ++i) {
                if (as<std::string>(intercepts[i]) == key) return true;
        }
        return false;
}

// Scalar intercept block
arma::mat scalar_intercept_block(int N, int rc, bool drop) {
        if (drop) return arma::mat(N * rc, 0, fill::zeros);
        return arma::ones<arma::mat>(N * rc, 1);
}

// Row intercept block with optional drop of one row
arma::mat row_intercept_block(int N, int r, int c, bool drop_last) {
        arma::mat I_r = arma::eye(r, r);
        if (drop_last) I_r.shed_row(r - 1);
        arma::mat block = kron(I_r, arma::ones<mat>(c, 1)); // size rc x (r or r-1)
        arma::mat out(N * r * c, block.n_cols);
        for (int i = 0; i < N; ++i) {
                out.rows(i * r * c, (i + 1) * r * c - 1) = block;
        }
        return out;
}

// Column intercept block with optional drop of one column
arma::mat col_intercept_block(int N, int r, int c, bool drop_last) {
        arma::mat I_c = arma::eye(c, c);
        if (drop_last) I_c.shed_row(c - 1);
        arma::mat block = kron(arma::ones<mat>(r, 1), I_c); // size rc x (c or c-1)
        arma::mat out(N * r * c, block.n_cols);
        for (int i = 0; i < N; ++i) {
                out.rows(i * r * c, (i + 1) * r * c - 1) = block;
        }
        return out;
}

// [[Rcpp::export]]
arma::mat create_design_matrix_cpp(int N, int r, int c,
                                   CharacterVector intercepts,
                                   arma::mat covariates_block) {
        std::vector<arma::mat> blocks;
        int rc = r * c;

        // Determine drops based on combinations
        bool has_row = contains(intercepts, "row");
        bool has_col = contains(intercepts, "col");
        bool has_scalar = contains(intercepts, "scalar");

        bool drop_row = has_row && has_col;   // if both: drop 1 row
        bool drop_col = has_row && has_col;   // if both: drop 1 col
        bool drop_scalar = has_scalar && (has_row || has_col);

        if (has_scalar)
                blocks.push_back(scalar_intercept_block(N, rc, drop_scalar));

        if (has_row)
                blocks.push_back(row_intercept_block(N, r, c, drop_row));

        if (has_col)
                blocks.push_back(col_intercept_block(N, r, c, drop_col));

        if (covariates_block.n_cols > 0)
                blocks.push_back(covariates_block);

        // Combine all blocks
        arma::mat X = blocks[0];
        for (size_t i = 1; i < blocks.size(); ++i) {
                if (blocks[i].n_rows != X.n_rows) {
                        Rcpp::stop("Block %d has wrong number of rows: %d vs %d", i, blocks[i].n_rows, X.n_rows);
                }
                X = join_horiz(X, blocks[i]);
        }

        return X;
}
