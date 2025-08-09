//#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link.h"
#include "nuissance.h"
using namespace Rcpp;


//============================ covariance matrices -- cc =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List get_covariance_matrices_cc(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const arma::vec & weights_vector,
                                      const char * link,
                                      const char * family,
                                      const arma::vec & mu_vector,
                                      const arma::vec & eta_vector,
                                      const char * correlation_structure,
                                      const arma::vec & alpha_vector,
                                      const double & phi) {
        double sample_size = max(id_vector);
        double params_no = model_matrix.n_cols;
        arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
        arma::mat meat_matrix = arma::zeros(params_no, params_no);
        arma::vec delta_star_vector = mueta(link, eta_vector);
        arma::vec s_vector = y_vector - mu_vector;
        arma::mat correlation_matrix =
                get_correlation_matrix(correlation_structure,
                                       alpha_vector,
                                       max(repeated_vector));
        for(int i=1; i < sample_size + 1; i++){
                arma::uvec id_vector_i = find(id_vector == i);
                arma::mat d_matrix_i =
                        arma::diagmat(delta_star_vector(id_vector_i)) *
                        model_matrix.rows(id_vector_i);
                arma::mat d_matrix_trans_v_matrix_inverse_i =
                        trans(
                                solve(
                                        get_v_matrix_cc(family,
                                                        mu_vector(id_vector_i),
                                                        repeated_vector(id_vector_i),
                                                        phi,
                                                        correlation_matrix,
                                                        weights_vector(id_vector_i)),
                                                        d_matrix_i));
                arma::vec u_vector_i =
                        d_matrix_trans_v_matrix_inverse_i * s_vector(id_vector_i);
                naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
                meat_matrix += u_vector_i * trans(u_vector_i);
        }
        arma::mat naive_matrix = arma::pinv(naive_matrix_inverse);
        arma::mat naive_matrix_meat_matrix = solve(naive_matrix_inverse, meat_matrix);
        arma::mat robust_matrix = solve(naive_matrix_inverse,
                                        trans(naive_matrix_meat_matrix));
        double obs_no_total = model_matrix.n_rows;
        double kappa = ((obs_no_total - 1)/(obs_no_total - params_no)) *
                (sample_size / (sample_size - 1));
        double lambda = params_no / (sample_size - params_no);
        if (lambda > 0.5) lambda = 0.5;
        double ksi = arma::trace(naive_matrix_meat_matrix)/ params_no;
        if (ksi < 1.0) ksi = 1.0;
        arma::mat bc_matrix =
                kappa * robust_matrix + lambda * ksi * naive_matrix;
        Rcpp::List ans;
        ans["naive_covariance"] = naive_matrix;
        ans["robust_covariance"] = robust_matrix;
        ans["bc_covariance"] = bc_matrix;
        return ans;
}
//==============================================================================


