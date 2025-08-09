//#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link.h"
#include "variance.h"
#include "nuissance.h"
#include "utils.h"

//============================ estimating equations - cc =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec estimating_equations_gee_cc(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const arma::vec & weights_vector,
                                      const char* link,
                                      const char* family,
                                      const arma::vec & beta_vector,
                                      const arma::vec & mu_vector,
                                      const arma::vec & eta_vector,
                                      const char * correlation_structure,
                                      const arma::vec & alpha_vector,
                                      const double & phi) {
        int params_no = model_matrix.n_cols;
        int sample_size = max(id_vector);
        int repeated_max = max(repeated_vector);
        arma::mat ans = arma::zeros(params_no);
        arma::mat correlation_matrix =
                get_correlation_matrix(correlation_structure,
                                       alpha_vector,
                                       repeated_max);
        arma::vec delta_vector = mueta(link, eta_vector);
        arma::vec s_vector = y_vector - mu_vector;
        for(int i=1; i < sample_size + 1; i++){
                arma::uvec id_vector_i = find(id_vector == i);
                arma::mat d_matrix_i =
                        arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
                arma::mat v_matrix_inverse_s_matrix_i =
                        solve(
                                get_v_matrix_cc(family,
                                                mu_vector(id_vector_i),
                                                repeated_vector(id_vector_i),
                                                phi,
                                                correlation_matrix,
                                                weights_vector(id_vector_i)),
                                                s_vector(id_vector_i));
                ans +=
                        trans(d_matrix_i) * v_matrix_inverse_s_matrix_i;
        }
        return ans;
}
//==============================================================================



