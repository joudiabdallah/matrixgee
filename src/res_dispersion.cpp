#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List residuals_dispersion_helper(NumericVector model_residuals,
                                 int sample_size,
                                 int rows_no,
                                 int cols_no,
                                 int parameters_no) {


        double residuals_sum_squares = sum(pow(model_residuals, 2.0));


        int a = sample_size * rows_no * cols_no - parameters_no;


        double dispersion_parameter = residuals_sum_squares / a;

        return List::create(
                Named("residuals") = model_residuals,
                Named("dispersion_param") = dispersion_parameter
        );
}
