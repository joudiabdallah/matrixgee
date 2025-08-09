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
#include "correlation_calc.h"
#include "correlation_updated.h"
#include "helper_data.h"


using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List matrixgee_cpp_nokron(NumericMatrix data,
                          SEXP covariates,
                          std::string intercept,
                          int sample_size,
                          int rows_no,
                          int cols_no,
                          int max_iter,
                          double tol,
                          std::string corstr) {

        DataFrame long_format_dataset = matrix_to_long_dataframe_cpp(data,
                                                                     sample_size,
                                                                     rows_no,
                                                                     cols_no);


        Rcout << "matrix_to_long_dataframe_cpp" << std::endl;

        Rcout << "Number of rows in long_format_dataset: "
              << long_format_dataset.nrows() << std::endl;

        Rcout << "Columns in long_format_dataset: ";
        CharacterVector colnames_1 = long_format_dataset.names();
        for (int i = 0; i < colnames_1.size(); ++i) {
                Rcout << as<std::string>(colnames_1[i]) << " ";
        }
        Rcout << std::endl;
        //Rcout << "matrix_to_long_dataframe_cpp" << std::endl;

        //Rcout << "Number of rows in long_format_dataset: "
        //<< long_format_dataset.nrows() << std::endl;

        //Rcout << "Columns in long_format_dataset: ";
        // CharacterVector colnames_1 = long_format_dataset.names();
        //for (int i = 0; i < colnames_1.size(); ++i) {
        // Rcout << as<std::string>(colnames_1[i]) << " ";
        //}
        //Rcout << std::endl;

        //IntegerVector id = long_format_dataset["id"];
        //IntegerVector rowcol_index = long_format_dataset["rowcol_index"];
        //NumericVector id_num = as<NumericVector>(long_format_dataset["id"]);
        //arma::vec id = as<arma::vec>(id_num);
        //NumericVector rowcol_index_num = as<NumericVector>(long_format_dataset["rowcol_index"]);
        //arma::vec rowcol_index = as<arma::vec>(rowcol_index_num);
        //NumericVector id_num = as<NumericVector>(long_format_dataset["id"]);
        //NumericVector rowcol_index_num = as<NumericVector>(long_format_dataset["rowcol_index"]);
        CharacterVector colnames_final = long_format_dataset.names();
        bool has_y = false;
        for (int i = 0; i < colnames_final.size(); ++i) {
                if (as<std::string>(colnames_final[i]) == "y") {
                        has_y = true;
                }
        }
        if (!has_y) {
                Rcout << "[ERROR] 'y' column missing from final dataset before model.matrix()\n";
        }

        Function cbind("cbind");
        int parameters_no;

        //dataset and id is now ready!!

        //define ones_c and I_r matrices
        arma::mat ones_c = arma::ones<arma::mat>(cols_no, 1);
        arma::mat I_R = arma::eye<arma::mat>(rows_no, rows_no);
        //get the covariates_fct
        //List covariates_fct = covariates_cpp(covariates, sample_size, rows_no, cols_no,intercept);
        List covariates_fct = covariates_upd(covariates, sample_size, rows_no, cols_no);
        // Extract covparams_no
        //int covparams_no = as<int>(covariates_fct["paramscov_no"]);
        // Extract covariates_block
        arma::mat covariates_block = covariates_fct["design_covblock"];
        //get the number of columns in covariates block
        int ncols = covariates_block.n_cols;
        //get the number of rows in covariates block
        int nrows = covariates_block.n_rows;
        //create an empty list covariates_df of length ncols
        List covariates_df(ncols);
        Rcout << "covariates_block size: "
              << covariates_block.n_rows << " x " << covariates_block.n_cols << std::endl;

        for (int j = 0; j < ncols; ++j) {
                NumericVector col(nrows);
                for (int i = 0; i < nrows; ++i) {
                        col[i] = covariates_block(i, j);
                }
                covariates_df[j] = col;
        }

        //create colnames of length ncols
        CharacterVector colnames(ncols);
        for (int j = 0; j < ncols; ++j) {
                colnames[j] = "cov_" + std::to_string(j + 1);
        }

        covariates_df.attr("names") = colnames;
        covariates_df.attr("class") = "data.frame";
        covariates_df.attr("row.names") = seq(1, nrows);

        //turn a list into dataframe
        DataFrame covariates_block_df = as<DataFrame>(covariates_df);
        //long_format_dataset = DataFrame(Rcpp::wrap(Rcpp::List(long_format_dataset).attr("names")));
        long_format_dataset = Rcpp::as<DataFrame>(Rcpp::Function("cbind")(long_format_dataset, covariates_block_df));

        // Call cbind() and store its result in a temporary R object
        //RObject cbind_result = Rcpp::Function("cbind")(long_format_dataset, covariates_block_df);

        // Print the type of the result
        //Rcout << "[DEBUG] cbind_result SEXP type: " << cbind_result.sexp_type() << std::endl;

        // Convert to DataFrame and inspect dimensions
        //DataFrame combined_df = Rcpp::as<DataFrame>(cbind_result);

        //Rcout << "[DEBUG] combined_df row count: " << combined_df.nrows() << std::endl;
        //Rcout << "[DEBUG] combined_df column count: " << combined_df.size() << std::endl;

        // Assign it back to long_format_dataset
        //long_format_dataset = combined_df;





        arma::mat model_matrix;

        if (intercept == "row_intercept") {
                arma::mat ones_I = arma::ones<arma::mat>(sample_size, 1);
                arma::mat ones_C = arma::ones<arma::mat>(cols_no, 1);
                arma::mat I_R = arma::eye(rows_no, rows_no);
                arma::mat intblock = arma::kron(ones_I, arma::kron(ones_C, I_R));
                int ncols_int = intblock.n_cols;
                int nrows_int = intblock.n_rows;
                Rcout << " intblock size "
                      << nrows_int << " x " << ncols_int << std::endl;
                List intblock_df(ncols_int);
                for (int j = 0; j < ncols_int; ++j) {
                        NumericVector col(nrows_int);
                        for (int i = 0; i < nrows_int; ++i) {
                                col[i] = intblock(i, j);
                        }
                        intblock_df[j] = col;
                }

                Rcout << " intblock_df " << ncols_int << " columns and " << nrows_int << " rows." << std::endl;


                CharacterVector rowintercept_names(ncols_int);
                for (int j = 0; j < ncols_int; ++j) {
                        rowintercept_names[j] = "row_index" + std::to_string(j + 1);
                }

                intblock_df.attr("names") = rowintercept_names;
                intblock_df.attr("class") = "data.frame";
                intblock_df.attr("row.names") = seq(1, nrows_int);


                Function as_data_frame("as.data.frame");
                DataFrame intblock_df_final = as_data_frame(intblock_df);


                Function cbind("cbind");
                long_format_dataset = cbind(long_format_dataset, intblock_df_final);
                // after cbind possibly strips column names:



                Function names("names");
                CharacterVector colnames_check = names(long_format_dataset);
                Rcout << " column names in long_format_dataset:\n";
                for (int i = 0; i < colnames_check.size(); ++i) {
                        Rcout << as<std::string>(colnames_check[i]) << ", ";
                }
                Rcout << std::endl;
                CharacterVector covariate_names = covariates_block_df.attr("names");
                std::string formula_str = "y ~ ";
                for (int j = 0; j < rowintercept_names.size(); ++j) {
                        formula_str += std::string(rowintercept_names[j]) + " + ";
                }

                for (int j = 0; j < covariate_names.size(); ++j) {
                        formula_str += std::string(covariate_names[j]) + " + ";
                }

                formula_str = formula_str.substr(0, formula_str.size() - 3) + " - 1";

                Function as_formula("as.formula");
                Function model_matrix_f("model.matrix");
                SEXP formula = as_formula(formula_str);

                Rcout << "formula string: " << formula_str << std::endl;
                Rcout << "Names in long_format_dataset before model.matrix:\n";
                CharacterVector names_before = long_format_dataset.attr("names");
                for (int i = 0; i < names_before.size(); ++i) {
                        Rcout << as<std::string>(names_before[i]) << ", ";
                }
                Rcout << std::endl;
                RObject mm_r = model_matrix_f(formula, long_format_dataset);
                Rcout << "SEXP type of model.matrix output: " << TYPEOF(mm_r) << std::endl;
                model_matrix = Rcpp::as<arma::mat>(model_matrix_f(formula, long_format_dataset));
                Rcout << " model_matrix size: "
                      << model_matrix.n_rows << " x " << model_matrix.n_cols << std::endl;

                Rcpp::Rcout << "Column means of model_matrix:\n";
                Rcpp::Rcout << arma::mean(model_matrix).t() << std::endl;

                arma::mat XtX = model_matrix.t() * model_matrix;

                //parameters_no = model_matrix.n_cols;

        }

        else if (intercept == "column_intercept") {
                arma::mat i_c = arma::eye<arma::mat>(cols_no, cols_no);
                arma::mat ones_r = arma::ones<arma::mat>(rows_no, 1);
                arma::mat ones_I = arma::ones<arma::mat>(sample_size, 1);
                //arma::mat intblock = arma::kron(i_c, ones_r);
                arma::mat intblock = arma::kron(ones_I, arma::kron(i_c, ones_r));

                int ncols = intblock.n_cols;
                int nrows = intblock.n_rows;
                List intblock_df(ncols);

                for (int j = 0; j < ncols; ++j) {
                        NumericVector col(nrows);
                        for (int i = 0; i < nrows; ++i) {
                                col[i] = intblock(i, j);
                        }
                        intblock_df[j] = col;
                }


                CharacterVector colnames(ncols);
                for (int j = 0; j < ncols; ++j) {
                        colnames[j] = "intblock_" + std::to_string(j + 1);
                }

                intblock_df.attr("names") = colnames;
                intblock_df.attr("class") = "data.frame";
                intblock_df.attr("row.names") = seq(1, nrows);


                DataFrame intblock_df_final = as<DataFrame>(intblock_df);
                Function cbind("cbind");
                long_format_dataset = cbind(long_format_dataset, intblock_df_final);


                CharacterVector covariate_names = covariates_block_df.attr("names");
                std::string formula_str = "y ~ ";


                for (int j = 0; j < colnames.size(); ++j) {
                        formula_str += std::string(colnames[j]) + " + ";
                }


                for (int j = 0; j < covariate_names.size(); ++j) {
                        formula_str += std::string(covariate_names[j]) + " + ";
                }


                formula_str = formula_str.substr(0, formula_str.size() - 3) + " - 1";


                Function as_formula("as.formula");
                Function model_matrix_f("model.matrix");
                SEXP formula = as_formula(formula_str);

                model_matrix = Rcpp::as<arma::mat>(model_matrix_f(formula, long_format_dataset));
                //parameters_no = model_matrix.n_cols;

        }


        else if (intercept == "matrix_intercept") {
                arma::mat ones_I = arma::ones<arma::mat>(sample_size, 1);
                arma::mat i_c = arma::eye<arma::mat>(cols_no, cols_no);
                arma::mat i_r = arma::eye<arma::mat>(rows_no, rows_no);
                //arma::mat intblock = arma::kron(i_c, i_r);
                arma::mat intblock = arma::kron(ones_I, arma::kron(i_c, i_r));



                int ncols = intblock.n_cols;
                int nrows = intblock.n_rows;
                List intblock_df(ncols);


                for (int j = 0; j < ncols; ++j) {
                        NumericVector col(nrows);
                        for (int i = 0; i < nrows; ++i) {
                                col[i] = intblock(i, j);
                        }
                        intblock_df[j] = col;
                }


                CharacterVector colnames(ncols);
                for (int j = 0; j < ncols; ++j) {
                        colnames[j] = "intblock_" + std::to_string(j + 1);
                }

                intblock_df.attr("names") = colnames;
                intblock_df.attr("class") = "data.frame";
                intblock_df.attr("row.names") = seq(1, nrows);


                DataFrame intblock_df_final = as<DataFrame>(intblock_df);
                Function cbind("cbind");
                long_format_dataset = cbind(long_format_dataset, intblock_df_final);


                CharacterVector covariate_names = covariates_block_df.attr("names");
                std::string formula_str = "y ~ ";

                for (int j = 0; j < colnames.size(); ++j) {
                        formula_str += std::string(colnames[j]) + " + ";
                }

                for (int j = 0; j < covariate_names.size(); ++j) {
                        formula_str += std::string(covariate_names[j]) + " + ";
                }


                formula_str = formula_str.substr(0, formula_str.size() - 3) + " - 1";


                Function as_formula("as.formula");
                Function model_matrix_f("model.matrix");
                SEXP formula = as_formula(formula_str);
                model_matrix = Rcpp::as<arma::mat>(model_matrix_f(formula, long_format_dataset));
                //parameters_no = model_matrix.n_cols;



        }


        else if (intercept == "scalar_intercept") {
                arma::mat ones_I = arma::ones<arma::mat>(sample_size, 1);
                arma::mat one_c = arma::ones<arma::mat>(cols_no, 1);
                arma::mat one_r = arma::ones<arma::mat>(rows_no, 1);
                //arma::mat intblock = arma::kron(one_c, one_r);
                arma::mat intblock = arma::kron(ones_I, arma::kron(one_c, one_r));



                int ncols = intblock.n_cols;
                int nrows = intblock.n_rows;
                List intblock_df(ncols);

                for (int j = 0; j < ncols; ++j) {
                        NumericVector col(nrows);
                        for (int i = 0; i < nrows; ++i) {
                                col[i] = intblock(i, j);
                        }
                        intblock_df[j] = col;
                }


                CharacterVector colnames(ncols);
                for (int j = 0; j < ncols; ++j) {
                        colnames[j] = "scalar_intercept";
                }

                intblock_df.attr("names") = colnames;
                intblock_df.attr("class") = "data.frame";
                intblock_df.attr("row.names") = seq(1, nrows);


                DataFrame intblock_df_final = as<DataFrame>(intblock_df);
                Function cbind("cbind");
                long_format_dataset = cbind(long_format_dataset, intblock_df_final);


                CharacterVector covariate_names = covariates_block_df.attr("names");
                std::string formula_str = "y ~ scalar_intercept + ";

                for (int j = 0; j < covariate_names.size(); ++j) {
                        formula_str += std::string(covariate_names[j]) + " + ";
                }


                formula_str = formula_str.substr(0, formula_str.size() - 3) + " - 1";


                Function as_formula("as.formula");
                Function model_matrix_f("model.matrix");
                SEXP formula = as_formula(formula_str);
                model_matrix = Rcpp::as<arma::mat>(model_matrix_f(formula, long_format_dataset));
                //parameters_no = model_matrix.n_cols;


        }

        parameters_no = model_matrix.n_cols;
        //We now have the formula_2 and the nb of parameters
        //+ the correct  format of the dataset.///////////

        //Initialize the convergence criteria which will be the same for different
        // conditions, + independent model will be needed.
        // After all intercept and covariate blocks are added
        NumericVector id_num = long_format_dataset["id"];
        NumericVector rowcol_index_num = long_format_dataset["rowcol_index"];
        arma::vec id = as<arma::vec>(id_num);
        arma::vec rowcol_index = as<arma::vec>(rowcol_index_num);
        NumericVector id_check = long_format_dataset["id"];
        Rcout << "ID min: " << min(id_check) << ", max: " << max(id_check) << std::endl;
        bool converged = false;
        int iteration = 0;
        NumericVector y_col = long_format_dataset["y"];
        y_col.attr("class") = R_NilValue;
        long_format_dataset["y"] = y_col;
//arma::vec y = as<arma::vec>(long_format_dataset["y"]);
        arma::vec y = as<arma::vec>(long_format_dataset["y"]);



        //arma::vec id_vector = as<arma::vec>(long_format_dataset["id"]);
        //arma::vec repeated_vector = as<arma::vec>(long_format_dataset["rowcol_index"]);


        CharacterVector df_names = long_format_dataset.names();
        std::vector<std::string> design_colnames;
        for (int i = 0; i < df_names.size(); ++i) {
                std::string name = as<std::string>(df_names[i]);
                if (name.rfind("cov_", 0) == 0 || name.rfind("intblock_", 0) == 0 || name.rfind("row_index", 0) == 0) {
                        design_colnames.push_back(name);
                }
        }

        int N = long_format_dataset.nrows();

        arma::vec offset = arma::zeros(N);
        arma::vec weights = arma::ones(N);
        //arma::vec beta_start = arma::zeros(p);
        //arma::vec beta_start = arma::zeros(parameters_no);
        arma::vec beta_start = arma::zeros(parameters_no);
        arma::vec alpha_vector = arma::ones(1);
        double phi = 1.0;
        int mdependence = 1;

        //arma::vec
        id = Rcpp::as<arma::vec>(Rcpp::wrap(id));
        //arma::vec
        rowcol_index = Rcpp::as<arma::vec>(Rcpp::wrap(rowcol_index));
        //arma::vec
        weights = Rcpp::as<arma::vec>(Rcpp::wrap(weights));
        arma::vec alpha_vec = Rcpp::as<arma::vec>(Rcpp::wrap(alpha_vector));

        Rcpp::Rcout << "First 5 rows of model_matrix:\n";
        model_matrix.rows(0,4).print();
        //Rcpp::Rcout << " Full covariate column:\n";
        //model_matrix.col(rows_no).t().print();
        //Rcpp::Rcout << "Row means of model_matrix covariates:\n";
        //Rcpp::Rcout << arma::mean(model_matrix.cols(10, 10), 1).t() << std::endl;

        Rcpp::Rcout << "First 5 entries of y:\n";
        y.subvec(0, 4).t().print();




        Rcpp::List independence_model = fit_geesolver_cc(
                y,
                model_matrix,
                id,
                rowcol_index,
                weights,
                "identity",
                "gaussian",
                beta_start,
                offset,
                25,
                1e-6,
                10,
                1,
                0.0,
                "gee",
                0,
                alpha_vector,
                1,
                "independence",
                mdependence,
                phi,
                1
        );
        arma::vec estimate_a = independence_model["beta_hat"];
        arma::vec residuals = independence_model["residuals"];
        double dispersion_param = independence_model["phi"];


        //Rcout << "First 5 residuals: " << residuals.head(5).t();
        //Rcout << "First 5 beta_hat: " << estimate_a.head(5).t();
        //Rcout << "Dispersion param: " << dispersion_param << std::endl;


        Rcpp::NumericVector residuals_vec(residuals.begin(), residuals.end());
        Rcout << "residuals size: " << residuals.size() << std::endl;


        NumericMatrix correlation = correlation_general(
                residuals_vec,
                sample_size,
                rows_no,
                cols_no,
                parameters_no,
                corstr,
                dispersion_param
        );




        arma::mat correlation_arma = as<arma::mat>(correlation);
        // extract the upper triangular part of the Kroncker matrix
        arma::vec upper_triangular_vec = extract_upper_from_matrix(correlation_arma);


        Rcpp::List fit_result;
        // start from independence estimates
        //arma::vec beta_start = estimate_a;
        //double     phi        = dispersion_param;

        while (!converged && iteration < max_iter) {
                iteration += 1;

                // 1) Refit with current fixed correlation and current beta_start
                fit_result = fit_geesolver_cc(
                        y,
                        model_matrix,
                        id,
                        rowcol_index,
                        weights,
                        "identity",
                        "gaussian",
                        beta_start,            // <-- start at current beta
                        offset,
                        25,
                        1e-6,
                        10,
                        1,
                        0.0,
                        "gee",
                        0,
                        upper_triangular_vec,  // <-- fixed working correlation vector (upper-tri)
                        1,
                        "fixed",
                        mdependence,
                        phi,
                        1
                );

                // 2) Pull updates from the fit
                arma::vec new_beta    = Rcpp::as<arma::vec>(fit_result["beta_hat"]);
                arma::vec new_resids  = Rcpp::as<arma::vec>(fit_result["residuals"]);
                double    dispersion_param     = fit_result["phi"];

                if (new_beta.n_elem == 0) {
                        Rcpp::stop("Error: null coefficient");
                }

                // 3) Convergence check (max |Δβ|)
                double diff = arma::max(arma::abs(new_beta - estimate_a));
                if (diff < tol) {
                        estimate_a = new_beta;
                        //phi        = new_;
                        converged  = true;
                        break;
                }

                // 4) Not converged: update state for next iteration
                estimate_a = new_beta;
                residuals  = new_resids;
                //phi        = new_phi;
                beta_start = estimate_a;

                // 5) Re-estimate the general correlation FROM updated residuals
                arma::mat correlation_new = as<arma::mat>(
                        correlation_general_updated(
                                Rcpp::NumericVector(new_resids.begin(), new_resids.end()),
                                sample_size,
                                rows_no,
                                cols_no,
                                parameters_no,
                                corstr,
                                dispersion_param
                        )
                );

                // 6) Refresh the fixed-correlation vector for the next fit
                upper_triangular_vec = extract_upper_from_matrix(correlation_new);
        }

        // (optional) annotate result and sync final params
        fit_result["converged"]  = Rcpp::wrap(converged);
        fit_result["iterations"] = Rcpp::wrap(iteration);
        fit_result["beta_hat"]   = Rcpp::wrap(estimate_a);
        fit_result["phi"]        = Rcpp::wrap(phi);

        return fit_result;}







