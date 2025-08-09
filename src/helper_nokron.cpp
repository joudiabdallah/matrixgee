// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;





// [[Rcpp::export]]
DataFrame vector_to_long_dataframe_cpp(NumericMatrix data, int sample_size, int rows_no, int cols_no) {
        //total_cells is the total number of elements inside the matrix of data
        int total_cells = sample_size * rows_no * cols_no;
        // Define all factors that will be needed
        NumericVector y(total_cells);
        IntegerVector id(total_cells);
        IntegerVector row_index(total_cells);
        IntegerVector col_index(total_cells);
        IntegerVector rowcol_index(total_cells);

        // Fill y with the matrix data values: it works col by col and assume
        //// that each matrix is rxc matrix and will be stacked vertically
        int pos = 0;
        for (int j = 0; j < data.ncol(); j++) {
                for (int i = 0; i < data.nrow(); i++) {
                        y[pos++] = data(i, j);
                }
        }


        for (int i = 0; i < sample_size; i++) {
                for (int c = 0; c < cols_no; c++) {
                        for (int r = 0; r < rows_no; r++) {
                                int index = i * (rows_no * cols_no) + c * rows_no + r;
                                id[index] = i + 1;
                                row_index[index] = r + 1;
                                col_index[index] = c + 1;
                                rowcol_index[index] = c * rows_no + r + 1;
                        }
                }
        }

        return DataFrame::create(
                Named("y") = y,
                Named("id") = id,
                Named("row_index") = row_index,
                Named("col_index") = col_index,
                Named("rowcol_index") = rowcol_index
        );
}


