// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


#include "RcppArmadillo.h"
#include <algorithm>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
// rcppaarma_fill_matrix
// it just fills a matrix with its attributes, adding a new column to it.
// @param matrix a matrix
// @param rowR rowindex in C style
// @param colR colindex in C style
// @param value the double value of Z
arma::mat rcpparma_fill_matrix(arma::mat matrix, arma::uword row,
                      arma::uword col, double value)
{
    matrix.insert_rows(matrix.n_rows, 1);
    matrix(matrix.n_rows-1, 0) = row;
    matrix(matrix.n_rows-1, 1) = col;
    matrix(matrix.n_rows-1, 2) = value;
    return(matrix);
}
//' rcpparma_get_disjoint_max_win
//' Computes the disjoint max_win matrix.
//' @param z0 a matrix.
//' @param sigwin sigwin.
//' @param zthresh zthresh.
//' @param nmax nmax.
//' @param verbose verbose.
//' @return a matrix of three columns (bin_idx, win_idx, z_val) idxs in C style.
//' @keywords internal
// [[Rcpp::export]]
arma::mat  rcpparma_get_disjoint_max_win(arma::mat z0,  int sigwin=10,
                                            double zthresh=10,
                                            arma::uword nmax=9999999,
                                            bool verbose=true)
{

    if(verbose) Rcpp::Rcout << "Maximizing with zthresh: " << zthresh <<
                                "\tsigwin: " << sigwin << "\n" ;

    arma::uword maxwin=z0.n_cols;
    int i=1;
    arma::mat s=arma::mat(0, 3);
    while (TRUE)
    {
        // findind row and col indexes of max value
        arma::umat maxv=arma::index_max(z0);
        arma::mat maxv1=arma::mat(1, maxv.n_elem);

        if(maxv1(1) == -std::numeric_limits<float>::infinity()) {break;}

        for(arma::uword j=0; j<maxv.n_elem; j++)
        {
            maxv1(j) = z0(maxv(j), j);
        }
        arma::vec maxv3 = arma::conv_to< arma::vec >::from(maxv1);
        arma::uword maxv2 = arma::index_max(maxv3);
        arma::uword maxr = maxv(maxv2);
        arma::uword maxc = maxv2;
        // end of maxvalues
        if (z0(maxr, maxc) < zthresh) {break;}
        // Rcpp::Rcout << "maxr: " << maxr << "\tmaxc: " << maxc << "\n" ;
        s=rcpparma_fill_matrix(s, maxr, maxc, z0(maxr, maxc));
        if(verbose) {
             if((i % 100) == 0) Rcpp::Rcout << "Maximizing window: " << maxr <<
                        ",\t" << maxc << "\tScore=" << z0(maxr, maxc) << "\n" ;
        }
        i++;
        int val = (maxr - sigwin - maxwin + 1);
        int st=0;
        if(val > 0)
        {
            st=val;
        }else{
            st=0;
        }
        int ed=std::min(maxr + maxc + sigwin, z0.n_rows-1);
        // Rcpp::Rcout << "st: " << st << "\ted:" << ed << "\n" ;
        for(int j=st; j<=ed; j++)
        {
            for(arma::uword k=0; k<z0.n_cols; k++)
            {
                z0(j,k) = -std::numeric_limits<float>::infinity();
            }
        }
        if (s.n_rows >= nmax) break;
    }
    return(s);
}
