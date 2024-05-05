#include <iostream> 
#include <fstream> 
#include <armadillo> 
#include "matgen.h"
#include "cholesky.h"

#include <chrono> 
#include <random> 

void chol_testing(arma::Mat<double> & A, bool pivot);

double chol_err(arma::Mat<double> & A, bool pivot); 
arma::vec forward_sub(arma::Mat<double> & L, arma::vec&b);
arma::vec backward_sub(arma::Mat<double> & L, arma::vec&y); 
arma::vec solve_lin(arma::Mat<double> & A, arma::vec& b, bool pivot); 
double calc_diff(arma::Mat<double> & A,  arma::vec& b,  arma::vec& x) ;
double inverse_testing(arma::Mat<double> & A, bool pivot); 

double chol_timing(arma::Mat<double> & A, bool pivot) ;

double LU_timing(arma::Mat<double> & A, bool pivot);

double arma_timing(arma::Mat<double> & A);
double LU_decomp_err(arma::Mat<double> & A, bool pivot);