#include <iostream> 
#include <fstream> 
#include <armadillo> 
#include "matgen.h"
#include "cholesky.h"

#include <chrono> 
#include <random> 

void chol_testing(int n, bool pivot);

double chol_err(arma::mat & A, bool pivot); 
arma::vec forward_sub( arma::mat&L, arma::vec&b);
arma::vec backward_sub( arma::mat&L, arma::vec&y); 
arma::vec solve_lin(arma::mat & A, arma::vec& b, bool pivot); 
double calc_diff( arma::mat& A,  arma::vec& b,  arma::vec& x) ;
double inverse_testing(arma::mat & A, bool pivot); 

double chol_timing(int n, bool pivot) ;

double LU_timing(int n, bool pivot);

double arma_timing(int n);
