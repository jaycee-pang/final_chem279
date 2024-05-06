#include <iostream> 
#include <fstream> 
#include <armadillo> 
#include "matgen.h"
#include "cholesky.h"
#include <chrono> 
#include <random> 


// Cholesky reconstruction error 
double chol_err(arma::Mat<double> & A, bool pivot); 

// forward substitution 
arma::vec forward_sub(arma::Mat<double> & L, arma::vec&b);
// backward substitution 
arma::vec backward_sub(arma::Mat<double> & L, arma::vec&y); 
// linear solve 
arma::vec solve_lin(arma::Mat<double> & A, arma::vec& b, bool pivot); 

// find error of inverse 
double inverse_testing(arma::Mat<double> & A, bool pivot); 

// time Cholesky 
double chol_timing(arma::Mat<double> & A, bool pivot) ;
// time LU decomposition 
double LU_timing(arma::Mat<double> & A, bool pivot);
// time arma Cholesky  
double arma_timing(arma::Mat<double> & A);
// Find the reconstruction error from LU decomposition. 
double LU_decomp_err(arma::Mat<double> & A, bool pivot);


void chol_testing(arma::Mat<double> & A, bool pivot);
double calc_diff(arma::Mat<double> & A,  arma::vec& b,  arma::vec& x) ;