#include <iostream> 
#include <random>
#include <armadillo> 
#include <cmath> 
#include <exception>

/**
 * Make a given matrix symmetric. 
 *
 * 1/2(A+A.t) 
 * @param A: arma::mat A
 * 
 */
void make_symmetric(arma::Mat<double> & A);

/**
 * Generate a nearly singular matrix.  
 *
 * Compute eigenvalues and eigenvecs of A and make them small by multiplying by a small value.  
 * @param n: int, size of matrix (square matrix)
 * @returns A: symmetric, nearly singular matrix. 
 */
arma::Mat<double> gen_symmetric(int n);

/**
 * Generate a symmetric positive definite matrix. 
 *
 * Average with transpose to make it symmetric, then 
 * add multiple of identity matrix to ensure all eigenvalues are positive. 
 * @param n: int, size of matrix (square matrix)
 * @returns pd: symmetric positive definite matrix of size nxn. 
 */
arma::Mat<double> gen_sympd(int n);

/**
 * Generate a nearly singular matrix.  
 *
 * Compute eigenvalues and eigenvecs of A and make them small by multiplying by a small value.  
 * @param n: int, size of matrix (square matrix)
 * @returns A: symmetric, nearly singular matrix. 
 */
arma::Mat<double> gen_singular(int n); 
