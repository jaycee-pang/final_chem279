
#include <iostream> 
#include <random>
#include <armadillo> 
#include <cmath> 
#include <exception>
#include "matgen.h"

/**
 * Make a given matrix symmetric. 
 *
 * 1/2(A+A.t) 
 * @param A: arma::mat A
 * 
 */
void make_symmetric(arma::Mat<double> & A) {
    
    A = 0.5*(A+A.t()); // symmetric 
    if (!A.is_symmetric()|| A.n_rows != A.n_cols)  {
        throw std::logic_error("Matrix could not be made symmetric"); 
    }
}


/**
 * Generate a symmetric matrix.
 *
 * 1/2(A+A.t) 
 * @param n: int, size of matrix (square matrix)
 * @returns symmetric: symmetric matrix of size nxn. 
 */
arma::Mat<double> gen_symmetric(int n) {
    if (n <=0) {
        throw std::invalid_argument("Invalid matrix size requested.");
    }
    arma::Mat<double> A = arma::randu(n,n);
    arma::Mat<double> symmetric = arma::symmatu(A); 
    return symmetric; 

}


/**
 * Generate a symmetric positive definite matrix. 
 *
 * Average with transpose to make it symmetric, then 
 * add multiple of identity matrix to ensure all eigenvalues are positive. 
 * @param n: int, size of matrix (square matrix)
 * @returns pd: symmetric positive definite matrix of size nxn. 
 */
arma::Mat<double> gen_sympd(int n) {
    if (n <=0) {
        throw std::invalid_argument("Invalid matrix size requested.");
    }
    arma::Mat<double> A = arma::randu<arma::mat>(n, n); 
    A = 0.5*(A+A.t());  // symmetric 
    arma::Mat<double> pd = A + n*arma::eye<arma::mat>(n, n); // p-d
    return pd ;
    
}

/**
 * Generate a nearly singular matrix.  
 *
 * Compute eigenvalues and eigenvecs of A and make them small by multiplying by a small value.  
 * @param n: int, size of matrix (square matrix)
 * @returns A: symmetric, nearly singular matrix. 
 */
arma::Mat<double> gen_singular(int n) {
    if (n <=0) {
        throw std::invalid_argument("Invalid matrix size requested.");
    }
    arma::Mat<double> A(n,n,arma::fill::randu); 
    A = A*A.t();
    arma::vec eigval; 
    arma::Mat<double> eigvec; 
    arma::eig_sym(eigval, eigvec, A);  // eigen decomposition to ensure symmetry 
    eigval*=1e-10; 
    return A; 
   
}
