
#include <iostream> 
#include <random>
#include <armadillo> 
#include <cmath> 
#include <exception>
#include "matgen.h"
/*
.is_trimatu / .is_trimatl	 	check whether matrix is upper/lower triangular
.is_diagmat	 	check whether matrix is diagonal
.is_square	 	check whether matrix is square sized
.is_symmetric	 	check whether matrix is symmetric
.is_hermitian	 	check whether matrix is hermitian
.is_sympd	 	check whether matrix is symmetric/hermitian positive definite
*/
// arma has symmatu 
// arma has symmatu / symmatl	 	generate symmetric matrix from given matrix


void make_symmetric(arma::Mat<double> & A) {
    
    A = 0.5*(A+A.t()); // symmetric 
    if (!A.is_symmetric()|| A.n_rows != A.n_cols)  {
        throw std::logic_error("Matrix could not be made symmetric"); 
    }
}

void make_sympd(arma::Mat<double>& A) {
    A = 0.5*(A*A.t()); // smmetric
    if (!A.is_symmetric() || A.n_rows != A.n_cols)  {
        throw std::logic_error("Matrix could not be made symmetric"); 
    }
    arma::Mat<double> sympd = A*A.t();
}

arma::Mat<double> gen_symmetric(int n) {
    arma::Mat<double> A = arma::randu(n,n);
    arma::Mat<double> symmetric = arma::symmatu(A); 
    return symmetric; 

}

arma::Mat<double> gen_sympd(int n) {
    arma::Mat<double> A = arma::randu<arma::mat>(n, n); 
    A = 0.5*(A+A.t());  // symmetric 

    arma::Mat<double> pd = A + n*arma::eye<arma::mat>(n, n); 
    return pd ;
    
}


arma::Mat<double> gen_singular(int n) {
    arma::Mat<double> A(n,n,arma::fill::randu); 
    A = A*A.t();
    arma::vec eigval; 
    arma::Mat<double> eigvec; 
    arma::eig_sym(eigval, eigvec, A);  // eigen decomposition to ensure symmetry 
    eigval*=1e-10; 
    return A; 
}
