
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


void make_symmetric(arma::mat & A) {
    
    A = 0.5*(A+A.t()); // symmetric 
    if (!A.is_symmetric()|| A.n_rows != A.n_cols)  {
        throw std::logic_error("Matrix could not be made symmetric"); 
    }
}

void make_sympd(arma::mat& A) {
    A = 0.5*(A*A.t()); // smmetric
    if (!A.is_symmetric() || A.n_rows != A.n_cols)  {
        throw std::logic_error("Matrix could not be made symmetric"); 
    }
    arma::mat sympd = A*A.t();
}

arma::mat gen_symmetric(int n) {
    arma::mat A = arma::randu(n,n);
    arma::mat symmetric = arma::symmatu(A); 
    return symmetric; 

}

arma::mat gen_sympd(int n) {
    // arma::mat A = arma::randu(n,n); 
    // arma::mat sympd = 0.5*(A*A.t());
    // A += 0.1 * arma::eye(n, n); // diagonal vals 
    // // arma::mat U = arma::trimatu(A); 
    // // arma::mat sympd = U*U.t(); 
    // return sympd; 
    // arma::arma_rng::set_seed_random();
    arma::mat A(n, n, arma::fill::randu); 
    arma::mat symm = 0.5 * (A + A.t()); 
    arma::mat pd = symm + 1e-6 * arma::eye(n, n); 
    pd = 0.5*(pd*pd.t());
    return pd;
}


arma::mat gen_singular(int n) {
    arma::mat A(n,n,arma::fill::randu); 
    A = A*A.t();
    arma::vec eigval; 
    arma::mat eigvec; 
    arma::eig_sym(eigval, eigvec, A); 
    eigval*=1e-10; 
    return A; 
}