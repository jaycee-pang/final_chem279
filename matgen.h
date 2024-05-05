#include <iostream> 
#include <random>
#include <armadillo> 
#include <cmath> 
#include <exception>

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


void make_symmetric(arma::Mat<double> & A);

void make_sympd(arma::Mat<double>& A); 

arma::Mat<double> gen_symmetric(int n);

arma::Mat<double> gen_sympd(int n);

arma::Mat<double> gen_singular(int n); 
