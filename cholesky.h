#include <iostream>
#include <cmath> 
#include <armadillo> 
#include <stdexcept>
#include <exception> 
/*
.is_trimatu / .is_trimatl	 	check whether matrix is upper/lower triangular
.is_diagmat	 	check whether matrix is diagonal
.is_square	 	check whether matrix is square sized
.is_symmetric	 	check whether matrix is symmetric
.is_hermitian	 	check whether matrix is hermitian
.is_sympd	 	check whether matrix is symmetric/hermitian positive definite
*/


/**
 * Find the pivot in the column A from the given starting index. 
 *
 * Look below the diagonal
 *
 * @param A: matrix A (matrix to be decomposed)
 * @return pivot row: row index of the maximum value to be used as the pivot. 
 */
int find_pivot(const arma::Mat<double>& A, int start) ;

/*
A  = LL.t 
L is lower triangular
*/
/* Following accuracy and stability of numerical algorithms 
for j=1:n // cols 
    for i=1:j-1  // rows 
        rij = (aij - sum from k=1 to i-1[rki*rkj])/rii
    rjj = (ajj - sum from k=1 to j-1[rkj^2])^1/2 
*/
// std::pair<arma::Mat<double>, arma::Mat<double>> cholesky(arma::mat& A); 

/*
pivoted: 
(P.t)AP=L(L.t)
*/
std::pair<arma::Mat<double>, arma::Mat<double>> pivoted_cholesky(arma::Mat<double>& A, bool pivot); 
std::pair<arma::Mat<double>,arma::Mat<double>> LU_decomp(arma::Mat<double> & A, bool pivot); 

std::pair<arma::Mat<double>,arma::Mat<double>> other_chol(arma::Mat<double>& A, bool pivot); 
std::pair<arma::Mat<double>, arma::Mat<double>> full_pivoted_cholesky(arma::Mat<double>& A); 
std::pair<int,int> find_full_pivot(arma::Mat<double>&A, int start); 
// void permute_cols(arma::mat & A,int j, int pivot_col);
// void permute_rows(arma::mat & A,int i, int pivot_row); 




