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
 * Look down a column to find the maximum value. 
 *
 * @param A: matrix A (matrix to be decomposed)
 * @return pivot row: row index of the maximum value in the start column to be used as the pivot. 
 */
int find_pivot(const arma::mat& A, int start) ;

/*
A  = LL.t 
diagonals: Lkk= sqrt(Akk- sum from j=1 to k-1 of [Lkj^2])
off-diagonals for i>k : Lik = (1/Lkk)(Aik - sum from j=1 to k-1 of [Lij*Lkj])
return lower cholesky factor 
pivoted: 
(P.t)AP=L(L.t)

L is lower triangular
*/
/* Following accuracy and stability of numerical algorithms 
for j=1:n
    for i=1:j-1 
        rij = (aij - sum from k=1 to i-1[rki*rkj])/rii
    rjj = (ajj - sum from k=1 to j-1[rkj^2])^1/2 
*/
std::pair<arma::mat, arma::mat> cholesky(arma::mat& A); 

std::pair<arma::mat, arma::mat> pivoted_cholesky(arma::mat& A, bool pivot); 
std::pair<arma::mat, arma::mat> other_chol(arma::mat& A, bool pivot); 
double chol_err(arma::mat & A, bool pivot); 
arma::vec solve(arma::mat & A, arma::vec& b, bool pivot); 
double calc_diff(const arma::mat& A, const arma::vec& b, const arma::vec&x); 
// void permute_cols(arma::mat & A,int j, int pivot_col);
// void permute_rows(arma::mat & A,int i, int pivot_row); 




