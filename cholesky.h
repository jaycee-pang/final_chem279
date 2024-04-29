#include <iostream>
#include <cmath> 
#include <armadillo> 
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
int find_pivot(const arma::mat& A, int start) {
    int n = A.n_rows;
    int pivot_row = start; 
    double max_val = std::abs(A(start, start)); 
    // find the max val 
    // compute elements below the idagonal 
    for (int i=start+1; i<n ; i++) {
        if (std::abs(A(i,start)) > max_val) {
            max_val = std::abs(A(i, start));
            pivot_row = i; 

        }
        
    }
    return pivot_row;
}

/*
A  = LL.t 
diagonals: Lkk= sqrt(Akk- sum from j=1 to k-1 of [Lkj^2])
off-diagonals for i>k : Lik = (1/Lkk)(Aik - sum from j=1 to k-1 of [Lij*Lkj])
return lower cholesky factor 
pivoted: 
(P.t)AP=L(L.t)

L is lower triangular
*/


bool pivoted_cholesky(arma::mat& A, arma::mat& L, bool pivot) {
    int n = A.n_rows;
    L.zeros(); 
    arma::uvec p(n); 
    p.zeros();

    for (int i = 0; i < n; i++) {
        if (pivot) {
            int pivot_row = find_pivot(A, i); 
            if (pivot_row != i) {
                // use arma::swap rows 
                A.swap_rows(i,pivot_row);
                L.swap_rows(i,pivot_row);
                p(i) = pivot_row;
            }

        }
        else {
            p(i) = i;
        }
        
        for (int j = i; j<n; j++) {
            double sum = 0.0; 
            for (int k = 0; k<i; k++) {
                sum += L(j,k)*L(i, k);
            }
            if (i == j) { 
                // check for positive definiteness 
                if (A(i,i) - sum <= 0) {return false; }
                L(i,j) = std::sqrt(A(i,i)-sum);
            } else {
                if (L(j,j) == 0) {return false; } // Zero division
                L(j,i) = (1.0 /L(j,j)) * (A(j,i) - sum);
            }
        }
    }
    return true;
}