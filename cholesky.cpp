// #include <iostream>
// #include <cmath> 
// #include <armadillo> 
// /*
// .is_trimatu / .is_trimatl	 	check whether matrix is upper/lower triangular
// .is_diagmat	 	check whether matrix is diagonal
// .is_square	 	check whether matrix is square sized
// .is_symmetric	 	check whether matrix is symmetric
// .is_hermitian	 	check whether matrix is hermitian
// .is_sympd	 	check whether matrix is symmetric/hermitian positive definite
// */

#include <iostream>
#include <cmath> 
#include <armadillo> 
#include <stdexcept>
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
/* Following accuracy and stability of numerical algorithms 
for j=1:n // cols 
    for i=1:j-1  // rows 
        rij = (aij - sum from k=1 to i-1[rki*rkj])/rii
    rjj = (ajj - sum from k=1 to j-1[rkj^2])^1/2 
*/
std::pair<arma::mat, arma::mat> cholesky(arma::mat& A) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n =A.n_rows; 
    arma::mat L(n,n);
    L.zeros();
    arma::uvec p(n);
    for (int j=0; j<n; j++) {
        for (int i = 0; i<=j; i++) {
            double sum = 0.0;

            for (int k = 0; k<i; k++) {
                sum += L(i,k) * L(j, k);
            }

            if (i == j) {
                L(i,j) = std::sqrt(A(i,j)-sum);
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    
     
    // arma::mat Lt = L.t();  // upper triangular 
    arma::mat Lt = arma::trans(L);
    // tolerance?
    arma::mat reconstructed = L*Lt;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L,Lt};
    
}


std::pair<arma::mat, arma::mat> pivoted_cholesky(arma::mat& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    L.zeros();
    arma::mat P(n,n, arma::fill::eye); 
    int pivot_row; 
    

    for (int j=0; j<n; j++) {
        if (pivot) {
            pivot_row = find_pivot(A, j); 

        }
        else {pivot_row = j; }
        if (pivot_row != j) {
            // use arma::swap rows 
            A.swap_rows(pivot_row, j);
            P.swap_rows(pivot_row, j);
            // p(j) = pivot_row;
        }
        for (int i = 0; i<=j; i++) {
            double sum = 0.0;
            
    
            for (int k = 0; k<i; k++) {
                sum += L(i,k) * L(j, k);
            }

            if (i == j) {
                L(i,j) = std::sqrt(A(i,j)-sum);
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    arma::mat Lt = arma::trans(L);
    arma::mat reconstructed = L*Lt;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L,Lt};

}

std::pair<arma::mat, arma::mat> full_pivoted_cholesky(arma::mat& A) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    L.zeros();
    arma::mat P(n,n, arma::fill::eye); 
    arma::mat Pt(n,n, arma::fill::eye); 

    for (int j = 0; j < n; j++) {
        auto max_it = std::max_element(A.begin_col(j) + j, A.end_col(j));
        int pivot_row = std::distance(A.begin_col(j), max_it);
        int pivot_col = j + std::distance(A.begin_col(j) + j, max_it);
        

        A.swap_rows(j, pivot_row);
        A.swap_cols(j, pivot_col);
        P.swap_rows(j, pivot_row);
        Pt.swap_cols(j, pivot_col);

        for (int i = 0; i<=j; i++) {
            double sum = 0.0;
            
    
            for (int k = 0; k<i; k++) {
                sum += L(i,k) * L(j, k);
            }

            if (i == j) {
                L(i,j) = std::sqrt(A(i,j)-sum);
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    arma::mat Lt = arma::trans(L);
    arma::mat reconstructed = L*Lt;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L,Lt};
}