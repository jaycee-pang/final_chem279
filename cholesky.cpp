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

#include "cholesky.h"


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
    int pivot_col = start; 
    std::cout << "Starting from col " << start<<std::endl;
    // goes through rows below current row, starts with column index start
    double max_val = A(start, start);
    // find the max val 
    // compute elements below the idagonal rows below start
    for (int i=start+1; i<n ; i++) {
        double diagonal = A(i,i);
        if (diagonal > max_val) {
            max_val = diagonal;
            pivot_col = i; 

        }
        
    }
     
    std::cout << "max val: "<< max_val << std::endl;
    std::cout << "pivot col: " << pivot_col << std::endl;
    return pivot_col;

}
void permute_cols(arma::mat & A,int j, int pivot_col) {

    if (j < 0 || j >= A.n_cols ||pivot_col < 0 ||pivot_col >= A.n_cols) {
        throw std::invalid_argument("Column indices are out of bounds.");
    }
    
    for (int i = 0; i < A.n_rows; ++i) {
        double temp = A(i, j);
        A(i,j) = A(i,pivot_col);
        A(i,j) = temp;
    }
}

/*
A  = LL.t 
diagonals: Lkk= sqrt(Akk- sum from j=1 to k-1 of [Lkj^2])
off-diagonals for i>k : Lik = (1/Lkk)(Aik - sum from j=1 to k-1 of [Lij*Lkj])
return lower cholesky factor 


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

/*
pivoted: 
(P.t)AP=L(L.t)
select the column with the largest diagonal element as the pivot col 
permute columns of hthe matrix to make this pivot column the current one 


*/
std::pair<arma::mat, arma::mat> pivoted_cholesky(arma::mat& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    // arma::mat R(n,n); 
    L.zeros();
    // R.zeros(); 
    arma::mat P(n,n, arma::fill::eye);  
    int pivot_col; 
    // columns j 
    for (int j=0; j<n; j++) {
        if (pivot) {
            pivot_col = find_pivot(A, j); 

        }

        else {pivot_col = j; }
        std::cout << "Before swapping:" << std::endl;
        std::cout << "A:" << std::endl;
        A.print();
        std::cout << "P:" << std::endl;
        P.print();
        if (pivot_col != j) {
            // A.swap_cols(j,pivot_col);
            // P.swap_cols(j,pivot_col);
            permute_cols(A, j, pivot_col); 
            permute_cols(P, j, pivot_col);
            for (int i = 0; i < P.n_cols; ++i) {
                double temp = P(j, i);
                P(j, i) = P(pivot_col, i);
                P(pivot_col, i) = temp;
            }
            // L.swap_cols(j, pivot_col);
     
        }
        else {
            arma::vec temp = P.col(j); // Store the column to be swapped
            P.shed_col(j);
            P.insert_cols(j, temp);

        }
        std::cout << "After swapping:" << std::endl;
        std::cout << "A:" << std::endl;
        A.print();
        std::cout << "P:" << std::endl;
        P.print();
        // rows i 
        for (int i = 0; i<=j; i++) {
            double sum = 0.0;
            
            for (int k = 0; k<i; k++) {
                sum += L(i,k) * L(j, k);
            }

            if (i == j) {
                // diagonal 
                double diagonal = A(i,j) - sum;
                
                if (diagonal <= 1e-10) {
                    // throw std::runtime_error("Numerical instability.");
                    std::cout << "Warning. Numerical instability."<<std::endl;
                    diagonal+= 1e-10; 
                }
   
                L(i,j) = std::sqrt(diagonal);
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }

    arma::mat Lt = arma::trans(L);
    arma::mat reconstructed = L*Lt;
    // arma::mat reconstructed = P.t()*L*L.t()*P;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    // L = P*L; 

    return {L, Lt};

}

