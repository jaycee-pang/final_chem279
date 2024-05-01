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
    int pivot_row = start; 
    std::cout << "Starting from row" << start<<std::endl;
    // goes through rows below current row, starts with column index start
    // double max_val = std::abs(A(start, start)); // diagonal
    double max_val = A(start, start);
    // find the max val 
    // compute elements below the idagonal rows below start
    // compare values in the column below start (exclude diagonal)
    for (int i=start+1; i<n ; i++) {
        if (A(i,start) > max_val) {
            max_val = A(i, start);
            pivot_row = i; 

        }
        
    }
    std::cout << "max val: "<< max_val << std::endl;
    return pivot_row;

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

*/
std::pair<arma::mat, arma::mat> pivoted_cholesky(arma::mat& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    arma::mat R(n,n); 
    L.zeros();
    R.zeros(); 

    arma::mat P(n,n, arma::fill::eye); 
    int pivot_row; 

    for (int j=0; j<n; j++) {
        if (pivot) {
            pivot_row = find_pivot(A, j); 
            std::cout << "pivot row: " << pivot_row << std::endl;
            

        }
        else {pivot_row = j; }
        if (pivot_row != j) {
            // // use arma::swap rows 
            
            A.swap_rows(j, pivot_row);
            L.swap_cols(j, pivot_row);
            P.swap_rows(j, pivot_row);
            

            
     
        }
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
                // L(i,j) = std::sqrt(A(i,j)-sum);
                L(i,j) = std::sqrt(diagonal);
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    // arma::mat Pt = arma::trans(P);
    // arma::mat reconstructed = L*Lt;
    arma::mat Lt = arma::trans(L);
    // arma::mat reconstructed = L*Lt;
    arma::mat reconstructed = P.t()*L*L.t()*P;
    // reconstructed.print();
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {P, L};

}


std::pair<arma::umat, arma::mat> full_pivoted_cholesky(arma::mat& A, bool temp) {
    int n = A.n_rows;
    arma::mat R(n,n, arma::fill::zeros);
    arma::umat P = arma::eye<arma::umat>(n, n); 

    for (int k = 0; k < n; ++k) {
 
        arma::mat B = A.submat(k, k, n-1, n-1);
        std::cout << "B size submat: " << B.size() << std::endl;
        arma::uword l;
        B.max(l);
        l += k;
        std::cout <<"l: " << l << std::endl;

        // rows and cols 
        A.swap_rows(k, l);
        A.swap_cols(k, l);
        R.swap_cols(k, l);
        P.swap_rows(k, l);


        R(k, k) = sqrt(A(k, k));
        R.submat(k, k+1, n-1, n-1) = A.submat(k, k+1, n-1, n-1) / R(k, k);
        

        A.submat(k+1, k+1, n-1, n-1) -= R.submat(k, k+1, n-1, n-1) * trans(R.submat(k, k+1, n-1, n-1));
    }


    // Reconstruction based on Cholesky factors (uncomment if needed)
    // arma::mat reconstructed = P * (R * trans(R)) * trans(P);
    // if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
    //     std::cout << "Cholesky successful." << std::endl; 
    // } else {
    //     std::cout << "Cholesky failed." << std::endl;
    // }
    
    return {P,R};
}
