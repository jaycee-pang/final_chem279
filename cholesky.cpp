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
int find_pivot(const arma::mat& A, int start) {
    int n = A.n_rows;
    int pivot_row = start; 
    // std::cout << "Starting from col " << start<<std::endl;
    double max_val = std::abs(A(start, start));
    // find the max val 
    // compute elements below the idagonal rows below start
    for (int i=start+1; i<n ; i++) {
        double val = A(i, start);
        if (std::abs(val) > max_val) {
            max_val = val; 
            pivot_row = i; 
            // row index of maximum val in a column

        }
        
    }
    std::cout << "max val: "<< max_val << std::endl;
    std::cout << "pivot row: " << pivot_row << std::endl;
    return pivot_row;

}


std::pair<arma::mat, arma::mat> pivoted_cholesky(arma::mat& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    L.zeros();
    arma::mat P(n,n, arma::fill::eye);  
    int pivot_col; 
    int pivot_row; 
    // columns j 
    for (int j=0; j<n; j++) {
        if (pivot) {
            pivot_row = find_pivot(A, j);
        }

        else {pivot_row = j; }
        // std::cout << "Before swapping:" << std::endl;
        // std::cout << "A:" << std::endl;
        // A.print();
        // std::cout << "P:" << std::endl;
        // P.print();

        if (pivot_row != j) {

            A.swap_rows(j, pivot_row); 
            P.swap_rows(j,pivot_row);
     
        }
    
   
        // std::cout << "After swapping:" << std::endl;
        // std::cout << "A:" << std::endl;
        // A.print();
        // std::cout << "P:" << std::endl;
        // P.print();
        // rows i 
        for (int i = 0; i<=j; i++) {
            double sum = 0.0;
            
            for (int k = 0; k<i; k++) {
                sum += L(i,k) * L(j, k);
            }


            if (i == j) {
                // diagonal 
                double diagonal = A(j,j) - sum;
                // std::cout << "diagonal: " << diagonal << std::endl;
                diagonal = A(j,j) - sum + 1e-10; 
                
                if (diagonal <= 1e-10) {
                    std::cout << "Warning. Numerical instability."<<std::endl;
                    diagonal+= 1e-10; 
                    
                    // L(j,j) = std::sqrt(diagonal);
                }
   
                L(j,j) = std::sqrt(std::abs(diagonal));
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    L = P*L; // apply permutations we kept track of in P
    arma::mat Lt = arma::trans(L); // U is upper traingular 
    arma::mat reconstructed = L*Lt;
    // reconstructed.print("reconstructed A"); 
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    
    return {L, Lt};

}



std::pair<arma::mat, arma::mat> other_chol(arma::mat& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::mat L(n,n);
    L.zeros();
    arma::mat P(n,n, arma::fill::eye); // row swaps 
    for (int j = 0; j < n; ++j) {
        double pivot = A(j,j);
        int max_idx = j;
        for (int i = j + 1; i < n; ++i) {
            if (std::abs(A(i,j)) > std::abs(A(max_idx, j))) {
                max_idx = i;
            }
        }
        if (max_idx != j) {
            A.swap_rows(j, max_idx);
            L.swap_rows(j, max_idx);
            P.swap_rows(j, max_idx);
        }

        // Update pivot after potential swap
        pivot = A(j, j);
        std::cout << "pivot: " << pivot << std::endl;

        L(j,j) = std::sqrt(std::abs(pivot));
        for (int i = j+1; i<n; i++) {
            L(i,j) =A(i,j)/L(j,j);
            for (int k = j+1; k <= i; k++) {
                A(i, k) -= L(i,j)*L(k,j);
            }
        }
    }
    L = P*L;
    arma::mat Lt = arma::trans(L); // U upper triangular
    arma::mat reconstructed = L*Lt;
    // if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
    //     std::cout << "Cholesky successful." << std::endl; 
        
    // }
    // else {
    //     std::cout << "Cholesky failed." << std::endl;
    // }
    
    return {L,Lt};

}




std::pair<arma::mat, arma::mat> LU_decomp(arma::mat & A, bool pivot) {
    int n = A.n_rows;
    arma::mat L(n,n,arma::fill::eye); 
    arma::mat U = A; 
    arma::uvec P = arma::regspace<arma::uvec>(0,n-1); 
    int pivot_row; 
    for (int k = 0; k<n-1; k++) {
        if (pivot) {
            pivot_row = find_pivot(U, k);
            U.swap_rows(k, pivot_row); 
            L.swap_rows(k, pivot_row);
            std::swap(P(k), P(pivot_row));
        }
        else {pivot_row = k;}
        
        for (int i = k+1; i<n; i++) {
            double factor = U(i,k) / U(k,k);
            L(i,k) = factor;
            for (int j=k; j<n; j++) {
                U(i,j) -= factor * U(k,j);
            }
        }
    }

    return {L,U};
}


