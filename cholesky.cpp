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
    // std::cout << "Starting from col " << start<<std::endl;
    double max_val = A(start, start);
    // find the max val 
    // compute elements below the idagonal rows below start
    for (int i=start+1; i<n ; i++) {
        double diagonal = A(i, start);
        if (diagonal > max_val) {
            max_val = diagonal; 
            pivot_row = i; 
            // row index of maximum val in a column

        }
        
    }
    // std::cout << "max val: "<< max_val << std::endl;
    // std::cout << "pivot row: " << pivot_row << std::endl;
    return pivot_row;

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
                double diagonal = A(i,j) - sum;
                // std::cout << "diagonal: " << diagonal << std::endl;
                
                if (diagonal <= 1e-10) {
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
    
    arma::mat Lt = arma::trans(L); // U upper triangular
    // tolerance?
    arma::mat reconstructed = L*Lt;
    // reconstructed.print("reconstructed A"); 
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        std::cout << "Cholesky successful." << std::endl; 
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L,Lt};
    
}

double chol_err(arma::mat & A, bool pivot) {
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
    arma::mat L = result.first; 
    arma::mat Lt = result.second; 
    arma::mat reconstructedA = L*Lt; 
    arma::mat diff = A - reconstructedA; 
    double error = arma::norm(diff, "fro"); 
    return error; 

}

// Ax=b
arma::vec solve(arma::mat & A, arma::vec& b, bool pivot) {
    
    
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
    arma::mat L = result.first; 
    arma::mat Lt = result.second; 
    
    arma::vec x = arma::solve(trimatl(L), arma::solve(trimatu(Lt), b));
    return x; 

}

double calc_diff(const arma::mat& A, const arma::vec& b, const arma::vec&x) {
    arma::vec residual = A*x -b; 
    double r_norm = arma::norm(residual); 
    return r_norm; 

}


std::pair<arma::mat, arma::mat> other_chol(arma::mat& A, bool pivot) {
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
   

        if (pivot_row != j) {

            A.swap_rows(j, pivot_row); 
            P.swap_rows(j,pivot_row);
     
        }
   
        for (int i = 0; i<=j; i++) {
            double sum = 0.0;
            if (i>j) {L(i,j) = 0.0; }
            else if (i == j) {
                for (int k = 0; k<i; k++) 
                    sum += L(i,k) * L(j, k);

                L(i,i) = std::sqrt(A(i,i)) - sum; 
         
            } 
            else if (i<j) {
                for (int k=0; k < i; k++) 
                    sum+= L(j,k)*L(j,k);
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



