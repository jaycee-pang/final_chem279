#include "cholesky.h"
/**
 * Find the pivot in the column A from the given starting index. 
 *
 * Look below the diagonal
 *
 * @param A: matrix A (matrix to be decomposed)
 * @return pivot row: row index of the maximum value to be used as the pivot. 
 */
int find_pivot(const arma::Mat<double>& A, int start) {
    int n = A.n_rows;
    int pivot_row = start; 
    // std::cout << "Starting from col " << start<<std::endl;
    double max_val = std::abs(A(start, start));
    // find the max val 
    // compute elements below the idagonal rows below start
    for (int i=start+1; i<n ; i++) {
        double val = A(i, start);
        if (std::abs(val) > std::abs(max_val)) {
            max_val = val; 
            pivot_row = i; 
            // row index of maximum val in a column

        }
        
    }
    // std::cout << "max val: "<< max_val << std::endl;
    // std::cout << "pivot row: " << pivot_row << std::endl;
    return pivot_row;

}

/**
 * Cholesky Decomposition of the matrix A
 *
 *
 * @param A: matrix A (matrix to be decomposed)
 * @param pivot (bool): whether to use pivoting 
 * @return L, Lt (arma::mats): decomposed L and Ltranspose 
 */
std::pair<arma::Mat<double>, arma::Mat<double>> pivoted_cholesky(arma::Mat<double>& A, bool pivot) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::Mat<double> L(n,n);
    L.zeros();
    arma::Mat<double> P(n,n, arma::fill::eye);  
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

            // A.swap_cols(j, pivot_row); 
            // P.swap_cols(j, pivot_row);
     
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
                double diagonal = A(j,j) - sum;
      
                
                // if (diagonal <= 1e-10) {
                //     std::cout << "Warning. Numerical instability."<<std::endl;
                //     diagonal+= 1e-10; 

                // }
   
                L(j,j) = std::sqrt(std::abs(diagonal));
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    L = P*L; // apply permutations we kept track of in P
    arma::Mat<double> Lt = arma::trans(L); // U is upper traingular 
    arma::Mat<double> reconstructed = L*Lt;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        // std::cout << "Cholesky successful." << std::endl; 
        return {L,Lt};
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L, Lt};

}

/**
 * Find the pivot in a given submatrix 
 *
 *
 * @param A: matrix A (matrix to be decomposed)
 * @param start: starting index
 * @return pivot row and column indicies: look from (start, start) to bottom right (last element) of A for the pivot 
 */ 
std::pair<int,int> find_full_pivot(arma::Mat<double>&A, int start) {
    int n = A.n_rows;
    int pivot_row = start; 
    int pivot_col = start; 
    double max_val = std::abs(A(start, start));
    // find the max val 
    // compute elements below the idagonal rows below start and cols! 
    for (int i = start; i < n; i++) {
        for (int j = start; j < n; j++) {
            double val = std::abs(A(i,j));
            if (val > max_val) {
                max_val = val; 
                pivot_row = i; 
                pivot_col = j; 
                // row index of maximum val in a column

        }
        }
       
        
        
    }
    // std::cout << "max val: "<< max_val << std::endl;
    // std::cout << "pivot row: " << pivot_row << std::endl;
    return {pivot_row, pivot_col};


}


/**
 * Cholesky Decomposition with full pivoting 
 *
 *
 * @param A: matrix A (matrix to be decomposed)
 * @return L,Lt: matrix factors after decomposing A
 */ 
std::pair<arma::Mat<double>, arma::Mat<double>> full_pivoted_cholesky(arma::Mat<double>& A) {
    if (!A.is_symmetric() || !A.is_sympd()) {
        throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
    }
    int n = A.n_rows;
    arma::Mat<double> L(n,n);
    L.zeros();
    arma::Mat<double> P(n,n, arma::fill::eye);  
    int pivot_col; 
    int pivot_row; 

    // columns j 
    for (int j=0; j<n; j++) {
        // std::cout << "Before swapping:" << std::endl;
        // std::cout << "A:" << std::endl;
        // A.print();
        // std::cout << "P:" << std::endl;
        // P.print();
        std::tie(pivot_row, pivot_col) = find_full_pivot(A, j);

        if (pivot_row != j) {

            A.swap_rows(j, pivot_row); 
            P.swap_rows(j,pivot_row);

     
        }
        if (pivot_col != j) {
            A.swap_cols(j, pivot_col);
            P.swap_cols(j, pivot_col);
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
                double diagonal = A(j,j) - sum;
                // if (diagonal <= 1e-10) {
                //     std::cout << "Warning. Numerical instability."<<std::endl;
                //     diagonal+= 1e-10; 

                // }
   
                L(j,j) = std::sqrt(std::abs(diagonal));
            } 
            else {
                L(j,i) = (A(j,i) - sum)/L(i,i);
            }
        }
    }
    L = P*L; // apply permutations we kept track of in P
    arma::Mat<double> Lt = arma::trans(L); // U is upper traingular 
    arma::Mat<double> reconstructed = L*Lt;
    if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
        // std::cout << "Cholesky successful." << std::endl; 
        return {L,Lt};
        
    }
    else {
        std::cout << "Cholesky failed." << std::endl;
    }
    
    return {L, Lt};

}


/**
 * LU Decomposition of the matrix A
 *
 *
 * @param A: matrix A (matrix to be decomposed)
 * @param pivot (bool): whether to use pivoting 
 * @return L,U: factors 
 */ 
std::pair<arma::Mat<double>, arma::Mat<double>> LU_decomp(arma::Mat<double> & A, bool pivot) {
    int n = A.n_rows;
    arma::Mat<double> L(n,n,arma::fill::eye); 
    arma::Mat<double> U = A; 
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

                    //////////////// END anything below doesn't work/isn't used/////////////////

// std::pair<arma::Mat<double>, arma::Mat<double>> other_chol(arma::Mat<double>& A, bool pivot) {
//     if (!A.is_symmetric() || !A.is_sympd()) {
//         throw std::invalid_argument("Input matrix is not symmetric/positive definite."); 
//     }
//     int n = A.n_rows;
//     arma::Mat<double> L(n,n);
//     L.zeros();
//     arma::Mat<double> P(n,n, arma::fill::eye); // row swaps 
//     for (int j = 0; j < n; ++j) {
//         double pivot = A(j,j);
//         int max_idx = j;
//         for (int i = j + 1; i < n; ++i) {
//             if (A(i,j) > A(max_idx, j)) {
//             // if (std::abs(A(i,j)) > std::abs(A(max_idx, j))) {
//                 max_idx = i;
//             }
//         }
//         if (max_idx != j) {
//             A.swap_rows(j, max_idx);
//             L.swap_rows(j, max_idx);
//             P.swap_rows(j, max_idx);
//         }

//         // Update pivot after potential swap
//         pivot = A(j, j);


//         L(j,j) = std::sqrt(std::abs(pivot));
//         for (int i = j+1; i<n; i++) {
//             L(i,j) =A(i,j)/L(j,j);
//             for (int k = j+1; k <= i; k++) {
//                 A(i, k) -= L(i,j)*L(k,j);
//             }
//         }
//     }
//     // L = P*L;
//     arma::Mat<double> Lt = arma::trans(L); // U upper triangular
//     arma::Mat<double> reconstructed = L*Lt;
//     if (arma::approx_equal(A, reconstructed, "absdiff", 1e-4)) {
//         return {L,Lt};
        
//     }
//     else {
//         std::cout << "Cholesky failed." << std::endl;
//     }
    
//     return {L,Lt};

// }




