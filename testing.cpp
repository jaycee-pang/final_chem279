#include "testing.h"
/**
 * Test inverse matrix construction using Cholesky decomposition. 
 *
 * Decompose the matrix A, compute inverse A as (L.t) inverse * L inverse 
 * Compare the Frobenius norm of the difference between arma's A inverse and 
 * our reconstructed solution
 *
 * @param A: matrix A (matrix to be decomposed)
 * @param pivot: boolean flag, use pivoting (true) or do not use pivoting (false)
 * @return L,Lt: matrix factors after decomposing A
 */ 
double inverse_testing(arma::Mat<double> & A, bool pivot) {
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::Mat<double> Ainv_arma = arma::inv(A);
    arma::Mat<double> Linv = arma::inv(L);
    arma::Mat<double> Ltinv = arma::inv(Lt);
    arma::Mat<double> Ainv_chol = Ltinv * Linv;
    arma::Mat<double> diff = Ainv_arma - Ainv_chol;
    

    double diff_norm = arma::norm(diff, "fro");
    
    return diff_norm;

}

/**
 * Find the reconstruction error from our matrix decomposition. 
 *
 * @param A: matrix A (matrix to be decomposed)
 * @param pivot: boolean flag, use pivoting (true) or do not use pivoting (false)
 * @return error: Frobenius norm of the difference between the original matrix and the 
 *                  reconstructed matrix. 
 */ 
double chol_err(arma::Mat<double> & A, bool pivot) {
    arma::Mat<double> og_A = A; 
    // A.print("original");
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::Mat<double> reconA = L*Lt;
    arma::Mat<double> diff =og_A -  reconA; 
    // reconA.print("reconstructed A"); 
    double error = arma::norm(diff, "fro"); 
    return error; 

}

/**
 * Find the reconstruction error from LU decomposition. 
 *
 * @param A: matrix A (matrix to be decomposed)
 * @return error: Frobenius norm of the difference between the original matrix and the 
 *                  reconstructed matrix. 
 */ 
double LU_decomp_err(arma::Mat<double> & A, bool pivot) {
    arma::Mat<double> originalA = A; 
    
    std::pair<arma::Mat<double>, arma::Mat<double>> result = LU_decomp(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::Mat<double> reconstructedA = L*Lt;
    arma::Mat<double> diff = originalA  - reconstructedA; 
    
    
    double error = arma::norm(diff, "fro"); 
    return error; 

}

/*
// solve Ly =b storing y in x 
for (int i=1; i <=n ; i++) {
    for (sum=b[i]; k=i-1; k>=1; k--) 
        sum -= A(i,k)*x[k];
        x[i] = sum/p[i];s

}
for (i=n; i >=1; i--) {
    // solve L.t*x=y
    for (sum=x[i], k=i+1; k<=n; k++) {
        sum -= A(k,i) * x[k];
        x[i]= sum/p[i]
    }
}
*/

/**
 * Forward substitution. Ly=b
 *
 * @param L: lower triangular matrix 
 * @param b: vec 
 * @return y from Lx=y; 
 */ 
arma::vec forward_sub(arma::Mat<double>&L, arma::vec&b) {
    arma::vec y(L.n_rows, arma::fill::zeros);
    for (int i = 0; i <L.n_rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j<i; j++) {
            sum += L(i,j)*y(j);
        }
        y(i) = (b(i) - sum)/L(i,i);
    }
    return y; 
    
}

/**
 * Forward substitution. L.t*x = y;
 *
 * @param L: lower triangular matrix 
 * @param y: from forward sub 
 * @return x: vector solution to the system 
 */ 
arma::vec backward_sub(arma::Mat<double>&L,  arma::vec&y) {
    arma::vec x(L.n_rows, arma::fill::zeros);
    for (int i = L.n_rows-1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i+1; j<L.n_rows; j++) {
            sum += L(j,i) * x(j);
        }

        x(i) = (y(i)-sum)/L(i, i);
    }

    return x;
}

/**
 * Solve system of linear equations Ax=b. 
 * Using Cholesky decomposition of A. 
 *
 * @param A: matrix to be decomposed 
 * @param b: vector
 * @return x: vector solution 
 */ 
arma::vec solve_lin(arma::Mat<double> & A,arma::vec& b, bool pivot) {
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::vec y = forward_sub(L, b); 
    arma::vec x = backward_sub(Lt, y);

    return x;  

}


/**
 * Time Cholesky decomposition
 *
 * @param A: matrix to be decomposed 
 * @param pivot: boolean flag for pivoting 
 * @reutrns: duration, time taken for the function 
 */ 
double chol_timing(arma::Mat<double> & A, bool pivot) {
    
    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();

}


/**
 * Time LU decomposition
 *
 * @param A: matrix to be decomposed 
 * @param pivot: boolean flag for pivoting 
 * @reutrns: duration, time taken for the function 
 */ 
double LU_timing(arma::Mat<double> & A, bool pivot) {
    // arma::Mat<double> A = gen_sympd(n);

    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::Mat<double>, arma::Mat<double>> result = LU_decomp(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();

}

/**
 * Time arma decomposition
 *
 * @param A: matrix to be decomposed 
 * @reutrns: duration, time taken for the function 
 */ 
double arma_timing(arma::Mat<double> & A) {

    arma::Mat<double> L; 
    auto start = std::chrono::high_resolution_clock::now(); 
    bool success = arma::chol(L,A);
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();


}



// void chol_testing(arma::Mat<double> & A, bool pivot) {
//     std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 


// }


// double calc_diff(arma::Mat<double>& A,  arma::vec& b, arma::vec& x) {
//     arma::vec residual = A*x-b; 
//     double r_norm = arma::norm(residual); 
//     return r_norm; 


// }




// double inverse_testing(arma::Mat<double> & A, bool pivot) {
//     std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
//     arma::Mat<double> L = result.first; 
//     arma::Mat<double> Lt = result.second; 
//     arma::Mat<double> Linv = arma::inv(L); 
//     arma::Mat<double> Ltinv = arma::inv(Lt); 
//     arma::Mat<double> Ainv = Ltinv * Linv; 
//     arma::Mat<double> I(A.n_rows, A.n_cols, arma::fill::eye);
//     arma::Mat<double> A_i = A*Ainv; 
//     arma::Mat<double> recon3 = L*Lt; 
//     double tolerance = 1e-6;
//     if (arma::approx_equal(A_i,I, "absdiff", tolerance)) {
//         std::cout << "Inverse of A matches." << std::endl;
//     } else {
//         std::cout << "Inverse of A not found." << std::endl;
//     }

//     arma::Mat<double> error = A_i - I; 
//     double norm_err = arma::norm(error, "fro"); 
//     return norm_err; 
    

// }