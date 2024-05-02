#include "testing.h"




double inverse_testing(arma::mat & A, bool pivot) {
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
    arma::mat L = result.first; 
    arma::mat Lt = result.second; 
    arma::mat Linv = arma::inv(L); 
    arma::mat Ltinv = arma::inv(Lt); 
    arma::mat Ainv = Ltinv * Linv; 
    arma::mat I(A.n_rows, A.n_cols, arma::fill::eye);
    arma::mat A_i = A*Ainv; 
    arma::mat recon3 = L*Lt; 
    double tolerance = 1e-6;
    if (arma::approx_equal(A_i,I, "absdiff", tolerance)) {
        std::cout << "Inverse of A matches." << std::endl;
    } else {
        std::cout << "Inverse of A not found." << std::endl;
    }

    arma::mat error = A_i - I; 
    double norm_err = arma::norm(error, "fro"); 
    return norm_err; 

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
arma::vec forward_sub( arma::mat&L, arma::vec&b) {
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

arma::vec backward_sub(arma::mat&L,  arma::vec&y) {
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
// Ax=b
arma::vec solve_lin(arma::mat & A, arma::vec& b, bool pivot) {
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
    arma::mat L = result.first; 
    arma::mat Lt = result.second; 
    // std::cout << "L rows L cols: " << L.n_rows <<","<< L.n_cols << std::endl;
    // std::cout << "Ltrows cols: " <<Lt.n_rows <<","<< Lt.n_cols << std::endl;
    
    arma::vec y = forward_sub(L, b); 
    arma::vec x = backward_sub(Lt, y);

    
    return x;  

}

double calc_diff(arma::mat& A,  arma::vec& b, arma::vec& x) {
    arma::vec residual = A*x-b; 
    double r_norm = arma::norm(residual); 
    return r_norm; 


}


void chol_testing(int n, bool pivot) {
    arma::mat A = gen_sympd(n);

    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 


}


double chol_timing(int n, bool pivot) {
    arma::mat A = gen_sympd(n);
    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
    return duration.count();

}


double LU_timing(int n, bool pivot) {
    arma::mat A = gen_sympd(n);
    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::mat, arma::mat> result = LU_decomp(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
    return duration.count();

}

double arma_timing(int n) {
    arma::mat A = gen_sympd(n); 
    arma::mat L; 
    auto start = std::chrono::high_resolution_clock::now(); 
    bool success = arma::chol(L,A);
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
    return duration.count();


}