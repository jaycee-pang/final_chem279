#include "testing.h"

double inverse_testing(arma::Mat<double> & A, bool pivot) {
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::Mat<double> Linv = arma::inv(L); 
    arma::Mat<double> Ltinv = arma::inv(Lt); 
    arma::Mat<double> Ainv = Ltinv * Linv; 
    arma::Mat<double> I(A.n_rows, A.n_cols, arma::fill::eye);
    arma::Mat<double> A_i = A*Ainv; 
    arma::Mat<double> recon3 = L*Lt; 
    double tolerance = 1e-6;
    if (arma::approx_equal(A_i,I, "absdiff", tolerance)) {
        std::cout << "Inverse of A matches." << std::endl;
    } else {
        std::cout << "Inverse of A not found." << std::endl;
    }

    arma::Mat<double> error = A_i - I; 
    double norm_err = arma::norm(error, "fro"); 
    return norm_err; 

}


double chol_err(arma::Mat<double> & A, bool pivot) {
    arma::Mat<double> og_A = A; 
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    arma::Mat<double> reconA = L*Lt;

    arma::Mat<double> diff =og_A -  reconA; 
    double error = arma::norm(diff, "fro"); 
    return error; 

}
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
// Ax=b
arma::vec solve_lin(arma::Mat<double> & A,arma::vec& b, bool pivot) {
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    arma::Mat<double> L = result.first; 
    arma::Mat<double> Lt = result.second; 
    // std::cout << "L rows L cols: " << L.n_rows <<","<< L.n_cols << std::endl;
    // std::cout << "Ltrows cols: " <<Lt.n_rows <<","<< Lt.n_cols << std::endl;
    
    arma::vec y = forward_sub(L, b); 
    arma::vec x = backward_sub(Lt, y);

    
    return x;  

}

double calc_diff(arma::Mat<double>& A,  arma::vec& b, arma::vec& x) {
    arma::vec residual = A*x-b; 
    double r_norm = arma::norm(residual); 
    return r_norm; 


}


void chol_testing(arma::Mat<double> & A, bool pivot) {
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 


}


double chol_timing(arma::Mat<double> & A, bool pivot) {
    
    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();

}


double LU_timing(arma::Mat<double> & A, bool pivot) {
    // arma::Mat<double> A = gen_sympd(n);

    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::Mat<double>, arma::Mat<double>> result = LU_decomp(A, pivot); 
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();

}

double arma_timing(arma::Mat<double> & A) {

    arma::Mat<double> L; 
    auto start = std::chrono::high_resolution_clock::now(); 
    bool success = arma::chol(L,A);
    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    return duration.count();


}