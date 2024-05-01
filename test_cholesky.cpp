#include <iostream> 
#include "cholesky.h"
#include "matgen.h"
#include <exception> 
bool is_posdef(const arma::mat& A) {
    int n = A.n_rows;
    for (int i = 1; i <= n; ++i) {
        arma::mat minor = A.submat(0, 0, i - 1, i - 1);
        if (arma::det(minor) <= 0) {
            return false;
        }
    }
    return true;
}

int main(void) {
    // arma::mat A = gen_sympd(5);
    arma::mat A = {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    }; 
    std::cout << is_posdef(A) << std::endl;
    std::cout << A.is_sympd() << std::endl;

    try {
        std::pair<arma::mat , arma::mat> result = cholesky(A);
        // std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, true);  
        std::cout <<std::endl;
        A.print("A");
        std::cout << "\n L:"<< std::endl;
        result.first.print();  
        std::cout << "\n L trnaspose:"<< std::endl; // upper triangular
        result.second.print(); 
        arma::mat reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
       

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    std::cout << "\nPivoted Cholesky" << std::endl;
    A = {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    }; 
    try {
        A.print("A");
        // A = P.t*L*L.t*P
        std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, true);  
        std::cout <<std::endl;
        
        std::cout << "\n Permutation matrix:"<< std::endl;
        result.first.print(); 
        std::cout << "\n L:"<< std::endl;
        result.second.print(); 
        arma::mat reconstructedA = result.first.t() * result.second * result.second.t() * result.first;
        reconstructedA.print("reconstructed"); 

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }


    std::cout << "\n new Pivoted Cholesky" << std::endl;
    A = {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    }; 
    try {
        A.print("A");

        std::pair<arma::umat, arma::mat> result = full_pivoted_cholesky(A, true);  
        std::cout <<std::endl;
        
        std::cout << "\n Permutation matrix:"<< std::endl;
        result.first.print(); 
        std::cout << "P size: " << result.first.size() << std::endl;

        std::cout << "\n L:"<< std::endl;
        result.second.print(); 
        std::cout << "L size: " << result.second.size() << std::endl;
        
        arma::mat reconstructedA = result.first.t() * result.second * result.second.t() * result.first;
        reconstructedA.print("reconstructed"); 

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    

    return 0; 
}