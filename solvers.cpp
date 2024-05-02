#include <iostream>
#include <stdexcept> 
#include "cholesky.h"
#include "matgen.h"
#include <fstream>
void chol_testing(int n, bool pivot) {
    arma::mat A = gen_sympd(n);

    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 


}
int main(void) {
    int n = 10; 
    arma::mat A = gen_sympd(n); 

    arma::vec b = arma::randu(n);
    std::cout << "Linear solve with Cholesky" << std::endl;
    arma::vec x = solve_lin(A, b, false);
    // double residual_norm = calc_diff(A, b, x);
    // std::cout << "residual norm: " << residual_norm << std::endl;
    std::vector<int> sizes = {1,5,10,25,50,100,250,500}; 
   

    arma::mat A2 = gen_sympd(n); 
    arma::vec b2 = arma::randu(n);
    std::cout << "Linear solve with pivoted Cholesky" << std::endl;
    arma::vec x2 = solve_lin(A2, b2, true);
    // double residual_norm2 = calc_diff(A2, b2, x2);
    // std::cout << "residual norm: " << residual_norm2 << std::endl;
    return 0; 
}