#include <iostream>
#include <stdexcept> 
#include "cholesky.h"
#include "matgen.h"
#include "testing.h"
#include <fstream>


int main(void) {
    int n = 100; 
    arma::Mat<double> A = gen_sympd(n); 

    arma::vec b = arma::randu(n);
    std::cout << "Linear solve with Cholesky" << std::endl;
    arma::vec x = solve_lin(A, b, false);
    arma::vec residual = A*x-b;
    double residual_norm = arma::norm(residual);
    std::cout << "residual norm: " << residual_norm << std::endl;
    
    arma::Mat<double> A2 = gen_sympd(n); 
    arma::vec b2 = arma::randu(n);
    std::cout << "Linear solve with pivoted Cholesky" << std::endl;
    arma::vec x2 = solve_lin(A2, b2, true);
    arma::vec residual2 = A2*x2-b2;
    double residual_norm2 = arma::norm(residual2);
    std::cout << "residual norm: " << residual_norm2 << std::endl;


    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::cout << "testing inverse matrix solving for Chol no pivoting" << std::endl;
    std::ofstream outfile1("cholesky_inv_err.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int size:mat_sizes) {
        arma::Mat<double> A = gen_sympd(n); 
        double error_inv = inverse_testing(A, false);
        std::cout << "Error inverse from Chol no pivot: " << error_inv <<std::endl;
        outfile1 << size << "\t" <<error_inv << std::endl;
    }
    outfile1.close(); 

    std::cout << "testing inverse matrix solving for Chol with pivoting" << std::endl;
    std::ofstream outfile2("piv_cholesky_inv_err.txt");
    if (!outfile2.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int size:mat_sizes) {
        arma::Mat<double> A = gen_sympd(n); 
        double error_inv = inverse_testing(A, true);
        std::cout << "Error inverse from Chol pivot: " << error_inv <<std::endl;
        outfile2 << size << "\t" <<error_inv << std::endl;
    }
    outfile2.close(); 



    return 0; 
}