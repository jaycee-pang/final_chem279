#include <iostream>
#include <stdexcept> 
#include "cholesky.h"
#include "matgen.h"
#include "testing.h"
#include <fstream>


int main(void) {
    int n = 100; 
    
    std::cout << "Linear solve with Cholesky" << std::endl;
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::vector<arma::Mat<double>> matrices;
    std::vector<arma::vec> bvecs; 
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
     
        bvecs.push_back(arma::randu(size));
    }
    std::ofstream outfile1("cholesky_Axb_err.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> A_arma = A; 
        arma::vec b=bvecs[i];
        arma::vec b_arma = b; 
        arma::vec x = solve_lin(A, b, false);

        arma::vec residual = A*x-b;
        arma::vec x_arma = arma::solve(A_arma, b_arma);
        // double residual_norm = arma::norm(residual);
        // std::cout << "residual norm: " << residual_norm << std::endl;
        double error_norm = arma::norm(x-x_arma); 
       
        outfile1 << size << "\t" << error_norm << std::endl;
        

    }
    outfile1.close();


    std::cout << "\n\nLinear solve with Pivoted Cholesky" << std::endl;
    std::ofstream outfile2("cholesky_piv_Axb_err.txt");
    if (!outfile2.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> A_arma = A; 
        // arma::vec b = arma::randu(size);
        arma::vec b=bvecs[i];
        arma::vec b_arma = b; 
        arma::vec x = solve_lin(A, b, true);

        arma::vec residual = A*x-b;
        arma::vec x_arma = arma::solve(A_arma, b_arma);
        // double residual_norm = arma::norm(residual);
        // std::cout << "residual norm: " << residual_norm << std::endl;
        double error_norm = arma::norm(x-x_arma); 
       
        outfile2 << size << "\t" << error_norm << std::endl;
        

    }
    outfile2.close();

    
    
    // for (int i=0; i<mat_sizes.size(); i++) {
    //     int size = mat_sizes[i];
    //     arma::Mat<double> A = matrices[i]; 
    //     arma::Mat<double> A_arma = A; 
    //     arma::vec b=bvecs[i];
    //     arma::vec b_arma = b; 
    
    //     arma::vec x = solve_lin(A, b, true);
    //     arma::vec x_arma = arma::solve(A_arma, b_arma);
    //     // arma::vec residual = A*x-b;
    //     // double residual_norm = arma::norm(residual);
    //     // std::cout << "residual norm: " << residual_norm << std::endl;
    //     double error_norm = arma::norm(x-x_arma); 
    //     std::cout << "error norm: " << error_norm << std::endl;
        

    // }
    
    


    std::cout << "\n\ntesting inverse matrix solving for Chol no pivoting" << std::endl;
    std::ofstream outfile3("cholesky_inv_err.txt");
    if (!outfile3.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        double error_inv = inverse_testing(A, false);
        std::cout << "Error inverse from Chol no pivot: " << error_inv <<std::endl;
        outfile3 << size << "\t" <<error_inv << std::endl;
    }
    outfile3.close(); 




    std::cout << "\n\ntesting inverse matrix solving for Chol with pivoting" << std::endl;
    std::ofstream outfile4("piv_cholesky_inv_err.txt");
    if (!outfile4.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        double error_inv = inverse_testing(A, true);
        std::cout << "Error inverse from Chol pivot: " << error_inv <<std::endl;
        outfile4 << size << "\t" <<error_inv << std::endl;
    }
    outfile4.close(); 



    return 0; 
}