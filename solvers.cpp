#include <iostream>
#include <stdexcept> 
#include "cholesky.h"
#include "matgen.h"
#include "testing.h"
#include <fstream>
/*
Solve systems of linear equations and find matrix inverse using Cholesky and pivoted Cholesky. 
Error is evaluated as the norm of the difference between 
arma::solve for x and my Cholesky decomposition of matrix A into L and solving for x using forward and 
backward substituttion. 

Solve for inverse of A using Ltinv and Linv, then evaluate error between arma::inverse and 
A^-1 using (Lt inverse) *(L inverse). 
   
*/


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
        double error_norm = arma::norm(x_arma-x); 
       
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
        double error_norm = arma::norm(x_arma-x); 
       
        outfile2 << size << "\t" << error_norm << std::endl;
        

    }
    outfile2.close();


    std::cout << "\n\nLinear solve with Full Pivoted Cholesky" << std::endl;
    std::ofstream outfile3a("cholesky_full_piv_Axb_err.txt");
    if (!outfile3a.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> A_arma = A; 
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A); 
        arma::Mat<double> L = result.first; 
        arma::Mat<double> Lt = result.second; 
        arma::vec b=bvecs[i];
        arma::vec b_arma = b; 
        arma::vec y = forward_sub(L, b); 
        arma::vec x = backward_sub(Lt, y); 

        arma::vec residual = A*x-b;
        arma::vec x_arma = arma::solve(A_arma, b_arma);
        double error_norm = arma::norm(x_arma-x); 
       
        outfile3a << size << "\t" << error_norm << std::endl;
        

    }
    outfile3a.close();

    
    


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


    std::cout << "\n\ntesting inverse matrix solving for Chol with full pivoting" << std::endl;
    std::ofstream outfile4a("full_piv_cholesky_inv_err.txt");
    if (!outfile4a.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A); 
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
        double error_inv = arma::norm(error, "fro"); 
        std::cout << "Error inverse from Chol pivot: " << error_inv <<std::endl;
        outfile4a << size << "\t" <<error_inv << std::endl;
    }
    outfile4a.close(); 




    return 0; 
}