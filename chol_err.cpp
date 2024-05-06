#include <iostream>
#include <fstream> 
#include "cholesky.h"
#include "matgen.h"
#include <armadillo> 
#include "testing.h"
int main(void) {
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::vector<arma::Mat<double>> matrices;
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
    }
    std::ofstream outfile1("cholesky_error.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        double error = chol_err(A, false); 
        outfile1 << size << "\t" << error << std::endl;
        
        
    }
    outfile1.close();

    std::ofstream outfile2("pivoted_cholesky_error.txt"); 
    if (!outfile2.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        double error = chol_err(A, true); 
        outfile2 << size << "\t" << error << std::endl;
        
    }
    outfile2.close();

    std::ofstream outfile3("arma_error.txt"); 
    if (!outfile3.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        arma::Mat<double> L = arma::chol(A); 
        arma::Mat<double> reconstructedA = L*L.t(); 
        arma::Mat<double> diff = A - reconstructedA; 
        double error = arma::norm(diff, "fro"); 
        outfile3 << size << "\t" << error << std::endl;
        
    }
    outfile3.close();

    std::ofstream outfile4("LU_decomp_error.txt"); 
    if (!outfile4.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        double error = LU_decomp_err(A, false); 
        outfile4 << size << "\t" << error << std::endl;
        
    }
    outfile4.close();

    std::ofstream outfile5("LU_pivot_decomp_error.txt"); 
    if (!outfile5.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        double error = LU_decomp_err(A, true); 
        outfile5 << size << "\t" << error << std::endl;
        
    }
    outfile5.close();

    std::ofstream outfile6("pivoted_cholesky_err_sing.txt"); 
    if (!outfile6.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    std::vector<arma::Mat<double>> almost_singular;
    for (int size : mat_sizes) {
        almost_singular.push_back(gen_singular(size));
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = almost_singular[i]; 
        double error = chol_err(A, true); 
        outfile6 << size << "\t" << error << std::endl;
        
        
    }
    outfile6.close();

    std::ofstream outfile7("full_piv_cholesky_error.txt");
    if (!outfile7.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> og_A = A; 
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A); 
        // double error = chol_err(A, false); 
        arma::Mat<double> L = result.first; 
        arma::Mat<double> Lt = result.second; 
        arma::Mat<double> reconA = L*Lt;

        arma::Mat<double> diff =og_A -  reconA; 
        double error = arma::norm(diff, "fro"); 
        outfile7 << size << "\t" << error << std::endl;
        
        
    }
    outfile7.close();


    std::ofstream outfile8("full_pivoted_cholesky_err_sing.txt"); 
    if (!outfile8.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
   
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = almost_singular[i]; 
        arma::Mat<double> og_A = A; 
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A); 
        arma::Mat<double> L = result.first; 
        arma::Mat<double> Lt = result.second; 
        arma::Mat<double> reconA = L*Lt;
        arma::Mat<double> diff =og_A -  reconA; 
        double error = arma::norm(diff, "fro"); 
 
        outfile8 << size << "\t" << error << std::endl;
        
        
    }
    outfile8.close();


 


    return 0; 

}