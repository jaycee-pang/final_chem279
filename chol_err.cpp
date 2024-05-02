#include <iostream>
#include <fstream> 
#include "cholesky.h"
#include "matgen.h"
#include <armadillo> 
#include "testing.h"
int main(void) {
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::ofstream outfile1("cholesky_error.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int size:mat_sizes) {
        arma::mat A = gen_sympd(size); 
        double error = chol_err(A, false); 
        outfile1 << size << "\t" << error << std::endl;
        
        
    }
    outfile1.close();

    std::ofstream outfile2("pivoted_cholesky_error.txt"); 
    if (!outfile2.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        arma::mat A = gen_sympd(size); 
        double error = chol_err(A, true); 
        outfile2 << size << "\t" << error << std::endl;
        
    }
    outfile2.close();

    std::ofstream outfile3("arma_error.txt"); 
    if (!outfile3.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        arma::mat A = gen_sympd(size); 
        // double error = chol_err(A, false); 
        arma::mat L = arma::chol(A); 
        arma::mat reconstructedA = L*L.t(); 
        arma::mat diff = A - reconstructedA; 
        double error = arma::norm(diff, "fro"); 
        outfile3 << size << "\t" << error << std::endl;
        
    }
    outfile3.close();


    return 0; 

}