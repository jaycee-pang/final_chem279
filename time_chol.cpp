#include <iostream> 
#include "matgen.h"
#include "cholesky.h"
#include <armadillo> 
#include <chrono> 
#include <fstream> 
double chol_timing(int n, bool pivot) {
    arma::mat A = gen_sympd(n);
    auto start = std::chrono::high_resolution_clock::now(); 
    std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A, pivot); 
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


int main(void) {
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::ofstream outfile1("cholesky_times.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    for (int size:mat_sizes) {
        double time = chol_timing(size, false); 
        outfile1 << size << "\t" << time << std::endl;
        

    }
    outfile1.close();

    std::ofstream outfile2("pivoted_cholesky_times.txt"); 
    if (!outfile2.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        double time = chol_timing(size, true); 
        outfile2 << size << "\t" << time << std::endl;
        
    }
    outfile2.close();

    std::ofstream outfile3("armachol_times.txt"); 
    if (!outfile3.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        double time = arma_timing(size); 
        outfile3 << size << "\t" << time << std::endl;
        
    }
    outfile3.close();


    
    return 0; 
}