#include <iostream> 
#include "matgen.h"
#include "cholesky.h"
#include <armadillo> 
#include <chrono> 
#include <fstream> 
#include "testing.h"
/*
Testing runtimes for matrix decomposition for increasing matrix sizes.
*/

int main(void) {
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::vector<arma::Mat<double>> matrices;
    // so we can reuse the same matrix for every test
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
    }
    std::vector<arma::Mat<double>> almost_singular;
    for (int size : mat_sizes) {
        almost_singular.push_back(gen_singular(size));
    }
    std::ofstream outfile1("cholesky_times.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
    std::cout << "Cholesky"<<std::endl;
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, false);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile1 << size << "\t" << time << std::endl;
        

    }
    outfile1.close();
    std::cout << "\nPivoted Cholesky"<<std::endl;
    std::ofstream outfile2("pivoted_cholesky_times.txt"); 
    if (!outfile2.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, true);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile2 << size << "\t" << time << std::endl;
        
        
    }
    outfile2.close();
    std::cout << "\nArmadillo Cholesky"<<std::endl;
    std::ofstream outfile3("armachol_times.txt"); 
    if (!outfile3.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> L; 
        auto start = std::chrono::high_resolution_clock::now(); 
        bool success = arma::chol(L,A);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile3 << size << "\t" << time << std::endl;

            
    }
    outfile3.close();

    std::cout << "\nLU"<<std::endl;
    std::ofstream outfile4("LU_times.txt"); 
    if (!outfile4.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        auto start = std::chrono::high_resolution_clock::now(); 
        std::pair<arma::Mat<double>, arma::Mat<double>> result = LU_decomp(A,false); 
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile4 << size << "\t" << time << std::endl;
        
    }
    outfile4.close();
    std::cout << "\nPivoted LU"<<std::endl;
    std::ofstream outfile5("LU_pivot_times.txt"); 
    if (!outfile5.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        auto start = std::chrono::high_resolution_clock::now(); 
        std::pair<arma::Mat<double>, arma::Mat<double>> result = LU_decomp(A,true); 
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile5 << size << "\t" << time << std::endl;
        
    }
    outfile5.close();

    std::cout << "\nPivoted Cholesky on nearly singular"<<std::endl;
    std::ofstream outfile6("pivoted_cholesky_times_sing.txt"); 
    if (!outfile6.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = almost_singular[i]; 
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, true);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile6 << size << "\t" << time << std::endl;
        
        
    }
    outfile6.close();

    std::cout << "\nFull Pivoted Cholesky"<<std::endl;
    std::ofstream outfile7("full_pivoted_cholesky_times.txt"); 
    if (!outfile7.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile7 << size << "\t" << time << std::endl;
        
        
    }
    outfile7.close();


    std::cout << "\nFull Pivoted Cholesky nearly singular"<<std::endl;
    std::ofstream outfile8("full_pivoted_cholesky_sing_times.txt"); 
    if (!outfile8.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = almost_singular[i]; 
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = full_pivoted_cholesky(A);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile8 << size << "\t" << time << std::endl;
        
        
    }
    outfile8.close();


    
    return 0; 
}