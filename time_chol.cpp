// #include <iostream> 
// #include "matgen.h"
// #include "cholesky.h"
// #include <armadillo> 
// #include <chrono> 
// #include <fstream> 
// #include "testing.h"


// int main(void) {
//     std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
//     std::ofstream outfile1("cholesky_times.txt");
//     if (!outfile1.is_open()) {
//         std::cerr <<"File error"<<std::endl;
//         return 1; 
//     }
//     for (int size:mat_sizes) {
//         double time = chol_timing(size, false); 
//         outfile1 << size << "\t" << time << std::endl;
        

//     }
//     outfile1.close();

//     std::ofstream outfile2("pivoted_cholesky_times.txt"); 
//     if (!outfile2.is_open()) {
//         std::cerr << "File error" << std::endl;
//         return 1;
//     }
//     for (int size:mat_sizes) {
//         double time = chol_timing(size, true); 
//         outfile2 << size << "\t" << time << std::endl;
        
//     }
//     outfile2.close();

//     std::ofstream outfile3("armachol_times.txt"); 
//     if (!outfile3.is_open()) {
//         std::cerr << "File error" << std::endl;
//         return 1;
//     }
//     for (int size:mat_sizes) {
//         double time = arma_timing(size); 
//         outfile3 << size << "\t" << time << std::endl;
        
//     }
//     outfile3.close();


//     std::ofstream outfile4("LU_times.txt"); 
//     if (!outfile4.is_open()) {
//         std::cerr << "File error" << std::endl;
//         return 1;
//     }
//     for (int size:mat_sizes) {
//         double time = LU_timing(size, false); 
//         outfile4 << size << "\t" << time << std::endl;
        
//     }
//     outfile4.close();

//     std::ofstream outfile5("LU_pivot_times.txt"); 
//     if (!outfile5.is_open()) {
//         std::cerr << "File error" << std::endl;
//         return 1;
//     }
//     for (int size:mat_sizes) {
//         double time = LU_timing(size, true); 
//         outfile5 << size << "\t" << time << std::endl;
        
//     }
//     outfile5.close();


    
//     return 0; 
// }

#include <iostream> 
#include "matgen.h"
#include "cholesky.h"
#include <armadillo> 
#include <chrono> 
#include <fstream> 
#include "testing.h"


int main(void) {
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::vector<arma::Mat<double>> matrices;
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
    }
    std::ofstream outfile1("cholesky_times.txt");
    if (!outfile1.is_open()) {
        std::cerr <<"File error"<<std::endl;
        return 1; 
    }
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


    std::ofstream outfile6("pivoted_cholesky_times_sing.txt"); 
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
        auto start = std::chrono::steady_clock::now();
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A, true);
        auto end = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
        double time = duration.count();
        outfile6 << size << "\t" << time << std::endl;
        
        
    }
    outfile6.close();

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


    
    return 0; 
}