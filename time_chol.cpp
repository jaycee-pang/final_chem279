#include <iostream> 
#include "matgen.h"
#include "cholesky.h"
#include <armadillo> 
#include <chrono> 
#include <fstream> 
#include "testing.h"




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


    std::ofstream outfile4("LU_times.txt"); 
    if (!outfile4.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        double time = LU_timing(size, false); 
        outfile4 << size << "\t" << time << std::endl;
        
    }
    outfile4.close();

    std::ofstream outfile5("LU_pivot_times.txt"); 
    if (!outfile5.is_open()) {
        std::cerr << "File error" << std::endl;
        return 1;
    }
    for (int size:mat_sizes) {
        double time = LU_timing(size, true); 
        outfile5 << size << "\t" << time << std::endl;
        
    }
    outfile5.close();


    
    return 0; 
}