#include "cholesky.h"
#include "matgen.h"
#include "testing.h"
#include <iostream> 
int main(void) {
    arma::Mat<double> A = gen_sympd(1000);
    std::vector<int> mat_sizes = {100,200,500,1000,1500,2000,2500}; 
    std::vector<arma::Mat<double>> matrices;
    std::vector<arma::Mat<double>> sing_matrices; 
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
        sing_matrices.push_back(gen_singular(size));
    }
    std::cout << "Pivoted Cholesky on sympd" << std::endl;
    for (int i=0; i <mat_sizes.size(); i++) {
        arma::Mat<double> A = matrices[i]; 
        arma::Mat<double> og_A = A; 
        std::pair<arma::Mat<double> , arma::Mat<double>> result = pivoted_cholesky(A,true);
        arma::Mat<double> L = result.first; 
        arma::Mat<double> Lt = result.second;
        arma::Mat<double> reconA = L*Lt; 
        double norm = arma::norm(og_A - reconA, "fro");
        std::cout << "norm: " << norm << std::endl;
    }
    std::cout << "Cholesky on sympd" << std::endl;
    for (int i=0; i <mat_sizes.size(); i++) {
        arma::Mat<double> A2 = matrices[i]; 
        arma::Mat<double> og_A2 = A2; 
        std::pair<arma::Mat<double> , arma::Mat<double>> result2 = pivoted_cholesky(A2,false);
        arma::Mat<double> L2 = result2.first; 
        arma::Mat<double> Lt2 = result2.second;
        arma::Mat<double> reconA2 = L2*Lt2; 
        double norm2 = arma::norm(og_A2 - reconA2, "fro");
        std::cout << "norm: " << norm2 << std::endl;
    }
    std::cout << "Cholesky with nearly singular matrices" << std::endl;
    for (int i=0; i <mat_sizes.size(); i++) {
        arma::Mat<double> A3 = matrices[i]; 
        arma::Mat<double> og_A3 = A3; 
        std::pair<arma::Mat<double> , arma::Mat<double>> result3 = pivoted_cholesky(A3,false);
        arma::Mat<double> L3 = result3.first; 
        arma::Mat<double> Lt3 = result3.second;
        arma::Mat<double> reconA3 = L3*Lt3; 
        double norm3 = arma::norm(og_A3 - reconA3, "fro");
        std::cout << "norm: " << norm3 << std::endl;
    }
    std::cout << "Pivoted Cholesky with nearly singular matrices" << std::endl;
    for (int i=0; i <mat_sizes.size(); i++) {
        arma::Mat<double> A4 = matrices[i]; 
        arma::Mat<double> og_A4 = A4; 
        std::pair<arma::Mat<double> , arma::Mat<double>> result4 = pivoted_cholesky(A4,true);
        arma::Mat<double> L4 = result4.first; 
        arma::Mat<double> Lt4 = result4.second;
        arma::Mat<double> reconA4 = L4*Lt4; 
        double norm4 = arma::norm(og_A4 - reconA4, "fro");
        std::cout << "norm: " << norm4 << std::endl;
    }

    

    
    

    

    return 0; 

}